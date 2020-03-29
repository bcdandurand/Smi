// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.

#include <string>
using namespace std;

#include <cassert>
#include <iostream>

#include "CoinPragma.hpp"
#include "SmiScnModel.hpp"
#include "SmiScnData.hpp"
#include "OsiClpSolverInterface.hpp"
#include "OsiCpxSolverInterface.hpp"
#include "cplex.h"

//#define DATASTOCHASTICDIR "/homes/bcdandurand/COIN-OR/Data/Stochastic"
#define DATASTOCHASTICDIR "/homes/bcdandurand/Data/SIPLIB"

//forward declarations

void testingMessage(const char * const);
void SmpsIO(const char* const);
vector<vector<double> > getColSolutionsByStageForScn(SmiScnModel &smi, int ns);
vector<vector<int> > extractScnSPColsByStage(SmiScnModel &smi, int ns);
vector<vector<int> > extractScnSPRowsByStage(SmiScnModel &smi, int ns);
vector<double> extractScnProbTrace(SmiScnModel &smi, int ns);
vector<int> extractScnSubmatrix(OsiSolverInterface *smiOsi, CoinPackedMatrix &spMat, vector<vector<int> > &stg_cols, vector<vector<int> > &stg_rows);
vector<int> extractScnSP(SmiScnModel &smi, OsiSolverInterface *smiOsi, int ns, OsiSolverInterface *osi);

int main()
{

	testingMessage( "Model generation using SMPS files for Cambridge-Watson problems.\n" );
	SmpsIO(DATASTOCHASTICDIR"/SSLP/sslp_5_25_50");

	testingMessage( "*** Done! *** \n");

 	return 0;
}

void SmpsIO(const char * const name )
{
		SmiScnModel smi;

		// read SMPS model from files
		//	<name>.core, <name>.time, and <name>.stoch
		smi.readSmps(name);

		// generate OSI solver object
		// 	here we use OsiClp
		OsiCpxSolverInterface *clp = new OsiCpxSolverInterface();

		// set solver object for SmiScnModel
		smi.setOsiSolverHandle(*clp);

		// load solver data
		// 	this step generates the deterministic equivalent
		//	and returns an OsiSolver object
		OsiSolverInterface *osiStoch = smi.loadOsiSolverData();

		// set some nice Hints to the OSI solver
		osiStoch->setHintParam(OsiDoPresolveInInitial,true);
		osiStoch->setHintParam(OsiDoScale,true);
		osiStoch->setHintParam(OsiDoCrash,true);

		// solve
#if 1
		osiStoch->initialSolve();
        osiStoch->branchAndBound();

		// print results
		printf("Solved stochastic program %s\n", name);
		printf("Number of rows: %d\n",osiStoch->getNumRows());
		printf("Number of cols: %d\n",osiStoch->getNumCols());
		printf("Optimal value: %g\n",osiStoch->getObjValue());
#endif

		// print solution to file
		int numScenarios=smi.getNumScenarios();
#if 0
		for (int i=0 ; i<numScenarios; ++i) {
            vector< vector<double> > solnsByStage = getColSolutionsByStageForScn(smi, i);
            cout << "Scenario " << i << endl;
            for(size_t stg=0; stg<solnsByStage.size(); stg++){
                cout << "\tStage " << stg + 1 << " soln: ";
                for(size_t jj=0; jj<solnsByStage[stg].size(); jj++){
                    cout << "  " << solnsByStage[stg][jj];
                }
                cout << endl;
            }
		}
#endif

        double optval_sum=0.0;
        std::vector<OsiCpxSolverInterface*> osi_cpx_scns;
		for (int ii=0 ; ii<numScenarios; ++ii) {
            osi_cpx_scns.push_back(new OsiCpxSolverInterface());
            cout << "Extracting subproblem: " << ii << endl;
            vector<int> int_indx = extractScnSP(smi, osiStoch, ii, osi_cpx_scns[ii]);
            cout << "Done extracting subproblem: " << ii << endl;
            osi_cpx_scns[ii]->switchToMIP();
            //for (unsigned int nn = 0; nn < smi.getIntegerLen(); nn++) {
            for (unsigned int nn = 0; nn < int_indx.size(); nn++) {
                //osi_cpx_scns[ii]->setInteger((smi.getIntegerInd())[nn]);
                osi_cpx_scns[ii]->setInteger(int_indx[nn]);
            }
            if( ii==0 ){
                osi_cpx_scns[ii]->writeLp("extract_model.out");
            }
	        osi_cpx_scns[ii]->messageHandler()->setLogLevel(0); //Nice hack to suppress all output
		    osi_cpx_scns[ii]->initialSolve();
            osi_cpx_scns[ii]->branchAndBound();
		    printf("Scenario %d optimal value: %g\n",ii,osi_cpx_scns[ii]->getObjValue());
            cout << "Solution is: " << endl;
            for(int nn=0; nn<osi_cpx_scns[ii]->getNumCols(); nn++){
                cout << " " << (osi_cpx_scns[ii]->getColSolution())[nn];
            } 
            cout << endl;
		    optval_sum += osi_cpx_scns[ii]->getObjValue();
        }
        cout << "Average optimal value: " << optval_sum << endl;


}


// Display message on stdout and stderr
void testingMessage( const char * const msg )
{
//  std::cerr <<msg;
  cout <<endl <<"*****************************************"
       <<endl <<msg <<endl;
}

vector<vector<double> > getColSolutionsByStageForScn(SmiScnModel &smi, int ns){
    const double * osiSoln = (smi.getOsiSolverInterface())->getColSolution();
    int numcols=0;
    vector<vector<double> > solnsByStage;
    SmiScnNode *node = smi.getLeafNode(ns);
    while (node != NULL){
        solnsByStage.push_back(vector<double>());
        node = node->getParent();
    }


    size_t n_stages = solnsByStage.size();
    size_t stg = n_stages;
    node = smi.getLeafNode(ns);
    while (node != NULL){
        stg--;
        // copy entries
        // getColStart returns the starting index of node in OSI model
        for(int j=node->getColStart(); j<node->getColStart()+node->getNumCols(); ++j){
            // getCoreColIndex returns the corresponding Core index
            // in the original (user's) ordering
            solnsByStage[stg].push_back(osiSoln[j]);
        }
        // get parent of node
        node = node->getParent();
    }
    return solnsByStage;
}

vector<vector<int> > extractScnSPColsByStage(SmiScnModel &smi, int ns){
    vector<vector<int> > stg_cols;
    size_t n_stg = 0;
    SmiScnNode *node = smi.getLeafNode(ns);
    while (node != NULL){
        stg_cols.push_back(vector<int>());
        node = node->getParent();
        n_stg++;
    }
    size_t stg = n_stg;
    node = smi.getLeafNode(ns);
    while (node != NULL){
        stg--;
        for(int jj=node->getColStart(); jj<node->getColStart()+node->getNumCols(); ++jj){
            stg_cols[stg].push_back(jj); 
        }
        node = node->getParent();
    }
    return stg_cols;
}
vector<vector<int> > extractScnSPRowsByStage(SmiScnModel &smi, int ns){
    vector<vector<int> > stg_rows;
    size_t n_stg = 0;
    SmiScnNode *node = smi.getLeafNode(ns);
    while (node != NULL){
        stg_rows.push_back(vector<int>());
        node = node->getParent();
        n_stg++;
    }
    size_t stg = n_stg;
    node = smi.getLeafNode(ns);
    while (node != NULL){
        stg--;
        for(int jj=node->getRowStart(); jj<node->getRowStart()+node->getNumRows(); ++jj){
            stg_rows[stg].push_back(jj); 
        }
        node = node->getParent();
    }
    return stg_rows;
}
vector<double> extractScnProbTrace(SmiScnModel &smi, int ns){
    vector<double> stg_probs;
    SmiScnNode *node = smi.getLeafNode(ns);
    while (node != NULL){
        stg_probs.push_back(node->getProb());
        node = node->getParent();
    }
    return stg_probs;
}
vector<int> extractScnSubmatrix(OsiSolverInterface *smiOsi, CoinPackedMatrix &spMat, vector<vector<int> > &stg_cols, vector<vector<int> > &stg_rows){
    const CoinPackedMatrix *smiMat = smiOsi->getMatrixByCol();
    size_t n_stg = stg_cols.size();
    size_t n_scn_cols=0;
    size_t n_scn_rows=0;
    for(size_t stg=0; stg<n_stg; stg++){
        n_scn_cols += stg_cols[stg].size();
        n_scn_rows += stg_rows[stg].size();
    }
    int * indMajor = new int[n_scn_cols];
    int col_idx=0;
    for( size_t stg=0; stg<n_stg; stg++){
        for( size_t col=0; col<stg_cols[stg].size(); col++){
            indMajor[col_idx++] = stg_cols[stg][col];
        }
    }
    spMat.submatrixOf(*smiMat, n_scn_cols, indMajor);
    delete [] indMajor;
    vector<int> rows_to_delete;
#if 0
    int row_del_start=0;
    cout << "Deleting the rows: ";
    for( size_t stg=0; stg<n_stg; stg++){
        for( size_t row=row_del_start; row<stg_rows[stg][0]; row++){
            rows_to_delete.push_back(row);
            cout << " " << row;
        }
        cout << endl;
        size_t n_rows = stg_rows[stg].size();
        row_del_start=stg_rows[stg][n_rows-1]+1;
    }
#endif
    return rows_to_delete;
}
//Precondition: concrete class constructor called for osi.
vector<int> extractScnSP(SmiScnModel &smi, OsiSolverInterface *smiOsi, int ns, OsiSolverInterface *osi){
    CoinPackedMatrix spMat;

    vector<vector<int> > stg_cols = extractScnSPColsByStage(smi, ns);
    vector<vector<int> > stg_rows = extractScnSPRowsByStage(smi, ns);
    vector<int> int_indx;

    vector<int> rows_to_delete = extractScnSubmatrix(smiOsi, spMat, stg_cols, stg_rows);
    vector<double> scn_prob_trace = extractScnProbTrace(smi, ns);
    cout << "Scen "<< ns << " probs: ";
    for(size_t kk=0; kk<scn_prob_trace.size(); kk++){
        cout << " " << scn_prob_trace[kk];
    }
    cout << endl;
    
	// pass data to osiStoch
    osi->loadProblem(spMat,NULL,NULL,NULL,NULL,NULL);
    cout << "Num cols: " << osi->getNumCols() << endl;
    cout << "Num rows: " << osi->getNumRows() << endl;
    int sp_col_idx=0;
    cout << "Scn column indices: " << endl;
    for( size_t stg=0; stg<stg_cols.size(); stg++){
        for( size_t col=0; col<stg_cols[stg].size(); col++){
            int col_idx=stg_cols[stg][col];
            double prob = scn_prob_trace[stg];
            if(smiOsi->isInteger(col_idx)){
                int_indx.push_back(sp_col_idx);
            }
            cout << " " << col_idx;
            osi->setColLower(sp_col_idx,(smiOsi->getColLower())[col_idx]);
            osi->setColUpper(sp_col_idx,(smiOsi->getColUpper())[col_idx]);
            osi->setObjCoeff(sp_col_idx,prob*(smiOsi->getObjCoefficients())[col_idx]);
            sp_col_idx++;
        }
        cout << endl;
    }
    int sp_row_idx=0;
    cout << "Scn row indices: " << endl;
    for( size_t stg=0; stg<stg_rows.size(); stg++){
        for( size_t row=0; row<stg_rows[stg].size(); row++){
            int row_idx=stg_rows[stg][row];
            cout << " " << row_idx;
            osi->setRowLower(row_idx,(smiOsi->getRowLower())[row_idx]);
            osi->setRowUpper(row_idx,(smiOsi->getRowUpper())[row_idx]);
            sp_row_idx++;
        }
        cout << endl;
    }
    int n_rows_del = rows_to_delete.size();
    int *row_inds = new int[n_rows_del]; 
    for( size_t row_del=0; row_del<rows_to_delete.size(); row_del++){
        row_inds[row_del] = rows_to_delete[row_del];
    }
    osi->deleteRows(n_rows_del, row_inds);
    delete [] row_inds;
    // load integer values in solver
#if 0
    osi->switchToMIP();
    for (unsigned int nn = 0; nn < smi.getIntegerLen(); nn++) {
        osi->setInteger((smi.getIntegerInd())[nn]);
    }
#endif
    osi->setObjSense(smiOsi->getObjSense());  // Set objective sense, MIN=1.0
    return int_indx; 

}
