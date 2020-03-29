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
void loadOsiSolverDataForScenarioSP(SmiScnModel smi, OsiSolverInterface *osi, int scn);
int setQuadDiagCoeffs(OsiCpxSolverInterface *osi_cpx, double *diag_coeffs);
int solveAsQMIP(OsiCpxSolverInterface *osi_cpx);

int main()
{

	testingMessage( "Model generation using SMPS files for SIPLIB problems.\n" );
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
        std::vector<int> intIndices = smi.getIntIndices();
        cout << "Integer indices are: " << endl;
    //int* getIntegerInd(){ return core_->getIntegerIndices();}
    //int getIntegerLen(){return core_->getIntegerLength(); }
    //int* getBinaryInd(){return core_->getBinaryIndices();}
    //int getBinaryLen(){return core_->getBinaryLength();}
        cout << endl;
	    int first_var_idx = (smi.getCore())->getColStart(0);
	    int last_var_idx = (smi.getCore())->getNumCols(0);  
        int total_ncols = (smi.getCore())->getNumCols(0) + (smi.getCore())->getNumCols(1);   
        double *diag_els = new double[total_ncols];
        for (int kk=0; kk< total_ncols; kk++){
            if (kk < last_var_idx){
                diag_els[kk]=1.0;
            }
            else{
                diag_els[kk]=0.0;
            }
        }
        int status=0;
        int nThreads=1;  //Might want to experiment with this if you're using the more recent versions of CPLEX
        double MIP_TOL=1e-6;
        int mode = 0; // mode = 0 for MIP, mode = 1 for QMIP
        int opt_target;
        //int opt_target=CPX_OPTIMALITYTARGET_OPTIMALCONVEX; //For problems known a priori to be convex
        if (mode == 1){
            //opt_target=CPX_OPTIMALITYTARGET_OPTIMALGLOBAL;  //Use for nonconvex, e.g., QMIP
        }
		int numScenarios=smi.getNumScenarios();
        std::vector<OsiCpxSolverInterface*> osi_cpx_scns;

#if 0
        SmiDiscreteDistribution *smiDD = new SmiDiscreteDistribution(smi.getCore());
        smi.processDiscreteDistributionIntoScenarios(smiDD, false);
        int n_scens = 0;
	    int num_rv = smiDD->getNumRV(); 
        cout << "num_rv: " << num_rv << endl;
        for (int jj=0; jj<num_rv; jj++){
	        SmiDiscreteRV * discRV = smiDD->getDiscreteRV(jj); 
            if(discRV->getStage()==1){  // second-stage
	            int n_events = discRV->getNumEvents();
                cout << "rv " << jj <<": n_events=" << n_events << endl;
                for(int ee=0; ee<n_events; ee++){
                    osi_cpx_scns.push_back(new OsiCpxSolverInterface());
	                const CoinPackedMatrix &matrix = discRV->getEventMatrix(ee);
	                const CoinPackedVector &dclo = discRV->getEventColLower(ee);
	                const CoinPackedVector &dcup = discRV->getEventColUpper(ee); 
	                const CoinPackedVector &dobj = discRV->getEventObjective(ee);
	                const CoinPackedVector &drlo = discRV->getEventRowLower(ee); 
	                const CoinPackedVector &drup = discRV->getEventRowUpper(ee); 
	                double prob = discRV->getEventProb(ee); 
	                // pass data to osiStoch
                    osi_cpx_scns[n_scens]->loadProblem(matrix,NULL,NULL,NULL,NULL,NULL);
                    for(int ii=0; ii < dclo.getNumElements(); ii++){
                        osi_cpx_scns[n_scens]->setColLower((dclo.getIndices())[ii], (dclo.getElements())[ii]);
                    }
                    for(int ii=0; ii < dcup.getNumElements(); ii++){
                        osi_cpx_scns[n_scens]->setColUpper((dcup.getIndices())[ii], (dcup.getElements())[ii]);
                    }
                    for(int ii=0; ii < drlo.getNumElements(); ii++){
                        osi_cpx_scns[n_scens]->setRowLower((drlo.getIndices())[ii], (drlo.getElements())[ii]);
                    }
                    for(int ii=0; ii < drup.getNumElements(); ii++){
                        osi_cpx_scns[n_scens]->setRowUpper((drup.getIndices())[ii], (drup.getElements())[ii]);
                    }
                    for(int ii=0; ii < dobj.getNumElements(); ii++){
                        osi_cpx_scns[n_scens]->setObjCoeff((dobj.getIndices())[ii], (dobj.getElements())[ii]);
                    }
                    // load integer values in solver
                    for (unsigned int ii = 0; ii < intIndices.size(); ii++) {
                        osi_cpx_scns[n_scens]->setInteger(intIndices[ii]);
                    }
                    osi_cpx_scns[n_scens]->setObjSense(1.0);  // Set objective sense, MIN=1.0
                    n_scens++;
                }
            }
        }

        cout << "Number of subproblems generated is: " << n_scens << " while numScenarios is: " << numScenarios << endl;

#endif
        
		for (int ii=0 ; ii<numScenarios; ++ii){
            cout << "Subproblem: " << ii << endl;
            osi_cpx_scns.push_back(new OsiCpxSolverInterface());
            smi.loadOsiSolverDataForScenarioSP(osi_cpx_scns[ii], ii);
            osi_cpx_scns[ii]->setObjSense(1.0);  // Set objective sense, MIN=1.0
            osi_cpx_scns[ii]->switchToMIP();
            for (unsigned int nn = 0; nn < smi.getIntegerLen(); nn++) {
                osi_cpx_scns[ii]->setInteger((smi.getIntegerInd())[nn]);
            }
	        osi_cpx_scns[ii]->messageHandler()->setLogLevel(0); //Nice hack to suppress all output

            status=CPXsetintparam( osi_cpx_scns[ii]->getEnvironmentPtr(), CPXPARAM_Threads, 1);
	        status=CPXsetintparam( osi_cpx_scns[ii]->getEnvironmentPtr(), CPXPARAM_Parallel, CPX_PARALLEL_DETERMINISTIC); //no iteration message until solution
	        status=CPXsetintparam( osi_cpx_scns[ii]->getEnvironmentPtr(), CPXPARAM_Threads, nThreads); //no iteration message until solution
            // Tolerance relevant parameters
            #if 0
	        status=CPXsetdblparam( osi_cpx_scns[ii]->getEnvironmentPtr(), CPXPARAM_MIP_Tolerances_AbsMIPGap, MIP_TOL);
	        status=CPXsetdblparam( osi_cpx_scns[ii]->getEnvironmentPtr(), CPXPARAM_MIP_Tolerances_MIPGap, MIP_TOL*1e-3);
	        status=CPXsetdblparam( osi_cpx_scns[ii]->getEnvironmentPtr(), CPXPARAM_MIP_Tolerances_Integrality, MIP_TOL);
	        status=CPXsetdblparam( osi_cpx_scns[ii]->getEnvironmentPtr(), CPXPARAM_Simplex_Tolerances_Optimality, MIP_TOL);
	        status=CPXsetdblparam( osi_cpx_scns[ii]->getEnvironmentPtr(), CPXPARAM_Simplex_Tolerances_Feasibility, MIP_TOL);
            #endif
            #if 0
            status=CPXsetintparam( osi_cpx_scns[ii]->getEnvironmentPtr(), CPXPARAM_Simplex_Display, 0); //no iteration message until solution
            status=CPXsetintparam( osi_cpx_scns[ii]->getEnvironmentPtr(), CPXPARAM_Tune_Display, 0); //turn off display
            status=CPXsetintparam( osi_cpx_scns[ii]->getEnvironmentPtr(), CPX_PARAM_SCRIND, CPX_OFF); //turn off display
            status=CPXsetintparam( osi_cpx_scns[ii]->getEnvironmentPtr(), CPXPARAM_Barrier_Display, 0); //turn off display
            #endif
            #if 0
            #endif
            if (mode == 0){
		        osi_cpx_scns[ii]->initialSolve();
                osi_cpx_scns[ii]->branchAndBound();
            }
            else if (mode==1){
		        //osi_cpx_scns[ii]->initialSolve();
	            //status=CPXsetintparam( osi_cpx_scns[ii]->getEnvironmentPtr(), CPXPARAM_OptimalityTarget, opt_target);
                setQuadDiagCoeffs(osi_cpx_scns[ii], diag_els);
                solveAsQMIP(osi_cpx_scns[ii]);
            }
            

		    printf("Scenario %d optimal value: %g\n",ii,osi_cpx_scns[ii]->getObjValue());
            cout << "Solution is: " << endl;
            for(int nn=0; nn<osi_cpx_scns[ii]->getNumCols(); nn++){
                cout << " " << (osi_cpx_scns[ii]->getColSolution())[nn];
            } 
            cout << endl;

            //delete osi_cpx_scns[ii];
            cout << "End of subproblem: " << ii << endl;
		}
        //delete [] osi_cpx_scns;

}

// Display message on stdout and stderr
void testingMessage( const char * const msg )
{
//  std::cerr <<msg;
  cout <<endl <<"*****************************************"
       <<endl <<msg <<endl;
}
void loadOsiSolverDataForScenarioSP(SmiScnModel smi, OsiSolverInterface *osi, int scn){

}

int setQuadDiagCoeffs(OsiCpxSolverInterface *osi_cpx, double *diag_coeffs){
    int err;
    err = CPXchgprobtype(osi_cpx->getEnvironmentPtr(), osi_cpx->getLpPtr(), CPXPROB_MIQP );
    if ( err ) {
        cerr << "Failed to set problem type to MIQP." << endl;
    }
    err = CPXcopyqpsep(osi_cpx->getEnvironmentPtr(), osi_cpx->getLpPtr(OsiCpxSolverInterface::KEEPCACHED_PROBLEM), diag_coeffs);
    if ( err ) {
        cerr << "Failed to set diagonal quadratic coefficients." << endl;
    }
    return err;
}

//Precondition: setQuadDiagCoeffs() has been successfully called
int solveAsQMIP(OsiCpxSolverInterface *osi_cpx){
    int status = CPXmipopt(osi_cpx->getEnvironmentPtr(), osi_cpx->getLpPtr( OsiCpxSolverInterface::FREECACHED_RESULTS ));
    if ( status ) {
        cerr << "Failed to optimize MIQP." << endl;
    }
    return status;
}
