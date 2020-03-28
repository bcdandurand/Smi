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

		// generate OSI solver object
		// 	here we use OsiClp
#if 0
		OsiCpxSolverInterface *clp = new OsiCpxSolverInterface();
        smi.solveWS(clp, 1.0);


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
        OsiCpxSolverInterface *osi_cpx = reinterpret_cast<OsiCpxSolverInterface*>(osiStoch);
        CPXsetintparam( osi_cpx->getEnvironmentPtr(), CPXPARAM_Threads, 0);

		// solve
		osiStoch->initialSolve();
        osiStoch->branchAndBound();

		// print results
		printf("Solved stochastic program %s\n", name);
		printf("Number of rows: %d\n",osiStoch->getNumRows());
		printf("Number of cols: %d\n",osiStoch->getNumCols());
		printf("Optimal value: %g\n",osiStoch->getObjValue());

		// print solution to file
		string outfilename(name);
		const string suffix(".out");
		outfilename = outfilename + suffix;
		FILE *fp = fopen(outfilename.c_str(),"w");
#endif
        //OsiCpxSolverInterface *osi_cpx = new OsiCpxSolverInterface();
        //double result = smi.getWSValue(osi_cpx, 1.0 );
        //cout << "Result: " << result << endl;
		int numScenarios=smi.getNumScenarios();

        OsiSolverInterface **osi_cpx_scns = new OsiSolverInterface*[numScenarios];
		for (int ii=0 ; ii<numScenarios; ++ii) {
#if 0
			double *dsoln=NULL;
			int numCols=0;
			fprintf(fp,"Scenario %d \n",ii);
			dsoln = smi.getColSolution(ii,&numCols);
			for (int jj=0; jj<numCols; jj++)
				fprintf(fp,"%g \n",dsoln[jj]);
			free(dsoln);
#endif
            cout << "Subproblem: " << ii << endl;
            osi_cpx_scns[ii] = new OsiCpxSolverInterface();
            osi_cpx_scns[ii] = smi.loadOsiSolverDataForScenarioSP(osi_cpx_scns[ii], ii);
            //osi_cpx_scns[ii] = dynamic_cast<OsiCpxSolverInterface*>(osi_cpx_scns[ii]);
#if 0
        for (int mi = 0; mi < osi_cpx_scns[ii]->getMatrixByRow()->getMinorDim();mi++){
            std::cout << osi_cpx_scns[ii]->getObjCoefficients()[mi] << "*x" << mi << " + ";
        }
        std::cout << std::endl;
        printf("\n Objective printed \n");
        for (int mi = 0; mi < osi_cpx_scns[ii]->getMatrixByRow()->getMajorDim(); mi++){
            printf("\n%g <= ", osi_cpx_scns[ii]->getRowLower()[mi]);
            for (int mj = 0; mj < osi_cpx_scns[ii]->getMatrixByRow()->getMinorDim(); mj++){
                printf("%g ",osi_cpx_scns[ii]->getMatrixByRow()->getCoefficient(mi,mj));
            }
            printf("<= %g",osi_cpx_scns[ii]->getRowUpper()[mi]);
        }
        printf("\n Matrix printed \n");
        for (int mi = 0; mi < osi_cpx_scns[ii]->getNumCols(); mi++) {
            printf("\n%g <= x%d <= %g",osi_cpx_scns[ii]->getColLower()[mi],mi,osi_cpx_scns[ii]->getColUpper()[mi]);
        }
        printf("\n");
        printf("\n Var bounds printed \n");
#endif
            // Set objSense
            osi_cpx_scns[ii]->setObjSense(1.0);
            //char fname[128];
            //sprintf(fname,"LpFiles/lp%d",ii);
            //cout << fname << endl;
            //osi_cpx_scns[ii]->writeLp(fname);
            //CPXsetintparam( osi_cpx_scns[ii]->getEnvironmentPtr(), CPXPARAM_Threads, 0);
            //CPXsetintparam( osi_cpx_scns[ii]->getEnvironmentPtr(), CPXPARAM_ScreenOutput, 0);

		    osi_cpx_scns[ii]->initialSolve();
            osi_cpx_scns[ii]->branchAndBound();
		    printf("Scenario %d optimal value: %g\n",ii,osi_cpx_scns[ii]->getObjValue());

            //delete osi_cpx_scns[ii];
            cout << "End of subproblem: " << ii << endl;
		}
		//fclose(fp);
        //delete [] osi_cpx_scns;

}

// Display message on stdout and stderr
void testingMessage( const char * const msg )
{
//  std::cerr <<msg;
  cout <<endl <<"*****************************************"
       <<endl <<msg <<endl;
}
