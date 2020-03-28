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

		int numScenarios=smi.getNumScenarios();

        OsiCpxSolverInterface **osi_cpx_scns = new OsiCpxSolverInterface*[numScenarios];
		for (int ii=0 ; ii<numScenarios; ++ii){
            cout << "Subproblem: " << ii << endl;
            osi_cpx_scns[ii] = new OsiCpxSolverInterface();
            smi.loadOsiSolverDataForScenarioSP(osi_cpx_scns[ii], ii);
            osi_cpx_scns[ii]->setObjSense(1.0);  // Set objective sense, MIN=1.0
	        osi_cpx_scns[ii]->messageHandler()->setLogLevel(0); //Nice hack to suppress all output

            int status=0;
            int nThreads=1;  //Might want to experiment with this if you're using the more recent versions of CPLEX
            int MIP_TOL=1e-6;
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

		    osi_cpx_scns[ii]->initialSolve();
            osi_cpx_scns[ii]->branchAndBound();
		    printf("Scenario %d optimal value: %g\n",ii,osi_cpx_scns[ii]->getObjValue());

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
