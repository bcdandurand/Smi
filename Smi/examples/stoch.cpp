// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <cassert>
#include <iostream>


#include "SmiScnModel.hpp"
#include "SmiScnData.hpp"
#include "OsiClpSolverInterface.hpp"

// Plus infinity
#ifndef COIN_DBL_MAX
#define COIN_DBL_MAX DBL_MAX
#endif
#define INF COIN_DBL_MAX
#define SMI_INFINITY COIN_DBL_MAX


//forward declarations

void testingMessage(const char * const);
void SmpsIO(const char* const);
void ModelScenario(const char* const);
void ModelDiscrete(const char* const);
void ModelBug();
void SmpsBug();

int main()
{


	testingMessage("Model generation for simple model Bug");
	ModelBug();

	testingMessage("Read SMPS version of simple model Bug");
	SmpsBug();

	testingMessage( "Model generation using scenario tree construction methods.\n");
	ModelScenario("Dantzig-Ferguson Aircraft Allocation using Scenarios");

	testingMessage( "Model generation using discrete distribution specification methods.\n" );
	ModelDiscrete("Dantzig-Ferguson Aircraft Allocation using Discrete Distribution");

	testingMessage( "Model generation using SMPS files for Watson problems.\n" );
	SmpsIO("../../Mps/Stochastic/wat_10_C_32");		

	testingMessage( "*** Done! *** \n");
	

  return 0;
}

void ModelBug()
{

/*
NAME    BUG
ROWS
  N  obj
  G  C0
  G  C1
  G  C2
  G  C3
*/
	int nCoreRows=4;
	double dCoreRup[] = { INF, INF, INF, INF };
/*
COLUMNS
   x01   obj   1
   x01   C3    1
   x01   C1    1
   x01   C0    1
   x02   obj   1
   x02   C2    1
   x02   C1    1
   x02   C0    1
   x03   obj   1
   x03   C3    1
   x03   C2    1
   x03   C0    1
   x04   obj   0.5
   x04   C3    1
   x04   C1    1
   x05   obj   0.5
   x05   C2    1
   x05   C1    1
   x06   obj   0.5
   x06   C3    1
   x06   C2    1
*/
	int nCoreCols=6;
	double dCoreClo[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	double dCoreCup[] = {INF, INF, INF, INF, INF, INF};
	double dCoreObj[] = {1.0, 1.0, 1.0, 0.5, 0.5, 0.5};

	int nCoreElts=15;
	int iCoreColStarts[] = {0,
						  3,
						  6,
						  9,
						  11,
						  13,
						  15};
	int iCoreRowIndice[] = {3,1,0,
						  2,1,0,
						  3,2,0,
						  3,1,
						  2,1,
						  3,2};
	double dCoreMatEntries[] = {1.0, 1.0, 1.0,
							  1.0, 1.0, 1.0,
							  1.0, 1.0, 1.0,
							  1.0, 1.0,
							  1.0, 1.0,
							  1.0, 1.0 };						

/*
RHS
  RHS    C0    0
  RHS    C1    1
  RHS    C2    1
  RHS    C3    1
ENDATA
*/
	double dCoreRlo[] = { 0.0, 1.0, 1.0, 1.0 };
/*
-----------------------------------------------------------
bug.time:

TIME          BUG
PERIODS       LP
     x01       C0                      STG01
     x04       C1                      STG02 
ENDATA 
-----------------------------------------------------------
*/
	int nCoreStages = 2;
	int iColStages[] = {0,0,0,1,1,1};
	int iRowStages[] = {0,1,1,1};

/*
----------------------------------------------------------
bug.stoch file:

NAME          BUG
SCENARIOS     DISCRETE                REPLACE
  SC SCEN01    ROOT           0.500    STG02
     RHS       C1             1.000
     RHS       C2             1.000
     RHS       C3             0.000
  SC SCEN02    ROOT           0.500    STG02
     RHS       C1             0.000
     RHS       C2             1.000
     RHS       C3             0.000 
ENDATA
*/

	int nScenarios = 2;
	int iBranchStage[] = {1,1};

	int iAncestorScn[] = {0,0};
	double dProbScn[] = {0.5, 0.5};
	
	double drlo0[] = {1.0, 1.0, 0.0 };
	int indices[]  = {1, 2, 3};
	CoinPackedVector *rlo0 = new CoinPackedVector(3, indices,drlo0);
	double drlo1[] = {0.0, 1.0, 0.0 };
	CoinPackedVector *rlo1 = new CoinPackedVector(3, indices,drlo1);

	// generate Core model
	
	OsiClpSolverInterface osi;	
	osi.loadProblem(nCoreCols, nCoreRows,iCoreColStarts,iCoreRowIndice,dCoreMatEntries,dCoreClo,dCoreCup,dCoreObj,
								dCoreRlo,dCoreRup);
	SmiCoreData *osiCore = new SmiCoreData(&osi,nCoreStages,iColStages,iRowStages);

	// initialize SmiScnModel
	SmiScnModel *smiModel = new SmiScnModel();

	// Add Scenarios
	int	is = smiModel->generateScenario(osiCore,NULL,NULL,NULL,NULL,
										rlo0,NULL,iBranchStage[0],iAncestorScn[0],dProbScn[0]);
	
	is = smiModel->generateScenario(osiCore,NULL,NULL,NULL,NULL,
										rlo1,NULL,iBranchStage[1],iAncestorScn[1],dProbScn[1]);

	// Set Solver
	OsiClpSolverInterface osiClp;
	smiModel->setOsiSolverHandle(osiClp);
	
	// Load Scenarios into Deterministic Equivalent
	smiModel->loadOsiSolverData();

	// Solve
	OsiSolverInterface *osiStoch = smiModel->getOsiSolverInterface();
	
	osiStoch->initialSolve();

			// print results
		printf("Solved stochastic program %s\n", "BUG");
		printf("Number of rows: %d\n",osiStoch->getNumRows());
		printf("Number of cols: %d\n",osiStoch->getNumCols());
		printf("Optimal value: %g\n",osiStoch->getObjValue());
		
		const double *drlo = osiStoch->getRowLower();
		const double *drup = osiStoch->getRowUpper();

		printf("Det Equiv rows\n");
		for (int i=0; i< osiStoch->getNumRows(); i++) 
			printf("Row  %d\t Lower %g \t Upper %g \n",i,drlo[i],drup[i]);
}

// MyCombineRule -- handle infinities in the correct way for SMPS model BUG.
class MyCombineRule: public SmiCoreCombineReplace
{
public:
	virtual void Process(double *d1, int o1, const CoinPackedVector &cpv2, char *type=0);
	virtual CoinPackedVector * Process(CoinPackedVector *cpv1, CoinPackedVector *cpv2, char *type=0);

}; 

void SmpsBug()
{
		SmiScnModel smi;

		char name[]="bug";

		MyCombineRule myCombineRule;

		// read SMPS model from files
		//	<name>.core, <name>.time, and <name>.stoch
		// the argument myCombineRule overrides the combine rule specified in the Stoch file
		smi.readSmps(name, &myCombineRule);		

		// generate OSI solver object
		// 	here we use OsiClp
		OsiClpSolverInterface *clp = new OsiClpSolverInterface();

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
		osiStoch->initialSolve();		

		// print results
		printf("Solved stochastic program %s\n", name);
		printf("Number of rows: %d\n",osiStoch->getNumRows());
		printf("Number of cols: %d\n",osiStoch->getNumCols());
		printf("Optimal value: %g\n",osiStoch->getObjValue());		
}	

CoinPackedVector *  MyCombineRule::Process(CoinPackedVector *cr,CoinPackedVector *nr, char *type)
{	

	CoinPackedVector *newrow=NULL;

	if (!cr && nr) newrow = new CoinPackedVector(*nr);
	if (cr && !nr) newrow = new CoinPackedVector(*cr);
	
	if (cr && nr)
	{
		int maxentries = CoinMax(cr->getMaxIndex(),nr->getMaxIndex());
		double* dense = cr->denseVector(maxentries+1);
		double* elt_nr = nr->getElements();
		int* ind_nr = nr->getIndices();
		
		for (int j=0; j<nr->getNumElements(); ++j)
			dense[ind_nr[j]] = elt_nr[j];
			
		// generate new packed vector
		newrow = new CoinPackedVector();
		for (int i=0; i<maxentries+1; ++i)
		{
			if (dense[i] && dense[i] < SMI_INFINITY )	// SPECIAL: handle infinite entry
				newrow->insert(i,dense[i]);
		}

		delete [] dense;
	}

	return newrow;

	
}

void MyCombineRule::Process(double *d, int o, const CoinPackedVector &cpv, char *type)
{
	const double *cd = cpv.getElements();
	const int *ci = cpv.getIndices();
	
	for (int j=0; j < cpv.getNumElements(); ++j) 
		if (cd[j] < SMI_INFINITY)						// SPECIAL: handle infinite entry
			d[ci[j]-o] = cd[j];

}





void SmpsIO(const char * const name )
{
		SmiScnModel smi;	

		// read SMPS model from files
		//	<name>.core, <name>.time, and <name>.stoch
		smi.readSmps(name);		

		// generate OSI solver object
		// 	here we use OsiClp
		OsiClpSolverInterface *clp = new OsiClpSolverInterface();

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
		osiStoch->initialSolve();		

		// print results
		printf("Solved stochastic program %s\n", name);
		printf("Number of rows: %d\n",osiStoch->getNumRows());
		printf("Number of cols: %d\n",osiStoch->getNumCols());
		printf("Optimal value: %g\n",osiStoch->getObjValue());		
}	

void ModelScenario(const char * const name )
{
	// example of direct interfaces for scenario generation
	
    /* Model dimensions */
    int ncol=27, nrow=9, nels=44;
	
	/* Sparse matrix data...organized by row */
    int *mrow,cmrow[]={ 0, 0, 0, 0, 0,
		1, 1, 1, 1,
		2, 2, 2,
		3, 3, 3, 3, 3,
		4, 4, 4, 4,
		5, 5, 5, 5, 5, 5,
		6, 6, 6, 6, 6,
		7, 7, 7, 7, 7, 7,
		8, 8, 8, 8, 8, 8 };
	  int *mcol,cmcol[]={ 0, 1, 2, 3, 4,
		5, 6, 7, 8,
		9,10, 11, 
		12, 13, 14, 15, 16, 
		0,        12, 17, 18,
		1, 5, 9,  13, 19, 20,
		2, 6,     14, 21, 22,
		3, 7, 10, 15, 23, 24,
		4, 8, 11, 16, 25, 26 };

    double dels[] = { 1.0, 1.0, 1.0, 1.0, 1.0,
		1.0, 1.0, 1.0, 1.0,
		1.0, 1.0, 1.0,
		1.0, 1.0, 1.0, 1.0, 1.0,
		16.0,              9.0, -1.0, 1.0,
		15.0, 10.0,  5.0, 11.0, -1.0, 1.0,
		28.0, 14.0,       22.0, -1.0, 1.0,
		23.0, 15.0,  7.0, 17.0, -1.0, 1.0,
		81.0, 57.0, 29.0, 55.0, -1.0, 1.0 };
	
    /* Objective */
    double *dobj,cdobj[]={ 18.0, 21.0, 18.0, 16.0, 10.0, 15.0, 16.0, 14.0, 9.0,
		10.0,  9.0,  6.0, 17.0, 16.0, 17.0, 15.0, 10.0, 0.0,
		13.0,  0.0, 13.0,  0.0,  7.0,  0.0,  7.0,  0.0, 1.0 };
	
    /* Column bounds */
    double *dclo,cdclo[]={ 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
		0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
		0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0 };



    double *dcup,cdcup[]={ INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF,
		INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF,
		INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF };
	
    /* Row bounds */
    double *drlo,cdrlo[]={ -INF, -INF, -INF, -INF,  0.0, 4.0, 0.0, 8.0, 10.0 };
    double *drup,cdrup[]={ 10.0, 19.0, 25.0, 15.0,  0.0, 7.0, 0.0, 8.0, 90.0 };
	
    /* Stages */
    //int nstg=2;
    int n_first_stg_rows=4;
    int *rstg,crstg[]={ 0,0,0,0,1,1,1,1,2 };
    int *cstg,ccstg[]={ 0,0,0,0,0,0,0,0,0,
		0,0,0,0,0,0,0,0,1,
		1,1,1,1,1,1,1,2,2 };
	
    /* Stochastic data */
    int nindp=5;
    int nsamp[]={ 5, 2, 5, 5, 3 };
    double demand[]={ 200, 220, 250, 270, 300,
		50, 150,
		140, 160, 180, 200, 220,
		10, 50, 80, 100, 340,
		580, 600, 620 };
    double dprobs[]={ 0.2, 0.05, 0.35, 0.2, 0.2,
		0.3, 0.7,
		0.1, 0.2, 0.4, 0.2, 0.1,
		0.2, 0.2, 0.3, 0.2, 0.1,
		0.1, 0.8, 0.1 };
	
	/* scramble */

	int irow[]={ 1,2,7,8,0,3,4,5,6};
    int icol[]={ 9,2,3,4,5,6,7,8,1,
		19,21,23,25,0,26,24,22,20,
		10,11,12,13,14,15,16,17,18 };

    /* local variables */
    int ns=1,ii,iii,jj,*indx,*incr;
    double dp=1.0;

    for (ii=0;ii<nindp;ii++) ns *= nsamp[ii];     /* Compute number of scenarios */

	// debug small sample
	// ns = 3;
	
	// initialize SmiModel
	SmiScnModel *smiModel = new SmiScnModel();


	OsiClpSolverInterface *osiClp1 = new OsiClpSolverInterface();

	smiModel->setOsiSolverHandle(*osiClp1);
	
	/* scramble LP entries */
	mrow = (int*)malloc(nels*sizeof(int));
	mcol = (int*)malloc(nels*sizeof(int));
	for (ii=0;ii<nels;ii++)
	{
		mcol[ii] = icol[cmcol[ii]];
		mrow[ii] = irow[cmrow[ii]];
	}
	
	drlo = (double *)malloc(nrow*sizeof(double));
	drup = (double *)malloc(nrow*sizeof(double));
	rstg = (int *)malloc(nrow*sizeof(int));
	for (ii=0;ii<nrow;ii++)
	{
		drlo[irow[ii]] = cdrlo[ii];
		drup[irow[ii]] = cdrup[ii];
		rstg[irow[ii]] = crstg[ii];
	}
	
	dclo = (double *)malloc(ncol*sizeof(double));
	dcup = (double *)malloc(ncol*sizeof(double));
	dobj = (double *)malloc(ncol*sizeof(double));
	cstg = (int *)malloc(ncol*sizeof(int));
	for (ii=0;ii<ncol;ii++)
	{
		dclo[icol[ii]] = cdclo[ii];
		dcup[icol[ii]] = cdcup[ii];
		dobj[icol[ii]] = cdobj[ii];
		cstg[icol[ii]] = ccstg[ii];
	}

	// this is a fake, just to demo that we can also handle stochastic matrix entries
	int corenels = nels - 4;

	// set core model using Osi interface
	OsiClpSolverInterface ocsi;
	ocsi.loadProblem(CoinPackedMatrix( 1,mrow,mcol,dels,corenels),dclo,dcup,dobj,drlo,drup);
	SmiCoreData *osiCore = new SmiCoreData(&ocsi,3,cstg,rstg);
	
	// Coin structures for scenario updates to right hand sides
	CoinPackedVector cpv_rlo;
	CoinPackedVector cpv_rup;

	// Coin structure for scenario "updates" to core matrix
	// ..row-ordered
	CoinPackedMatrix *cpm_mat = new CoinPackedMatrix(false,mrow+corenels,mcol+corenels,dels+corenels,nels-corenels);
		
    // initialize right hand side data for first scenario
    indx = (int *) malloc( (1+nindp)*sizeof(int) );
    memset( indx,0,(1+nindp)*sizeof(int));
    for (jj=0;jj<nindp;jj++) {
		indx[jj+1] += indx[jj] + nsamp[jj];
		dp *= dprobs[ indx[jj] ];

		drlo[irow[n_first_stg_rows + jj]] = demand[ indx[jj] ];
		drup[irow[n_first_stg_rows + jj]] = demand[ indx[jj] ];
		
		cpv_rlo.insert(irow[n_first_stg_rows + jj],demand[ indx[jj] ]);
		cpv_rup.insert(irow[n_first_stg_rows + jj],demand[ indx[jj] ]);
    }
	
	// first scenario
	int anc = 0;
	int branch = 1;
	int	is = smiModel->generateScenario(osiCore,cpm_mat,NULL,NULL,NULL,
									&cpv_rlo,&cpv_rup,branch,anc,dp);	
	

	
	/***** ...main loop to generate scenarios from discrete random variables
		For each scenario index ii:
        If the sample size nsamp[jj] divides the scenario index ii,
		reverse the increment direction incr[jj]
		and increase the random variable index jj by 1.
        Increment the jj'th random variable by incr[jj]
		and generate new sample data.
    ***** */
	
    /* sample space increment initialized to 1 */
    incr = (int *) malloc( nindp*sizeof(int) );
    for (jj=0;jj<nindp;jj++) incr[jj] = 1;
	
    for (int iss=1;iss<ns;iss++) {
		iii=iss; jj=0;
		while ( !(iii%nsamp[jj]) ) {
			iii /= nsamp[jj];
			incr[jj] = -incr[jj];
			jj++;
		}
		dp /= dprobs[ indx[jj] ];
		indx[jj] += incr[jj];
		dp *= dprobs[ indx[jj] ];

		// set data
		drlo[irow[n_first_stg_rows + jj]] = demand[ indx[jj] ];
		drup[irow[n_first_stg_rows + jj]] = demand[ indx[jj] ];

		cpv_rlo.setElement(cpv_rlo.findIndex(irow[n_first_stg_rows + jj]),demand[ indx[jj] ]);
		cpv_rup.setElement(cpv_rup.findIndex(irow[n_first_stg_rows + jj]),demand[ indx[jj] ]);
		
		// genScenario
		is = smiModel->generateScenario(osiCore,cpm_mat,NULL,NULL,NULL,
			&cpv_rlo,&cpv_rup,branch,anc,dp);	
		
		
	}
	
	assert(ns==smiModel->getNumScenarios());

	
	// load problem data into OsiSolver
	smiModel->loadOsiSolverData();
	// get Osi pointer
	OsiSolverInterface *smiOsi = smiModel->getOsiSolverInterface();
	// set some parameters
	smiOsi->setHintParam(OsiDoPresolveInInitial,true);
	smiOsi->setHintParam(OsiDoScale,true);
	smiOsi->setHintParam(OsiDoCrash,true);
	// solve using Osi Solver
	smiOsi->initialSolve();
	// test optimal value
    	assert(fabs(smiOsi->getObjValue()-1566.042)<0.01);

	// test solutions
	const double *dsoln = smiOsi->getColSolution();
	double objSum = 0.0;

	/* The canonical way to traverse the tree:
	   For each scenario, get the leaf node.
	   Then get the parent.  Repeat until parent is NULL.
	   (Only the root node has a NULL parent.)
	 */
	

	for(is=0; is<ns; ++is)
	{
		/* this loop calculates the scenario objective value */
		double scenSum = 0.0;

		// start with leaf node
		SmiScnNode *node = smiModel->getLeafNode(is);

		// leaf node probability is the scenario probability
		double scenprob = node->getModelProb();
	
		while (node != NULL)
		{

			// getColStart returns the starting index of node in OSI model
			for(int j=node->getColStart(); j<node->getColStart()+node->getNumCols(); ++j)
			{
				// getCoreColIndex returns the corresponding Core index
				// in the original (user's) ordering
				scenSum += dobj[node->getCoreColIndex(j)]*dsoln[j];	
				

			}			
			// get parent of node
			node = node->getParent();
		}
		objSum += scenSum*scenprob;
	}

	assert(fabs(smiOsi->getObjValue()-objSum) < 0.01);

		// print results
		printf("Solved stochastic program %s\n", name);
		printf("Number of rows: %d\n",smiOsi->getNumRows());
		printf("Number of cols: %d\n",smiOsi->getNumCols());
		printf("Optimal value: %g\n",smiOsi->getObjValue());

}

void ModelDiscrete(const char * const name)
{
	// example of direct interfaces for discrete distribution
	
    /* Model dimensions */
    int ncol=27, nrow=9, nels=44;
	
	/* Sparse matrix data...organized by row */
    int *mrow,cmrow[]={ 0, 0, 0, 0, 0,
		1, 1, 1, 1,
		2, 2, 2,
		3, 3, 3, 3, 3,
		4, 4, 4, 4,
		5, 5, 5, 5, 5, 5,
		6, 6, 6, 6, 6,
		7, 7, 7, 7, 7, 7,
		8, 8, 8, 8, 8, 8 };
	  int *mcol,cmcol[]={ 0, 1, 2, 3, 4,
		5, 6, 7, 8,
		9,10, 11, 
		12, 13, 14, 15, 16, 
		0,        12, 17, 18,
		1, 5, 9,  13, 19, 20,
		2, 6,     14, 21, 22,
		3, 7, 10, 15, 23, 24,
		4, 8, 11, 16, 25, 26 };

    double dels[] = { 1.0, 1.0, 1.0, 1.0, 1.0,
		1.0, 1.0, 1.0, 1.0,
		1.0, 1.0, 1.0,
		1.0, 1.0, 1.0, 1.0, 1.0,
		16.0,              9.0, -1.0, 1.0,
		15.0, 10.0,  5.0, 11.0, -1.0, 1.0,
		28.0, 14.0,       22.0, -1.0, 1.0,
		23.0, 15.0,  7.0, 17.0, -1.0, 1.0,
		81.0, 57.0, 29.0, 55.0, -1.0, 1.0 };
	
    /* Objective */
    /* Objective */
    double *dobj,cdobj[]={ 18.0, 21.0, 18.0, 16.0, 10.0, 15.0, 16.0, 14.0, 9.0,
		10.0,  9.0,  6.0, 17.0, 16.0, 17.0, 15.0, 10.0, 0.0,
		13.0,  0.0, 13.0,  0.0,  7.0,  0.0,  7.0,  0.0, 1.0 };
	
    /* Column bounds */
    double *dclo,cdclo[]={ 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
		0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
		0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0 };



    double *dcup,cdcup[]={ INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF,
		INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF,
		INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF,  INF };
	
    /* Row bounds */
    double *drlo,cdrlo[]={ -INF, -INF, -INF, -INF,  0.0, 4.0, 0.0, 8.0, 10.0 };
    double *drup,cdrup[]={ 10.0, 19.0, 25.0, 15.0,  0.0, 7.0, 0.0, 8.0, 90.0 };
	
    /* Stages */
    //int nstg=2;
    int n_first_stg_rows=4;
    int *rstg,crstg[]={ 0,0,0,0,1,1,1,1,2 };
    int *cstg,ccstg[]={ 0,0,0,0,0,0,0,0,0,
		0,0,0,0,0,0,0,0,1,
		1,1,1,1,1,1,1,2,2 };
	
    /* Stochastic data */
    int nindp=5;
    int nsamp[]={ 5, 2, 5, 5, 3 };
    double demand[]={ 200, 220, 250, 270, 300,
		50, 150,
		140, 160, 180, 200, 220,
		10, 50, 80, 100, 340,
		580, 600, 620 };
    double dprobs[]={ 0.2, 0.05, 0.35, 0.2, 0.2,
		0.3, 0.7,
		0.1, 0.2, 0.4, 0.2, 0.1,
		0.2, 0.2, 0.3, 0.2, 0.1,
		0.1, 0.8, 0.1 };
	
	/* scramble */

	int irow[]={ 1,2,7,8,0,3,4,5,6};
    int icol[]={ 9,2,3,4,5,6,7,8,1,
		19,21,23,25,0,26,24,22,20,
		10,11,12,13,14,15,16,17,18 };

    /* local variables */
    int ii,jj;


	// initialize SmiModel
	SmiScnModel *smiModel = new SmiScnModel();


	OsiClpSolverInterface *osiClp1 = new OsiClpSolverInterface();

	smiModel->setOsiSolverHandle(*osiClp1);
	
	/* scramble LP entries */
	mrow = (int*)malloc(nels*sizeof(int));
	mcol = (int*)malloc(nels*sizeof(int));
	for (ii=0;ii<nels;ii++)
	{
		mcol[ii] = icol[cmcol[ii]];
		mrow[ii] = irow[cmrow[ii]];
	}
	
	drlo = (double *)malloc(nrow*sizeof(double));
	drup = (double *)malloc(nrow*sizeof(double));
	rstg = (int *)malloc(nrow*sizeof(int));
	for (ii=0;ii<nrow;ii++)
	{
		drlo[irow[ii]] = cdrlo[ii];
		drup[irow[ii]] = cdrup[ii];
		rstg[irow[ii]] = crstg[ii];
	}
	
	dclo = (double *)malloc(ncol*sizeof(double));
	dcup = (double *)malloc(ncol*sizeof(double));
	dobj = (double *)malloc(ncol*sizeof(double));
	cstg = (int *)malloc(ncol*sizeof(int));
	for (ii=0;ii<ncol;ii++)
	{
		dclo[icol[ii]] = cdclo[ii];
		dcup[icol[ii]] = cdcup[ii];
		dobj[icol[ii]] = cdobj[ii];
		cstg[icol[ii]] = ccstg[ii];
	}

	// this to test the matrix update stanza in genScenario
	int corenels = nels - 4;
	
	
	// set core model using Osi interface
	OsiClpSolverInterface ocsi;
	ocsi.loadProblem(CoinPackedMatrix( 1,mrow,mcol,dels,corenels),dclo,dcup,dobj,drlo,drup);
	
	// core model with 3 stages
	SmiCoreData *smiCore = new SmiCoreData(&ocsi,3,cstg,rstg);

	// Coin structure for scenario "updates" to core matrix
	// ..row-ordered
	CoinPackedMatrix *cpm_mat = new CoinPackedMatrix(false,mrow+corenels,mcol+corenels,dels+corenels,nels-corenels);



	// Create discrete distribution
	SmiDiscreteDistribution *smiDD = new SmiDiscreteDistribution(smiCore);

	int index=0;
	for (jj=0;jj<nindp;jj++)
	{
		SmiDiscreteRV *smiRV = new SmiDiscreteRV(1);
		for (ii=0;ii<nsamp[jj];ii++)
		{			
			CoinPackedVector empty_vec;
			CoinPackedVector cpv_rlo ;
			CoinPackedVector cpv_rup ;
			cpv_rlo.insert(irow[n_first_stg_rows + jj], demand[index+ii]);
			cpv_rup.insert(irow[n_first_stg_rows + jj], demand[index+ii]);
			smiRV->addEvent(*cpm_mat,empty_vec,empty_vec,empty_vec,cpv_rlo,cpv_rup,dprobs[index+ii]);
			cpv_rlo.clear();
			cpv_rup.clear();
		}
		assert(smiRV->getNumEvents()==nsamp[jj]);
		for (ii=0;ii<nsamp[jj];ii++)
		{
			assert(smiRV->getEventColLower(ii).getNumElements()==0);
			assert(smiRV->getEventColUpper(ii).getNumElements()==0);
			assert(smiRV->getEventObjective(ii).getNumElements()==0);
			assert(smiRV->getEventMatrix(ii).getNumElements()==cpm_mat->getNumElements());
			assert(smiRV->getEventRowLower(ii).getElements()[0] == demand[index+ii]);
			assert(smiRV->getEventRowLower(ii).getIndices()[0] == irow[n_first_stg_rows + jj]);
			assert(smiRV->getEventRowUpper(ii).getElements()[0] == demand[index+ii]);
			//printf("event prob %g\n",smiRV->getEventProb(ii));
			assert(fabs(smiRV->getEventProb(ii) - dprobs[index+ii]) < 0.0000001);
		}
		smiDD->addDiscreteRV(smiRV);
		index+=nsamp[jj];
	}

	assert(smiDD->getNumRV() == nindp);

	smiModel->processDiscreteDistributionIntoScenarios(smiDD);

	
	// load problem data into OsiSolver
	smiModel->loadOsiSolverData();
	// get Osi pointer
	OsiSolverInterface *smiOsi = smiModel->getOsiSolverInterface();
	// set some parameters
	smiOsi->setHintParam(OsiDoPresolveInInitial,true);
	smiOsi->setHintParam(OsiDoScale,true);
	smiOsi->setHintParam(OsiDoCrash,true);
	// solve using Osi Solver
	smiOsi->initialSolve();
	// test optimal value
    	assert(fabs(smiOsi->getObjValue()-1566.042)<0.01);

	// test solutions
	const double *dsoln = smiOsi->getColSolution();
	double objSum = 0.0;

	/* The canonical way to traverse the tree:
	   For each scenario, get the leaf node.
	   Then get the parent.  Repeat until parent is NULL.
	   (Only the root node has a NULL parent.)
	 */

	for(int is=0; is<smiModel->getNumScenarios(); ++is)
	{
		/* this loop calculates the scenario objective value */
		double scenSum = 0.0;

		// start with leaf node
		SmiScnNode *node = smiModel->getLeafNode(is);

		// leaf node probability is the scenario probability
		double scenprob = node->getModelProb();
	
		while (node != NULL)
		{
			// getColStart returns the starting index of node in OSI model
			for(int j=node->getColStart(); j<node->getColStart()+node->getNumCols(); ++j)
			{
				// getCoreColIndex returns the corresponding Core index
				// in the original (user's) ordering
				scenSum += dobj[node->getCoreColIndex(j)]*dsoln[j];	
				

			}			
			
			// get parent of node
			node = node->getParent();
		}
		objSum += scenSum*scenprob;
	}

	assert(fabs(smiOsi->getObjValue()-objSum) < 0.01);

		// print results
		printf("Solved stochastic program %s\n", name);
		printf("Number of rows: %d\n",smiOsi->getNumRows());
		printf("Number of cols: %d\n",smiOsi->getNumCols());
		printf("Optimal value: %g\n",smiOsi->getObjValue());


}

// Display message on stdout and stderr
void testingMessage( const char * const msg )
{
//  std::cerr <<msg;
  cout <<endl <<"*****************************************"
       <<endl <<msg <<endl;
}