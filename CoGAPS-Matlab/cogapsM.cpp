// cogapsM.cpp

// =============================================================================
// This is the main code for Cogaps. (7th Sep, 2013)
// =============================================================================

#include <iostream>       // for use with standard I/O
#include <fstream>        // for output to files
#include <limits>         // for extracting numerical limits of C++

// ------ incorporated to use Cogaps_options ------------
#include <vector>
#include <iomanip>
//#include <boost/program_options.hpp>
//#include "Cogaps_options.hpp"
// ------------------------------------------------------
#include "randgen.h"   // for incorporating a random number generator.
#include "Matrix.h"    // for incorporating a Matrix class
#include "AtomicSupport.h"  // for incorporating an Atomic class
#include "GAPSNorm.h"  // for incorporating calculation of statistics in cogaps.
#include "GibbsSampler.h" // for incorporating the GibbsSampler which
                           // does all the atomic space to matrix conversion
                           // and sampling actions.
#include "mex.h"
// ------------------------------------------------------

//namespace bpo = boost::program_options;
using namespace std;
using namespace gaps;
using std::vector;

boost::mt19937 rng(43);

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
//int main(int ac, char* av[]){

  // global objects for the program:
  /*Cogaps_options_class Cogaps_options(ac, av);
  try {
    // Cogaps_options_class Cogaps_options(ac, av); // WS -- change from Sasha's
    if (Cogaps_options.to_help){
      Cogaps_options.help(cout);
      return 0;
    }
    //cout << Cogaps_options; 
  } 
  catch( const exception & e) {
    cerr <<e.what()<<endl;
    return 1;
  }
*/
  // ===========================================================================
  // Initialization of the random number generator.
  // Different seeding methods:
  // --- fixed seed 
  //std::vector<unsigned long> ve(2);
  //ve[0]=198782;ve[1]=89082;
  //boost::random::seed_seq seq(ve);
  //rng.seed(seq);
  // --- seeded with time
  rng.seed(static_cast<boost::uint32_t>(std::time(0)));
  //---------------------

  // ===========================================================================
  // Part 1) Initialization: 
  // In this section, we read in the system parameters from the paremter file
  // parameter.txt, and matrices D and S from datafile.txt. 
  // Then we initialize A and P in both their atomic domains and 
  // matrix forms.
  // ===========================================================================
	
	//Matlab pointers for inputs
    mxArray* DArrayElems;
	mxArray* SArrayElems;
	mxArray* InputColsPtr;
	mxArray* InputRowsPtr;

	//Matlab pointers for config 
	mxArray* nFactorPtr;
	mxArray* simulationIDPtr;
	mxArray* nEquilPtr;
	mxArray* nSamplePtr;
	mxArray* QAtomicPtr;
	mxArray* alphaAPtr;
	mxArray* nIterAPtr;
	mxArray* nMaxAPtr;
	mxArray* maxGibbsAPtr;
	mxArray* lambdaAScalePtr;
	mxArray* alphaPPtr;
	mxArray* nIterPPtr;
	mxArray* nMaxPPtr;
	mxArray* maxGibbsPPtr;
	mxArray* lambdaPScalePtr;
	
	//Transferring inputs from mex gateway to pointers
	DArrayElems = mxGetCell(prhs[0], 0);
	SArrayElems = mxGetCell(prhs[0], 1);
	InputColsPtr = mxGetCell(prhs[0], 2);
	InputRowsPtr = mxGetCell(prhs[0], 3);
	nFactorPtr = mxGetCell(prhs[0], 4);
	simulationIDPtr = mxGetCell(prhs[0], 5);
	nEquilPtr = mxGetCell(prhs[0], 6);
	nSamplePtr = mxGetCell(prhs[0], 7);
	QAtomicPtr = mxGetCell(prhs[0], 8);
	alphaAPtr = mxGetCell(prhs[0], 9);
	nIterAPtr = mxGetCell(prhs[0], 10);
	nMaxAPtr = mxGetCell(prhs[0], 11);
	maxGibbsAPtr = mxGetCell(prhs[0], 12);
	lambdaAScalePtr = mxGetCell(prhs[0], 13);
	alphaPPtr = mxGetCell(prhs[0], 14);
	nIterPPtr = mxGetCell(prhs[0], 15);
	nMaxPPtr = mxGetCell(prhs[0], 16);
	maxGibbsPPtr = mxGetCell(prhs[0], 17);
	lambdaPScalePtr = mxGetCell(prhs[0], 18);
	
	
	//creating native C style data types
	double* DArray;
	double* SArray;
	double* numInputRows;
	double* numInputCols;
	
	//Establish the C types for the Config Information
	double* nFactorC = mxGetPr(nFactorPtr);
	
	int ptrLength; 
	char* simulationIDC;
	ptrLength = mxGetNumberOfElements(simulationIDPtr) + 1;
	simulationIDC = (char *)(mxCalloc(ptrLength, sizeof(char)));
	mxGetString(simulationIDPtr, simulationIDC, ptrLength);
	
	double* nEquilC = mxGetPr(nEquilPtr);
	double* nSampleC = mxGetPr(nSamplePtr);
	double* QAtomicC = mxGetPr(QAtomicPtr);
	double* alphaAC = mxGetPr(alphaAPtr);
	double* nIterAC = mxGetPr(nIterAPtr);
	double* nMaxAC = mxGetPr(nMaxAPtr);
	double* maxGibbsAC = mxGetPr(maxGibbsAPtr);
	double* lambdaAScaleC = mxGetPr(lambdaAScalePtr);
	double*	alphaPC = mxGetPr(alphaPPtr);
	double*	nIterPC = mxGetPr(nIterPPtr);
	double* nMaxPC = mxGetPr(nMaxPPtr);
	double* maxGibbsPC = mxGetPr(maxGibbsPPtr);
	double* lambdaPScaleC = mxGetPr(lambdaPScalePtr);
	
	int nIRows;
	int nICols;
	
	
	
	//converting from pointer to matlab type then to C array
	DArray = mxGetPr(DArrayElems);
	SArray = mxGetPr(SArrayElems);
	numInputCols = mxGetPr(InputRowsPtr);
	numInputRows = mxGetPr(InputColsPtr);
	
	//removing need for arrays for dimensions
	nIRows = numInputRows[0];
	nICols = numInputCols[0];
	
	//Putting C style arrays into C++ vectors to pass to the Gibbs Sampler
	vector< vector<double> > DVector(nIRows, vector<double>(nICols));
	vector< vector<double> > SVector(nIRows, vector<double>(nICols));
	
	for(int i = 0; i < nIRows; i++)
	{
		for(int j = 0; j < nICols; j++)
		{
			DVector[i][j] = DArray[i + (j*nIRows)];
		}
	}
	
	for(int i = 0; i < nIRows; i++)
	{
		for(int j = 0; j < nICols; j++)
		{
			SVector[i][j] = SArray[i + (j*nIRows)];
		}
	}
	
	/*Testing output of matrices
	cout << endl;
	for(int i = 0; i < nIRows; i++)
	{
		for(int j = 0; j < nICols; j++)
		{
			cout << DVector[i][j] << " ";
		}
		cout << endl;
	}
	
	cout << endl;
	for(int i = 0; i < nIRows; i++)
	{
		for(int j = 0; j < nICols; j++)
		{
			cout << SVector[i][j] << " ";
		}
		cout << endl;
	}
	*/
	
  // Parameters or data to be read in:
  unsigned long nEquil = nEquilC[0]; //Cogaps_options.nEquil;    // # outer loop iterations 
                                                   // for equilibration
  unsigned long nSample = nSampleC[0];//Cogaps_options.nSample;  // # outer loop iterations 
                                                   // for sampling
  unsigned int nFactor = nFactorC[0]; //Cogaps_options.nFactor;   // # patterns
  double alphaA = alphaAC[0]; //Cogaps_options.alphaA;
  double alphaP = alphaPC[0]; //Cogaps_options.alphaP;
  double nMaxA = nMaxAC[0];//Cogaps_options.nMaxA;             // max. number of atoms in A
  double nMaxP = nMaxPC[0];//Cogaps_options.nMaxP;             // number of atomic bins for P
  //string datafile = "Data2.txt";//Cogaps_options.datafile;        // File for D
  //string variancefile = "Noise2.txt";//Cogaps_options.variancefile; // File for S
  string simulation_id = simulationIDC;//Cogaps_options.simulation_id; // simulation id
  unsigned long nIterA = nIterAC[0];//Cogaps_options.nIterA;    // initial # of inner-loop iterations for A 
  unsigned long nIterP = nIterPC[0];//Cogaps_options.nIterP;    // initial # of inner-loop iterations for P 
  double max_gibbsmass_paraA = maxGibbsAC[0];//Cogaps_options.max_gibbsmass_paraA; 
                           // maximum gibbs mass parameter for A 
  double max_gibbsmass_paraP = maxGibbsPC[0];//Cogaps_options.max_gibbsmass_paraP; 
                           // maximum gibbs mass parameter for P 
  double lambdaA_scale_factor = lambdaAScaleC[0];//Cogaps_options.lambdaA_scale_factor;
                           // scale factor for lambdaA
  double lambdaP_scale_factor = lambdaPScaleC[0];//Cogaps_options.lambdaP_scale_factor;
                           // scale factor for lambdaP
  bool Q_output_atomic = QAtomicC[0];//Cogaps_options.Q_output_atomic;
                           // whether to output the atomic space info

  // Parameters or structures to be calculated or constructed:
  unsigned int nRow;       // number of items in observation (= # of genes)
  unsigned int nCol;       // number of observation (= # of arrays)
  unsigned int nBinsA;     // number of atomic bins for A
  unsigned int nBinsP;     // number of atomic bins for P
  double lambdaA;
  double lambdaP;
  //unsigned long nIterA;    // number of inner loop iterations for A
  //unsigned long nIterP;    // number of inner loop iterations for P  
  //atomic At, Pt;           // atomic space for A and P respectively
  unsigned long atomicSize; // number of atomic points

  char label_A = 'A';  // label for matrix A
  char label_P = 'P';  // label for matrix P
  char label_D = 'D';  // label for matrix D
  char label_S = 'S';// label for matrix S

  // Output parameters and computing info to files:
  /* Commented out for Matlab and R Versions
  char outputFilename[80];
  strcpy(outputFilename,simulation_id.c_str());
  strcat(outputFilename,"_computing_info.txt");
  
 
  ofstream outputFile;
  outputFile.open(outputFilename,ios::out);  // start by deleting previous content and 
                                             // rewriting the file
  outputFile << "Common parameters and info:" << endl;
  outputFile << "input data file: " << datafile << endl;
  outputFile << "input variance file: " << variancefile << endl;
  outputFile << "simulation id: " << simulation_id << endl;
  outputFile << "nFactor = " << nFactor << endl;
  outputFile << "nEquil = " << nEquil << endl;
  outputFile << "nSample = " << nSample << endl;
  outputFile << "Q_output_atomic (bool) = " << Q_output_atomic << endl << endl;

  outputFile << "Parameters for A:" << endl;
  outputFile << "alphaA = " << alphaA << endl;
  outputFile << "nMaxA = " << nMaxA << endl;
  outputFile << "nIterA = " << nIterA << endl;
  outputFile << "max_gibbsmass_paraA = " << max_gibbsmass_paraA << endl;
  outputFile << "lambdaA_scale_factor = " << lambdaA_scale_factor << endl << endl;

  outputFile << "Parameters for P:" << endl;
  outputFile << "alphaP = " << alphaP << endl;
  outputFile << "nMaxP = " << nMaxP << endl;
  outputFile << "nIterP = " << nIterP << endl;
  outputFile << "max_gibbsmass_paraP = " << max_gibbsmass_paraP << endl;
  outputFile << "lambdaP_scale_factor = " << lambdaP_scale_factor << endl << endl;

  outputFile.close();
  */



  // ---------------------------------------------------------------------------
  // Initialize the GibbsSampler.

  /*Regular Version
  GibbsSampler GibbsSamp(nEquil,nSample,nFactor,   // construct GibbsSampler and 
                         alphaA,alphaP,nMaxA,nMaxP,// Read in D and S matrices
                         nIterA,nIterP,
			 max_gibbsmass_paraA, max_gibbsmass_paraP, 
			 lambdaA_scale_factor, lambdaP_scale_factor,
                         atomicSize,
                         label_A,label_P,label_D,label_S,
			 datafile,variancefile,simulation_id);
	*/
	
	//R and Matlab Version
	GibbsSampler GibbsSamp(nEquil,nSample,nFactor,   // construct GibbsSampler and 
                         alphaA,alphaP,nMaxA,nMaxP,// Read in D and S matrices
                         nIterA,nIterP,
			 max_gibbsmass_paraA, max_gibbsmass_paraP, 
			 lambdaA_scale_factor, lambdaP_scale_factor,
                         atomicSize,
                         label_A,label_P,label_D,label_S,
			 DVector,SVector,simulation_id);

  // ---------------------------------------------------------------------------
  // Based on the information of D, construct and initialize for A and P both 
  // the matrices and atomic spaces.

  GibbsSamp.init_AMatrix_and_PMatrix(); // initialize A and P matrices
  GibbsSamp.init_AAtomicdomain_and_PAtomicdomain(); // intialize atomic spaces
                                                    // A and P
  GibbsSamp.init_sysChi2(); // initialize the system chi2 value

  // ===========================================================================
  // Part 2) Equilibration:
  // In this section, we let the system eqilibrate with nEquil outer loop 
  // iterations. Within each outer loop iteration, A is iterated nIterA times 
  // and P is iterated nIterP times. After equilibration, we update nIterA and 
  // nIterP according to the expected number of atoms in the atomic spaces 
  // of A and P respectively.
  // ===========================================================================
  
      // --------- temp for initializing output to chi2.txt
	  /*
      char outputchi2_Filename[80];
      strcpy(outputchi2_Filename,simulation_id.c_str());
      strcat(outputchi2_Filename,"_chi2.txt");
      ofstream outputchi2_File;
      outputchi2_File.open(outputchi2_Filename,ios::out);
      outputchi2_File << "chi2" << endl;
      outputchi2_File.close();
	  */
      // --------------
      


  double chi2;
  //Matlab control variables 
  double tempChiSq;
  double tempAtomA;
  double tempAtomP;
	int outCount = 0;
	int numOutputs = 100; 
	int totalChiSize = nSample + nEquil;
	
	//Establishing Matlab containers for atoms and chisquare
	double* nAEquil;
	mxArray* nAEMx;
	nAEMx = mxCreateDoubleMatrix(1, nEquil, mxREAL);
	nAEquil = mxGetPr(nAEMx);
	
	double* nASamp;
	mxArray* nASMx;
	nASMx = mxCreateDoubleMatrix(1, nSample, mxREAL);
	nASamp = mxGetPr(nASMx);
	
	double* nPEquil;
	mxArray* nPEMx;
	nPEMx = mxCreateDoubleMatrix(1, nEquil, mxREAL);
	nPEquil = mxGetPr(nPEMx);
	
	double* nPSamp;
	mxArray* nPSMx;
	nPSMx = mxCreateDoubleMatrix(1, nSample, mxREAL);
	nPSamp = mxGetPr(nPSMx);
	
	double* chiVect;
	mxArray* chiMx;
	chiMx = mxCreateDoubleMatrix(1, totalChiSize, mxREAL);
	chiVect = mxGetPr(chiMx);
	

  for (unsigned long ext_iter=1; ext_iter <= nEquil; ++ext_iter){
    GibbsSamp.set_iter(ext_iter);
    GibbsSamp.set_AnnealingTemperature();


    for (unsigned long iterA=1; iterA <= nIterA; ++iterA){
      GibbsSamp.update('A');
	  GibbsSamp.check_atomic_matrix_consistency('A');
      //GibbsSamp.check_atomic_matrix_consistency('A');
      //GibbsSamp.detail_check(outputchi2_Filename);
    }
    

    for (unsigned long iterP=1; iterP <= nIterP; ++iterP){
      GibbsSamp.update('P');
	  GibbsSamp.check_atomic_matrix_consistency('P');

      //GibbsSamp.check_atomic_matrix_consistency('P');
      //GibbsSamp.detail_check(outputchi2_Filename);
    }
	tempChiSq = GibbsSamp.get_sysChi2();
	chiVect[(ext_iter)-1] = tempChiSq;

	tempAtomA = GibbsSamp.getTotNumAtoms('A');
	tempAtomP = GibbsSamp.getTotNumAtoms('P');
	nAEquil[outCount] = tempAtomA;
	nPEquil[outCount] = tempAtomP;
	outCount++;
    // ----------- output computing info ---------
    if ( ext_iter % numOutputs == 0){
      //chi2 = 2.*GibbsSamp.cal_logLikelihood();
      
      cout << "Equil: " << ext_iter << " of " << nEquil << 
              " ,nA: " << tempAtomA <<
	      " ,nP: " << tempAtomP << 
	// " ,chi2 = " << chi2 <<
              " ,System Chi2 = " << tempChiSq << endl;
    }

    // -------------------------------------------
    // re-calculate nIterA and nIterP to the expected number of atoms 
    nIterA = (unsigned long) randgen('P',max((double) GibbsSamp.getTotNumAtoms('A'),10.));
    nIterP = (unsigned long) randgen('P',max((double) GibbsSamp.getTotNumAtoms('P'),10.));
    //nIterA = (unsigned long) randgen('P',(double) GibbsSamp.getTotNumAtoms('A')+10.);
    //nIterP = (unsigned long) randgen('P',(double) GibbsSamp.getTotNumAtoms('P')+10.);
    // --------------------------------------------

  }  // end of for-block for equilibration
 


  // ===========================================================================
  // Part 3) Sampling:
  // After the system equilibriates in Part 2, we sample the systems with an 
  // outer loop of nSample iterations. Within each outer loop iteration, A is 
  // iterated nIterA times and P is iterated nIterP times. After sampling, 
  // we update nIterA and nIterP according to the expected number of atoms in 
  // the atomic spaces of A and P respectively.
  // ===========================================================================
  


  unsigned int statindx = 0;
  outCount = 0;
  for (unsigned long ext_iter=1; ext_iter <= nSample; ++ext_iter){
    for (unsigned long iterA=1; iterA <= nIterA; ++iterA){
      GibbsSamp.update('A');
      //GibbsSamp.check_atomic_matrix_consistency('A');
      //GibbsSamp.detail_check(outputchi2_Filename);
    }
    GibbsSamp.check_atomic_matrix_consistency('A');

    for (unsigned long iterP=1; iterP <= nIterP; ++iterP){ 
      GibbsSamp.update('P');
      //GibbsSamp.check_atomic_matrix_consistency('P');
      //GibbsSamp.detail_check(outputchi2_Filename);
    }
    GibbsSamp.check_atomic_matrix_consistency('P');

    if (Q_output_atomic == true){
       GibbsSamp.output_atomicdomain('A',ext_iter);
       GibbsSamp.output_atomicdomain('P',ext_iter);
    }

     statindx += 1;
     GibbsSamp.compute_statistics_prepare_matrices(statindx);
	 
	tempChiSq = GibbsSamp.get_sysChi2();
	
	chiVect[(nEquil + ext_iter)-1] = tempChiSq;
	
	tempAtomA = GibbsSamp.getTotNumAtoms('A');
	tempAtomP = GibbsSamp.getTotNumAtoms('P');
	nASamp[outCount] = tempAtomA;
	nPSamp[outCount] = tempAtomP;
	outCount++;
    // ----------- output computing info ---------
    if ( ext_iter % numOutputs == 0){
  
      // chi2 = 2.*GibbsSamp.cal_logLikelihood();
      //GibbsSamp.output_atomicdomain('A',(unsigned long) statindx);
      //GibbsSamp.output_atomicdomain('P',(unsigned long) statindx);
		
      cout << "Samp: " << ext_iter << " of " << nSample << 
              " ,nA: " << tempAtomA <<
	      " ,nP: " << tempAtomP << 
	// " ,chi2 = " << chi2 <<
              " , System Chi2 = " << tempChiSq << endl;
		
      if (ext_iter == nSample){
         chi2 = 2.*GibbsSamp.cal_logLikelihood();
	 cout << " *** Check value of final chi2: " << chi2 << " **** " << endl; 
      }



    }

    // -------------------------------------------
    // re-calculate nIterA and nIterP to the expected number of atoms 
    nIterA = (unsigned long) randgen('P',max((double) GibbsSamp.getTotNumAtoms('A'),10.));
    nIterP = (unsigned long) randgen('P',max((double) GibbsSamp.getTotNumAtoms('P'),10.));
    //nIterA = (unsigned long) randgen('P',(double) GibbsSamp.getTotNumAtoms('A')+10.);
    //nIterP = (unsigned long) randgen('P',(double) GibbsSamp.getTotNumAtoms('P')+10.);
    // --------------------------------------------

  }  // end of for-block for Sampling

 

  // ===========================================================================
  // Part 4) Calculate statistics:
  // In this final section, we calculate all statistics pertaining to the final
  // sample and check the results.
  // ===========================================================================

  char outputAmean_Filename[80];
  strcpy(outputAmean_Filename,simulation_id.c_str());
  strcat(outputAmean_Filename,"_Amean.txt");

  char outputAsd_Filename[80];
  strcpy(outputAsd_Filename,simulation_id.c_str());
  strcat(outputAsd_Filename,"_Asd.txt");

  char outputPmean_Filename[80];
  strcpy(outputPmean_Filename,simulation_id.c_str());
  strcat(outputPmean_Filename,"_Pmean.txt");

  char outputPsd_Filename[80];
  strcpy(outputPsd_Filename,simulation_id.c_str());
  strcat(outputPsd_Filename,"_Psd.txt");

  char outputAPmean_Filename[80];
  strcpy(outputAPmean_Filename,simulation_id.c_str());
  strcat(outputAPmean_Filename,"_APmean.txt");



 /* GibbsSamp.compute_statistics(outputFilename,
                               outputAmean_Filename,outputAsd_Filename,
			       outputPmean_Filename,outputPsd_Filename,
			       outputAPmean_Filename,
                               statindx);          // compute statistics like mean and s.d.n */
							   
	vector<vector <double> > AMeanVector;
	vector<vector <double> > AStdVector;
	vector<vector <double> > PMeanVector;
	vector<vector <double> > PStdVector;
	
	GibbsSamp.compute_statistics(statindx,
							   AMeanVector, AStdVector, PMeanVector, PStdVector);   
	
	/*
	cout << endl;
	for(int i = 0; i < AMeanVector.size(); i++)
	{
		for(int j = 0; j < AMeanVector[0].size(); j++)
		{
			cout << AMeanVector[i][j] << " "; 
		}
		cout << endl;
	}
	*/
	
	/*double *x,*y;
	size_t mrows,ncols;
	mrows = mxGetM(prhs[0]);
	ncols = mxGetN(prhs[0]);
	plhs[0] = mxCreateDoubleMatrix((mwSize)mrows, (mwSize)ncols, mxREAL);
	x = mxGetPr(prhs[0]);
	y = mxGetPr(plhs[0]);*/
	
	//Get Dims of Matrices for matlab
	int nrowA, ncolA;
	int nrowP, ncolP;
	ncolA = AMeanVector.size();
	nrowA = AMeanVector[0].size();
	ncolP = PMeanVector.size();
	nrowP = PMeanVector[0].size();
	
	//Declare Matrices for matlab
	double* AMeanMatArray;
	double* AStdMatArray;
	double* PMeanMatArray;
	double* PStdMatArray;
	
	//Convertible Matlab datatypes
	mxArray* AMeanMx;
	mxArray* AStdMx;
	mxArray* PMeanMx;
	mxArray* PStdMx;
	
	//Allocate memory for the Matlab Matrices and establish their data types
	AMeanMx = mxCreateDoubleMatrix(ncolA, nrowA, mxREAL);
	AStdMx = mxCreateDoubleMatrix(ncolA, nrowA, mxREAL);
	PMeanMx = mxCreateDoubleMatrix(ncolP, nrowP, mxREAL);
	PStdMx = mxCreateDoubleMatrix(ncolP, nrowP, mxREAL);
	
	
	//Steps to Create the Cell Array to pass back the data matrices 
	mxArray* TempArry;
	const int* dims;
	TempArry = mxCreateDoubleMatrix(9, 9, mxREAL); //TODO, currently creating the dimensions of a cell Matrix by creating a temp matrix due to lack of understanding of mex datatypes, fix!
	dims = mxGetDimensions(TempArry);
	plhs[0] = mxCreateCellArray(1, dims);
	
	

	 //Fill the matrices individually 
	AMeanMatArray = mxGetPr(AMeanMx);
	double tempVectValue;
	for(int i = 0; i < ncolA; i++)
	{
		for(int j = 0; j < nrowA; j++)
		{
			tempVectValue = AMeanVector[i][j];
			AMeanMatArray[i + (j*ncolA)] = tempVectValue;
		}
	}
	
	AStdMatArray = mxGetPr(AStdMx);
	for(int i = 0; i < ncolA; i++)
	{
		for(int j = 0; j < nrowA; j++)
		{
			tempVectValue = AStdVector[i][j];
			AStdMatArray[i + (j*ncolA)] = tempVectValue;
		}
	}
	
	PMeanMatArray = mxGetPr(PMeanMx);
	for(int i = 0; i < ncolP; i++)
	{
		for(int j = 0; j < nrowP; j++)
		{
			tempVectValue = PMeanVector[i][j];
			PMeanMatArray[i + (j*ncolP)] = tempVectValue;
		}
	}
	
	PStdMatArray = mxGetPr(PStdMx);
	for(int i = 0; i < ncolP; i++)
	{
		for(int j = 0; j < nrowP; j++)
		{
			tempVectValue = PStdVector[i][j];
			PStdMatArray[i + (j*ncolP)] = tempVectValue;
		}
	}
	
	mxSetCell(plhs[0], 0, AMeanMx);
	mxSetCell(plhs[0], 1, AStdMx);
	mxSetCell(plhs[0], 2, PMeanMx);
	mxSetCell(plhs[0], 3, PStdMx);
	mxSetCell(plhs[0], 4, nAEMx);
	mxSetCell(plhs[0], 5, nASMx);
	mxSetCell(plhs[0], 6, nPEMx);
	mxSetCell(plhs[0], 7, nPSMx);
	mxSetCell(plhs[0], 8, chiMx);
	
}
