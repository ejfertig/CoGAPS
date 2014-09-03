// cogaps.cpp

// =============================================================================
// This is the main code for CogapsMAP. (7th Sep, 2013)
// This code ONLY incorporates CogapsMAP (for fixed patterns) (7th Aug 2014)
// =============================================================================


#include <iostream>       // for use with standard I/O
#include <fstream>        // for output to files
#include <limits>         // for extracting numerical limits of C++

// ------ incorporated to use Cogaps_options ------------
#include <vector>
#include <iomanip>
#include <boost/program_options.hpp>
#include "Cogaps_options.hpp"
// ------------------------------------------------------
#include "randgen.h";   // for incorporating a random number generator.
#include "Matrix.h";    // for incorporating a Matrix class
#include "AtomicSupport.h";  // for incorporating an Atomic class
#include "GAPSNorm.h";  // for incorporating calculation of statistics in cogaps.
#include "GibbsSampler.h"; // for incorporating the GibbsSampler which
                           // does all the atomic space to matrix conversion
                           // and sampling actions.
#include "GibbsSamplerMap.h"; // for using cogaps with fixed patterns
// ------------------------------------------------------

namespace bpo = boost::program_options;
using namespace std;
using namespace gaps;

boost::mt19937 rng(43);

int main(int ac, char* av[]){

  // global objects for the program:
  Cogaps_options_class Cogaps_options(ac, av);
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

  // Parameters or data to be read in:
  unsigned long nEquil = Cogaps_options.nEquil;    // # outer loop iterations 
                                                   // for equilibration
  unsigned long nSample = Cogaps_options.nSample;  // # outer loop iterations 
                                                   // for sampling
  unsigned int nFactor = Cogaps_options.nFactor;   // # patterns
  double alphaA = Cogaps_options.alphaA;
  double alphaP = Cogaps_options.alphaP;
  double nMaxA = Cogaps_options.nMaxA;             // max. number of atoms in A
  double nMaxP = Cogaps_options.nMaxP;             // number of atomic bins for P
  string datafile = Cogaps_options.datafile;        // File for D
  string variancefile = Cogaps_options.variancefile; // File for S
  string simulation_id = Cogaps_options.simulation_id; // simulation id
  unsigned long nIterA = Cogaps_options.nIterA;    // initial # of inner-loop iterations for A 
  unsigned long nIterP = Cogaps_options.nIterP;    // initial # of inner-loop iterations for P 
  double max_gibbsmass_paraA = Cogaps_options.max_gibbsmass_paraA; 
                           // maximum gibbs mass parameter for A 
  double max_gibbsmass_paraP = Cogaps_options.max_gibbsmass_paraP; 
                           // maximum gibbs mass parameter for P 
  double lambdaA_scale_factor = Cogaps_options.lambdaA_scale_factor;
                           // scale factor for lambdaA
  double lambdaP_scale_factor = Cogaps_options.lambdaP_scale_factor;
                           // scale factor for lambdaP
  bool Q_output_atomic = Cogaps_options.Q_output_atomic;
                           // whether to output the atomic space info
  bool runCogapsMap = Cogaps_options.runCogapsMap;
						   // whether to run cogaps with fixed patterns
  char fixedmatrix = Cogaps_options.fixedmatrix;
						   // the matrix to fix
  string fixedpatternsfile = Cogaps_options.fixedpatternsfile;
						   // the file containing the fixed patterns

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
  outputFile << "Q_output_atomic (bool) = " << Q_output_atomic << endl; 
  outputFile << "runCogapsMap (bool) = " << runCogapsMap << endl;
  if (runCogapsMap){
  outputFile << "Fixed Matrix = " << fixedmatrix << endl;
  outputFile << "Fixed patterns file: " << fixedpatternsfile << endl;
  }
  outputFile << endl;

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

// Running Cogaps with fixed mappings and GibbsSamplerMap class:

// Unique to Fixed Map Cogaps 
// Read the input file and initialize the vector of fixed patterns
    ifstream file;
    file.open(fixedpatternsfile.c_str());
    vector<vector<double> > allPattern;
    string line;
    while(getline(file, line)) {
        stringstream ss(line);
        istream_iterator<double> begin(ss), end;
        vector<double> onePattern(begin, end);
        allPattern.push_back(onePattern);
    }
    file.close();
  
 
  // ---------------------------------------------------------------------------
  // Initialize the GibbsSamplerMap.

  GibbsSamplerMap GibbsSampMap(nEquil,nSample,nFactor,   // construct GibbsSampler and 
                         alphaA,alphaP,nMaxA,nMaxP,// Read in D and S matrices
                         nIterA,nIterP,
			 max_gibbsmass_paraA, max_gibbsmass_paraP, 
			 lambdaA_scale_factor, lambdaP_scale_factor,
                         atomicSize,
                         label_A,label_P,label_D,label_S,
			 datafile,variancefile,simulation_id, allPattern, fixedmatrix);

  // ---------------------------------------------------------------------------
  // Based on the information of D, construct and initialize for A and P both 
  // the matrices and atomic spaces.

  GibbsSampMap.init_AMatrix_and_PMatrix(); // initialize A and P matrices
  GibbsSampMap.init_Mapped_Matrix(); //initialize the mapped matrix
  GibbsSampMap.init_AAtomicdomain_and_PAtomicdomain(); // intialize atomic spaces
                                                      // A and P
  GibbsSampMap.initialize_atomic_domain_map(); // keep the domains in line with matrix
  GibbsSampMap.init_sysChi2(); // initialize the system chi2 value
  GibbsSampMap.print_A_and_P(); // for debugging
  GibbsSampMap.check_atomic_matrix_consistency('A');
  GibbsSampMap.check_atomic_matrix_consistency('P');
  // ===========================================================================
  // Part 2) Equilibration:
  // In this section, we let the system equilibrate with nEquil outer loop 
  // iterations. Within each outer loop iteration, A is iterated nIterA times 
  // and P is iterated nIterP times. After equilibration, we update nIterA and 
  // nIterP according to the expected number of atoms in the atomic spaces 
  // of A and P respectively.
  // ===========================================================================
  
      // --------- temp for initializing output to chi2.txt
      char outputchi2_Filename[80];
      strcpy(outputchi2_Filename,simulation_id.c_str());
      strcat(outputchi2_Filename,"_chi2.txt");
      ofstream outputchi2_File;
      outputchi2_File.open(outputchi2_Filename,ios::out);
      outputchi2_File << "chi2" << endl;
      outputchi2_File.close();
      // --------------
      
  double chi2;

  for (unsigned long ext_iter=1; ext_iter <= nEquil; ++ext_iter){
    GibbsSampMap.set_iter(ext_iter);
    GibbsSampMap.set_AnnealingTemperature();


    for (unsigned long iterA=1; iterA <= nIterA; ++iterA){
      GibbsSampMap.mapUpdate('A');
	  }
    
    GibbsSampMap.check_atomic_matrix_consistency('A');

    for (unsigned long iterP=1; iterP <= nIterP; ++iterP){
      GibbsSampMap.mapUpdate('P');
	}
    GibbsSampMap.check_atomic_matrix_consistency('P');

    // ----------- output computing info ---------
    if ( ext_iter % 20 == 0){
      GibbsSampMap.output_computing_info(outputFilename,ext_iter,nEquil,0,nSample);
      cout << "Equil: " << ext_iter << " of " << nEquil << 
              " ,nA: " << GibbsSampMap.getTotNumAtoms('A') <<
	      " ,nP: " << GibbsSampMap.getTotNumAtoms('P') << 
              " ,System Chi2 = " << GibbsSampMap.get_sysChi2() << endl;
    }

    // -------------------------------------------
    // re-calculate nIterA and nIterP to the expected number of atoms 
    nIterA = (unsigned long) randgen('P',max((double) GibbsSampMap.getTotNumAtoms('A'),10.));
    nIterP = (unsigned long) randgen('P',max((double) GibbsSampMap.getTotNumAtoms('P'),10.));
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
  for (unsigned long ext_iter=1; ext_iter <= nSample; ++ext_iter){
  
    for (unsigned long iterA=1; iterA <= nIterA; ++iterA){
      GibbsSampMap.mapUpdate('A');
	  }
    
    GibbsSampMap.check_atomic_matrix_consistency('A');

    for (unsigned long iterP=1; iterP <= nIterP; ++iterP){ 
      GibbsSampMap.mapUpdate('P');
    }
    GibbsSampMap.check_atomic_matrix_consistency('P');

    if (Q_output_atomic == true){
       GibbsSampMap.output_atomicdomain('A',ext_iter);
       GibbsSampMap.output_atomicdomain('P',ext_iter);
    }

     statindx += 1;
     GibbsSampMap.compute_statistics_prepare_matrices(statindx);

    // ----------- output computing info ---------
    if ( ext_iter % 20 == 0){
	
      GibbsSampMap.output_computing_info(outputFilename,nEquil,nEquil,ext_iter,nSample);
      cout << "Samp: " << ext_iter << " of " << nSample << 
              " ,nA: " << GibbsSampMap.getTotNumAtoms('A') <<
	      " ,nP: " << GibbsSampMap.getTotNumAtoms('P') << 
              " , System Chi2 = " << GibbsSampMap.get_sysChi2() << endl;

      if (ext_iter == nSample){
         chi2 = 2.*GibbsSampMap.cal_logLikelihood();
	 cout << " *** Check value of final chi2: " << chi2 << " **** " << endl; 
      }


    }

    // -------------------------------------------
    // re-calculate nIterA and nIterP to the expected number of atoms 
    nIterA = (unsigned long) randgen('P',max((double) GibbsSampMap.getTotNumAtoms('A'),10.));
    nIterP = (unsigned long) randgen('P',max((double) GibbsSampMap.getTotNumAtoms('P'),10.));
    // --------------------------------------------

  }  // end of for-block for Sampling


  // ===========================================================================
  // Part 4) Calculate statistics:
  // In this final section, we calculate all statistics pertaining to the final
  // sample and check the results.
  // ===========================================================================
   GibbsSampMap.print_A_and_P();

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



  GibbsSampMap.compute_statistics(outputFilename,
                               outputAmean_Filename,outputAsd_Filename,
			       outputPmean_Filename,outputPsd_Filename,
			       outputAPmean_Filename,
                               statindx);          // compute statistics like mean and s.d.
  

  return 0;

}

