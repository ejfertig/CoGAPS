#ifndef _GIBBSSAMPLER_H_
#define _GIBBSSAMPLER_H_

#include<iostream>
#include<fstream>
#include <string>
#include <vector>
#include <algorithm>
#include "GAPSNorm.h"
#include "randgen.h"
#include "sub_func.h"
#include "Matrix.h"
#include "AtomicSupport.h"
#include<limits>

using namespace gaps;

class GibbsSampler
{
 protected:
  // Parameters or data that are read in:
  unsigned long _nEquil;    // # outer loop iterations for equilibration
  unsigned long _nSample;   // # outer loop iterations for sampling
  unsigned int _nFactor;    // # patterns
  double _alphaA;
  double _alphaP;
  double _nMaxA;             // max. number of atoms in A
  double _nMaxP;             // number of atomic bins for P
  unsigned long _nIterA;    // initial # of inner-loop iterations for A 
  unsigned long _nIterP ;    // initial # of inner-loop iterations for P 
  string _simulation_id;   // simulation id for the run
  double _max_gibbsmass_paraA; // max gibbs mass parameter for A
  double _max_gibbsmass_paraP; // max gibbs mass parameter for P
  double _lambdaA_scale_factor; // factor to rescale _lambdaA
  double _lambdaP_scale_factor; // factor to rescale _lambdaP

  // Parameters or structures to be calculated or constructed:
  unsigned int _nRow;       // number of items in observation (= # of genes)
  unsigned int _nCol;       // number of observation (= # of arrays)
  unsigned int _nBinsA;     // number of atomic bins for A
  unsigned int _nBinsP;     // number of atomic bins for P
  double _lambdaA;
  double _lambdaP;
  double _max_gibbsmassA;  // max gibbs mass for A
  double _max_gibbsmassP;  // max gibbs mass for P
  unsigned long _atomicSize; // number of atomic points

  char _label_A;  // label for matrix A
  char _label_P;  // label for matrix P
  char _label_D;  // label for matrix D
  char _label_S;  // label for matrix S

  unsigned long _iter;
  double _annealingTemperature;

  AtomicSupport _AAtomicdomain,_PAtomicdomain;
  Matrix _AMatrix,_PMatrix,_DMatrix,_SMatrix;

  map<unsigned long long, double> _atomicProposal;
  unsigned int _nChange_atomicProposal;
  char _oper_type;

  unsigned int _nChange_matrixElemChange;
  vector<unsigned int> _Row_changed;
  vector<unsigned int> _Col_changed;
  vector<double> _mass_changed;
  vector<boost::tuple<unsigned int, unsigned int, double> > _matrixElemChange;
  // vector<boost::tuple<unsigned int, unsigned int, double> > _ElemChange;


  map<unsigned long long, double> _new_atomicProposal;
  unsigned int _new_nChange_atomicProposal;

  unsigned int _new_nChange_matrixElemChange;
  vector<unsigned int> _new_Row_changed;
  vector<unsigned int> _new_Col_changed;
  vector<double> _new_mass_changed;
  vector<boost::tuple<unsigned int, unsigned int, double> > _new_matrixElemChange;

  // for computing statistics with matrices A and P
  // unsigned long _statindx_A, _statindx_P;  // counter
  double ** _Amean;
  double ** _Asd; 
  double ** _Pmean; 
  double ** _Psd; 

 public:

  // ******************** CONSTRUCTOR ********************************************
  GibbsSampler(){};

  GibbsSampler(unsigned long nEquil, unsigned long nSample, unsigned int nFactor, 
	       double alphaA, double alphaP, double nMaxA, double nMaxP,
	       unsigned long nIterA, unsigned long nIterP, 
	       double max_gibbsmass_paraA, double max_gibbsmass_paraP, 
	       double lambdaA_scale_factor, double lambdaP_scale_factor, 
               unsigned long atomicSize,
	       char label_A,char label_P,char label_D,char label_S,
	       const string & datafile, const string & variancefile,
               const string & simulation_id);

  ~GibbsSampler(){};


  // *************** METHOS FOR INITIALIZATION, DISPALY, OUTPUT ***********************
  void init_AMatrix_and_PMatrix();

  void init_AAtomicdomain_and_PAtomicdomain();

  void clear_Proposal();

  void clear_new_Proposal();

  void display_matrix(char matrix_label);

  void display_atomicdomain(char atomic_label);

  void local_display_matrix(vector<vector<double> > Mat, 
                            unsigned int n_row, unsigned int n_col);

  void local_display_matrix2(double ** Mat_ptr, 
			     unsigned int n_row, unsigned int n_col);

  void local_display_matrix2F(ofstream& outputFile, double ** Mat_ptr, 
			      unsigned int n_row, unsigned int n_col);

  void check_results();

  void check_resultsF(ofstream& outputFile);

  void output_atomicdomain(char atomic_label,unsigned long Samp_cycle);

  void output_computing_info(char outputFilename[],
                             unsigned long Equil_cycle, unsigned long nEquil,
			     unsigned long Samp_cycle, unsigned long nSample, 
                             double chi2);

  // ********* METHODS TO GO BETWEEN ATOMIC SPACE AND MATRIX *****************

  unsigned int getRow( char matrix_label ,unsigned int iBin);

  unsigned int getCol(char matrix_label ,unsigned int iBin);

  unsigned int getTotNumAtoms(char matrix_label);

  // **************** METHODS FOR COMPUTING LIKELIHOOD FUNCTIONS *****************
  double cal_logLikelihood();

  void cal_delloglikelihood_example();

  double computeDeltaLL(char the_matrix_label,
			double const * const * D,
			double const * const * S,
			double const * const * A,
			double const * const * P,
			unsigned int the_nChange_matrixElemChange,
			const vector<boost::tuple<unsigned int, unsigned int, double> > the_matrixElemChange);

  double computeDeltaLL2(char the_matrix_label,
			double const * const * D,
			double const * const * S,
			double const * const * A,
			double const * const * P,
			unsigned int the_nChange_matrixElemChange,
			const vector<boost::tuple<unsigned int, unsigned int, double> > the_matrixElemChange);


  // *************** METHODS FOR MAKING PROPOSAL *********************************
  void update_example(char atomic_domain_label);

  vector<vector<double> > atomicProposal2Matrix(char atomic_domain_label,
  			     double const * const * origMatrix); 
 
  vector<vector<double> > atomicProposal2FullMatrix(char atomic_domain_label,
						    double const * const * origMatrix);

  void extract_atomicProposal(char the_matrix_label);

  void extract_new_atomicProposal(char the_matrix_label);

  void update(char the_matrix_label);

  void get_oper_type(char the_matrix_label);

  bool birth_death(char the_matrix_label, 	
				      double const * const * D,
				      double const * const * S,
				      double ** AOrig,
				      double ** POrig); 

  bool move_exchange(char the_matrix_label,
					double const * const * D,
					double const * const * S,
					double ** AOrig,
					double ** POrig);


  // ************ METHODS FOR LOOPING AND CONTROL ********************************
  void set_iter(unsigned long ext_iter);

  double get_AnnealingTemperature();

  void set_AnnealingTemperature(); 

  void check_atomic_matrix_consistency(char the_matrix_label);

  void compute_statistics_prepare_matrices(unsigned long statindx);

  void compute_statistics(char outputFilename[],
                          char outputAmean_Filename[],char outputAsd_Filename[],
                          char outputPmean_Filename[],char outputPsd_Filename[],
			  char outputAPmean_Filename[],
                          unsigned int Nstat);



  // --------------------------------------------------------------------------- 
  // Adapt directly from old codes:

  bool performUpdate(char the_matrix_label, double origMass, 
		     unsigned int iRow, unsigned int iCol,
		     double const * const * A, double const * const * P);

  bool performUpdateKill(char the_matrix_label, unsigned int iRow, unsigned int iCol,
		       double const * const * otherMatrix);


  double getMass(char the_matrix_label, double origMass,
		 unsigned int iRow,
		 unsigned int iCol,
		 double const * const * otherMatrix, 
		 const vector<vector<double> > currentChainMatrix,
		 double const * const * D, double const * const * S,
		 double rng);
  // ---------------------------------------------------------------------------

  void detail_check(char outputchi2_Filename[]);



};
#endif
