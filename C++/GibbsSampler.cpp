#include <iostream>
#include <cmath>
#include <limits>
#include <stdexcept>
#include "GibbsSampler.h"

using namespace std;
using namespace gaps;
using std::vector;

// -----------------------------------------------------------------------------
unsigned long atomicSize = std::numeric_limits<unsigned long>::max();  
const double DOUBLE_POSINF = std::numeric_limits<double>::max();
const double DOUBLE_NEGINF = -std::numeric_limits<double>::max();
const double epsilon = 1e-10;
// -----------------------------------------------------------------------------



// ******************** CONSTRUCTOR ********************************************
GibbsSampler:: GibbsSampler(unsigned long nEquil, unsigned long nSample, unsigned int nFactor, 
			    double alphaA, double alphaP, double nMaxA, double nMaxP,
			    unsigned long nIterA, unsigned long nIterP, 
                            double max_gibbsmass_paraA, double max_gibbsmass_paraP,
                            double lambdaA_scale_factor, double lambdaP_scale_factor,
                            unsigned long atomicSize,
			    char label_A,char label_P,char label_D,char label_S,
			    const string & datafile, const string & variancefile,
                            const string & simulation_id)
  :_DMatrix(datafile.c_str(),label_D),
   _SMatrix(variancefile.c_str(),label_S){
 
  _nEquil = nEquil;
  _nSample = nSample;
  _nFactor = nFactor;
  _alphaA = alphaA;
  _alphaP = alphaP;
  _nMaxA = nMaxA;
  _nMaxP = nMaxP;
  _nIterA = nIterA;
  _nIterP = nIterP;
  _max_gibbsmass_paraA = max_gibbsmass_paraA;
  _max_gibbsmass_paraP = max_gibbsmass_paraP;
  _lambdaA_scale_factor = lambdaA_scale_factor;
  _lambdaP_scale_factor = lambdaP_scale_factor;
  _simulation_id = simulation_id;
  _atomicSize = atomicSize;
  _label_A = label_A;
  _label_P = label_P;
  _iter = 1;  // tmp use
  _annealingTemperature = 0.4; // tmp use
  _sysChi2 = 0.0; // tmp use
 
}


// *************** METHOS FOR INITIALIZATION, DISPALY, OUTPUT ***********************
void GibbsSampler::init_AMatrix_and_PMatrix(){
  // extract information from D as parameters
  _nRow = _DMatrix.get_nRow();
  _nCol = _DMatrix.get_nCol();

  // initialize matrices A and p
  _AMatrix.born_matrix(_nRow,_nFactor,_label_A,_alphaA);
  _PMatrix.born_matrix(_nFactor,_nCol,_label_P,_alphaP);
}

void GibbsSampler::init_AAtomicdomain_and_PAtomicdomain(){
  // extract information from D as parameters
  _nRow = _DMatrix.get_nRow();
  _nCol = _DMatrix.get_nCol();
  double D_mean = _DMatrix.cal_mean();

  // calcuate #Bins and lambda for the atomic spaces
  _nBinsA = _nRow*_nFactor;
  _lambdaA = _alphaA*sqrt(_nFactor/D_mean) * _lambdaA_scale_factor;
  _nBinsP = _nFactor*_nCol;
  _lambdaP = _alphaP*sqrt(_nFactor/D_mean) * _lambdaP_scale_factor;

  // calculate the maximum gibbs mass for A and p
  _max_gibbsmassA = _max_gibbsmass_paraA / _lambdaA;
  _max_gibbsmassP = _max_gibbsmass_paraP / _lambdaP;

  // initialize the atomic spaces
  _AAtomicdomain.initializeAtomic(_nBinsA,atomicSize,_alphaA,_lambdaA,_label_A);
  _PAtomicdomain.initializeAtomic(_nBinsP,atomicSize,_alphaP,_lambdaP,_label_P);

  cout << "_lambdaA = " << _lambdaA << ", _max_gibbsmassA = " << _max_gibbsmassA << endl;
  cout << "_lambdaP = " << _lambdaP << ", _max_gibbsmassP = " << _max_gibbsmassP << endl << endl;

}

// clear all quantities related to the local matrix proposal
void GibbsSampler::clear_Proposal(){
  _Row_changed.clear();
  _Col_changed.clear();
  _mass_changed.clear();
  _atomicProposal.clear();
  _matrixElemChange.clear();

}

// clear all quantities related to the new local matrix proposal
void GibbsSampler::clear_new_Proposal(){
  _new_Row_changed.clear();
  _new_Col_changed.clear();
  _new_mass_changed.clear();
  _new_atomicProposal.clear();
  _new_matrixElemChange.clear();
}


void GibbsSampler::display_matrix(char matrix_label){
  switch(matrix_label){
  case 'D':
    {_DMatrix.display_matrix();break;}
  case 'S':
    {_SMatrix.display_matrix();break;}
  case 'A':
    {_AMatrix.display_matrix();break;}
  case 'P':
    {_PMatrix.display_matrix();break;}
  }
}



void GibbsSampler::display_atomicdomain(char atomic_label){
  switch(atomic_label){
  case 'A':
    {_AAtomicdomain.printAtomicInfo(); break;}
  case 'P':
    {_PAtomicdomain.printAtomicInfo(); break;}
  }
}

void GibbsSampler::local_display_matrix(vector<vector<double> > Mat, 
					unsigned int n_row, unsigned int n_col)
{
  cout << endl;
  for(unsigned int m=0;m<n_row;++m)
    {
      for (unsigned int n=0; n<n_col;++n)
	{
	  cout << std::setw(10) << std::right;
	  cout << Mat[m][n] << " ";
	}
      cout << endl;
    }
  cout << endl;
}


void GibbsSampler::local_display_matrix2(double ** Mat_ptr, 
					 unsigned int n_row, unsigned int n_col)
{
  cout << endl;
  for(unsigned int m=0;m<n_row;++m)
    {
      for (unsigned int n=0; n<n_col;++n)
	{
	  cout << std::setw(10) << std::right; 
	  cout << Mat_ptr[m][n] << " ";
	}
      cout << endl;
    }
  cout << endl;
}

void GibbsSampler::local_display_matrix2F(ofstream& outputFile, double ** Mat_ptr, 
					  unsigned int n_row, unsigned int n_col){
  //outputFile << endl;
  for(unsigned int m=0;m<n_row;++m)
    {
      for (unsigned int n=0; n<n_col;++n)
	{
	  outputFile << std::setw(10) << std::right; 
	  outputFile << Mat_ptr[m][n] << " ";
	}
      outputFile << endl;
    }
  //outputFile << endl;

}



// -----------------------------------------------------------------------------
void GibbsSampler::check_results(){
  double const * const * D = _DMatrix.get_matrix();
  double const * const * S = _SMatrix.get_matrix();
  double const * const * A = _AMatrix.get_matrix();
  double const * const * P = _PMatrix.get_matrix();

  vector<vector<double> > AP;
  AP.resize(_nRow,vector<double>(_nCol,0.0));

  for (unsigned int m=0; m < _nRow; ++m){
    for (unsigned int n=0; n < _nCol; ++n){
      for (unsigned int k=0; k < _nFactor; ++k){
	AP[m][n] += A[m][k]*P[k][n];
      }
    }
  }

  cout << "The product matrix AP = A*P is: " << endl;
  local_display_matrix(AP,_nRow,_nCol);

}

void GibbsSampler::check_resultsF(ofstream& outputFile){
  double const * const * D = _DMatrix.get_matrix();
  double const * const * S = _SMatrix.get_matrix();
  double const * const * A = _AMatrix.get_matrix();
  double const * const * P = _PMatrix.get_matrix();

  vector<vector<double> > AP;
  AP.resize(_nRow,vector<double>(_nCol,0.0));

  for (unsigned int m=0; m < _nRow; ++m){
    for (unsigned int n=0; n < _nCol; ++n){
      for (unsigned int k=0; k < _nFactor; ++k){
	AP[m][n] += A[m][k]*P[k][n];
      }
    }
  }

  outputFile << "The product matrix AP = A*P is: " << endl;
  outputFile << endl;
  for(unsigned int m=0; m < _nRow;++m)
    {
      for (unsigned int n=0; n< _nCol;++n)
	{
	  outputFile << setiosflags(ios::right) << setw(10) << AP[m][n] << " ";
	}
      outputFile << endl;
    }
  outputFile << endl;

}

// -----------------------------------------------------------------------------
void GibbsSampler::output_atomicdomain(char atomic_label,unsigned long Samp_cycle){

  char outputFilename[80]; 
  switch(atomic_label){
  case 'A':
    { 
      strcpy(outputFilename,_simulation_id.c_str());
      strcat(outputFilename,"_A_atomicdomain.txt");
      _AAtomicdomain.writeAtomicInfo(outputFilename,Samp_cycle);
      break;
    }
  case 'P':
    {     
      strcpy(outputFilename,_simulation_id.c_str());
      strcat(outputFilename,"_P_atomicdomain.txt");
      _PAtomicdomain.writeAtomicInfo(outputFilename,Samp_cycle);
      break;
    }
  }
}

// -----------------------------------------------------------------------------
void GibbsSampler::output_computing_info(char outputFilename[],
                                         unsigned long Equil_cycle, unsigned long nEquil,
					 unsigned long Samp_cycle, unsigned long nSample){

  ofstream outputFile;
  outputFile.open(outputFilename,ios::out|ios::app);

  outputFile << " *************************************************** " << endl;
  outputFile << " --------------- NEW ROUND ------------------------- " << endl;
  outputFile << " *************************************************** " << endl << endl;
  outputFile << "Equilibration cycle index = " << Equil_cycle << endl;
  outputFile << "Total number of equilibrating cycles to perform = " <<  nEquil << endl;
  outputFile << "Sampling cycle index = " << Samp_cycle << endl;
  outputFile << "Total number of sampling cycles to perform = " <<  nSample << endl;
  outputFile << "System Chi2-value = " << _sysChi2 << endl;

  //_AAtomicdomain.printAtomicInfoF(outputFile);
  //_PAtomicdomain.printAtomicInfoF(outputFile);
  _AMatrix.display_matrixF(outputFile);
  _PMatrix.display_matrixF(outputFile);
  check_resultsF(outputFile);

  outputFile.close();

}



// ********* METHODS TO GO BETWEEN ATOMIC SPACE AND MATRIX  ********************
unsigned int GibbsSampler::getRow( char matrix_label ,unsigned int iBin){
  switch(matrix_label){ 
  case 'A':  // A - horizontal addressing implicit
    { return (floor(iBin / _nFactor)); break;}
  case 'P':  // P - vertical addressing implicit
    { return iBin % _nFactor; break;}
  }
}

unsigned int GibbsSampler::getCol(char matrix_label ,unsigned int iBin){
  switch(matrix_label){
  case 'A':  // A - horizontal addressing implicit
    { return iBin % _nFactor; break;}
  case 'P':  // P - vertical addressing implicit
    { return floor(iBin / _nFactor); break;}
  }
}


unsigned int GibbsSampler::getTotNumAtoms(char matrix_label){

  switch(matrix_label){
  case 'A':
    { return _AAtomicdomain.getNAtom();
      break;}
  case 'P':
    { return _PAtomicdomain.getNAtom();
      break;}
  }

}



// ************* METHODS FOR COMPUTING LIKELIHOOD FUNCTIONS ******************
// ---------------------------------------------------------------------------
double GibbsSampler::cal_logLikelihood(){
  double ** D = _DMatrix.get_matrix();
  double ** S = _SMatrix.get_matrix();
  double ** A = _AMatrix.get_matrix();
  double ** P = _PMatrix.get_matrix();
  return GAPSNorm::calChi2(D,S,A,P,_nRow,_nCol,_nFactor)/2.;
}

void GibbsSampler::cal_delloglikelihood_example(){

}

// -----------------------------------------------------------------------------
// Comupte the change in likelihood, DeltaLL, by first getting the proposal from
// the atomic space, then it invokes the corresponding methods in GAPSNorm 
// according to the proposal size.
double GibbsSampler::computeDeltaLL(char the_matrix_label,
				    double const * const * D,
				    double const * const * S,
				    double const * const * A,
				    double const * const * P,
				    unsigned int the_nChange_matrixElemChange,
				    const vector<boost::tuple<unsigned int, unsigned int, double> > the_matrixElemChange){

  double DelLL;

 
  switch(the_matrix_label){
  case 'A':
    { 
      if (the_nChange_matrixElemChange == 0){
	DelLL = 0.0;
      } else if (the_nChange_matrixElemChange == 1){
	DelLL = GAPSNorm::calcDeltaLL1E('A',D,S,A,P,the_matrixElemChange,_nRow,
                                        _nCol,_nFactor);
      } else if (the_nChange_matrixElemChange == 2){
	DelLL = GAPSNorm::calcDeltaLL2E('A',D,S,A,P,the_matrixElemChange,_nRow,
                                        _nCol,_nFactor);
      } else {
	DelLL = GAPSNorm::calcDeltaLLGen('A',D,S,A,P,the_matrixElemChange,_nRow,
					 _nCol,_nFactor);
      } // end of if-block according to proposal.size()

      break;} // end of switch-block 'A'
  case 'P':
    { 
      if (the_nChange_matrixElemChange == 0){
	DelLL = 0.0;
      } else if (the_nChange_matrixElemChange == 1){
	DelLL = GAPSNorm::calcDeltaLL1E('P',D,S,A,P,the_matrixElemChange,_nRow,
                                        _nCol,_nFactor);
      } else if (the_nChange_matrixElemChange == 2){
	DelLL = GAPSNorm::calcDeltaLL2E('P',D,S,A,P,the_matrixElemChange,_nRow,
                                        _nCol,_nFactor);
      } else {
	DelLL = GAPSNorm::calcDeltaLLGen('P',D,S,A,P,the_matrixElemChange,_nRow,
					 _nCol,_nFactor);
      } // end of if-block according to proposal.size()

      break;} // end of switch-block 'P' 

  } // end of switch block

  return DelLL;

} // end of computeDeltaLL 


// ---------------- For checking against computeDeltaLL
double GibbsSampler::computeDeltaLL2(char the_matrix_label,
				     double const * const * D,
				     double const * const * S,
				     double const * const * A,
				     double const * const * P,
				     unsigned int the_nChange_matrixElemChange,
				     const vector<boost::tuple<unsigned int, unsigned int, double> > the_matrixElemChange){

  double DelLL;
  switch(the_matrix_label){
  case 'A':
    {
      DelLL = GAPSNorm::calcDeltaLLGen('A',D,S,A,P,the_matrixElemChange,_nRow,
				       _nCol,_nFactor);
      break;}
  case 'P':
    {
      DelLL = GAPSNorm::calcDeltaLLGen('P',D,S,A,P,the_matrixElemChange,_nRow,
				       _nCol,_nFactor);
      break;}
  }
  return DelLL;

} // end of computeDeltaLL2 
// -------------------------------------


// *************** METHODS FOR MAKING PROPOSAL ********************************
// -----------------------------------------------------------------------------
void GibbsSampler::update_example(char atomic_domain_label){

} // end of update_example



// -----------------------------------------------------------------------------
// Construct a "newMatrix" whose mass is given by
// newMatrix = proposed changes mapped from corresponding atomic space 
vector<vector<double> > GibbsSampler::atomicProposal2Matrix(char atomic_domain_label,
							    double const * const * origMatrix)
{  
  unsigned int bin;
  unsigned int chRow, chCol;
  
  switch(atomic_domain_label){
  case 'A':
    { 
      vector<vector<double> > newMatrix(_nRow,vector<double>(_nFactor,0));
      map<unsigned long long, double> proposal = _AAtomicdomain.getProposedAtoms();
      for (map<unsigned long long, double>::const_iterator
	     iter=proposal.begin(); iter != proposal.end(); ++ iter) {
	bin = _AAtomicdomain.getBin(iter->first);
	chRow = getRow('A',bin);
	chCol = getCol('A',bin);
	newMatrix[chRow][chCol] += iter->second;
      }
      return newMatrix;
      break;} // end of case 'A'

  case 'P':
    {
      vector<vector<double> > newMatrix(_nFactor,vector<double>(_nCol,0));
      map<unsigned long long, double> proposal = _PAtomicdomain.getProposedAtoms();
      for (map<unsigned long long, double>::const_iterator
	     iter=proposal.begin(); iter != proposal.end(); ++ iter) {
	bin = _PAtomicdomain.getBin(iter->first);
	chRow = getRow('P',bin);
	chCol = getCol('P',bin);
	newMatrix[chRow][chCol] += iter->second;
      }
      return newMatrix;
      break;} // end of case 'P'

  } // end of switch  
} // end of atomicProposal2Matrix

// -----------------------------------------------------------------------------
// Construct a "FullnewMatrix" whose mass is given by
// newMatrix = origMatrix + proposed changes from corresponding atomic space 
vector<vector<double> > GibbsSampler::atomicProposal2FullMatrix(char atomic_domain_label,
								double const * const * origMatrix)
{  
  unsigned int bin;
  unsigned int chRow, chCol;
  
  switch(atomic_domain_label){
  case 'A':
    { 
      vector<vector<double> > FullnewMatrix(_nRow,vector<double>(_nFactor,0));
      
      for (unsigned int iRow=0; iRow < _nRow; ++iRow){
	for (unsigned int iCol=0; iCol < _nFactor; ++iCol){
	  if (origMatrix[iRow][iCol] < epsilon){
	    FullnewMatrix[iRow][iCol]=0.;
	  } else {
	    FullnewMatrix[iRow][iCol]=origMatrix[iRow][iCol];
	  }
        }
      } // end of for block
      
      map<unsigned long long, double> proposal = _AAtomicdomain.getProposedAtoms();
      for (map<unsigned long long, double>::const_iterator
	     iter=proposal.begin(); iter != proposal.end(); ++ iter) {
	bin = _AAtomicdomain.getBin(iter->first);
	chRow = getRow('A',bin);
	chCol = getCol('A',bin);
	FullnewMatrix[chRow][chCol] += iter->second;
      }
      return FullnewMatrix;
      break;} // end of case 'A'

  case 'P':
    {
      vector<vector<double> > FullnewMatrix(_nFactor,vector<double>(_nCol,0));
      
      for (unsigned int iRow=0; iRow < _nFactor; ++iRow){
	for (unsigned int iCol=0; iCol < _nCol; ++iCol){
	  if (origMatrix[iRow][iCol] < epsilon){
	    FullnewMatrix[iRow][iCol]=0.;
	  } else {
	    FullnewMatrix[iRow][iCol]=origMatrix[iRow][iCol];
	  }
        }
      } // end of for block
      
      map<unsigned long long, double> proposal = _PAtomicdomain.getProposedAtoms();
      for (map<unsigned long long, double>::const_iterator
	     iter=proposal.begin(); iter != proposal.end(); ++ iter) {
	bin = _PAtomicdomain.getBin(iter->first);
	chRow = getRow('P',bin);
	chCol = getCol('P',bin);
	FullnewMatrix[chRow][chCol] += iter->second;
      }
      return FullnewMatrix;
      break;} // end of case 'P'

  } // end of switch  
} // end of atomicProposal2FullMatrix

// ----------------------------------------------------------------------------
// Extract information of the proposal made in the atomic space. 
// Assuming a _atomicProposal, this method instantiates the 
// corresponding member variables for _matrixElemChange. In the current version, 
// we store the new proposal (in matrix space) in two different class variables
// of GibbsSampler:
// 1. vector<unsigned int> _Row_changed; vector<unsigned int> _Col_changed;
//    vector<double> _mass_changed;
// 2. _matrixElemChange (this is a boost::tuple)
// 
void GibbsSampler::extract_atomicProposal(char the_matrix_label){
  unsigned int bin, chRow, chCol;
  double chmass;
  map<unsigned long long, double>::const_iterator iter;

  _nChange_matrixElemChange = 0;
  _nChange_atomicProposal = _atomicProposal.size(); 

  if (_nChange_atomicProposal == 0) {   // atomic proposal size = 0
    _nChange_matrixElemChange = 0;
  } 

  else if (_nChange_atomicProposal == 1){  // atomic proposal size = 1 
    iter = _atomicProposal.begin();
    switch(the_matrix_label){
    case 'A':
      {
	bin = _AAtomicdomain.getBin(iter->first);
	chRow = getRow('A',bin);	  
	chCol = getCol('A',bin);
	break;}
    case 'P':
      {
	bin = _PAtomicdomain.getBin(iter->first);
	chRow = getRow('P',bin);
	chCol = getCol('P',bin);
	break;}
    } // end of switch-block for atomic proposal size = 1
    chmass = iter->second;
    _Row_changed.push_back(chRow);
    _Col_changed.push_back(chCol);
    _mass_changed.push_back(chmass);
    _nChange_matrixElemChange = 1;
    _matrixElemChange.push_back(boost::make_tuple(chRow,chCol,chmass));

  } // end of if-block for atomic proposal size = 1

  else {     // atomic proposal size = 2 or above
    unsigned int count = 0;

    for (iter = _atomicProposal.begin(); iter != _atomicProposal.end(); ++ iter) {

      switch (the_matrix_label){
      case 'A':
	{
	  bin = _AAtomicdomain.getBin(iter->first);
	  chRow = getRow('A',bin);
	  chCol = getCol('A',bin);
	  break;}
      case 'P':
	{
	  bin = _PAtomicdomain.getBin(iter->first);
	  chRow = getRow('P',bin);
	  chCol = getCol('P',bin);
	  break;}
      } // end of switch-block
      chmass = iter->second;

      if (count == 0){    // nothing to check for the first count
	_Row_changed.push_back(chRow);
	_Col_changed.push_back(chCol);
	_mass_changed.push_back(chmass);
	count += 1;
	_nChange_matrixElemChange += 1;
      } else {
	for (unsigned int m = 0; m < count; ++m){
	  if (chRow == _Row_changed[m] && chCol == _Col_changed[m]){
	    _mass_changed[m] += chmass;
	    if (_mass_changed[m] == 0) {	 
	      _nChange_matrixElemChange -= 1;
	      _Row_changed.erase(_Row_changed.begin()+m);
	      _Col_changed.erase(_Col_changed.begin()+m);
	      _mass_changed.erase(_mass_changed.begin()+m);
	    }
	  } else {
	    _Row_changed.push_back(chRow);
	    _Col_changed.push_back(chCol);
	    _mass_changed.push_back(chmass);
	    _nChange_matrixElemChange += 1;
	  } // end of if-block when chRow and chRol refers to new matrix elements
	} // end of for-block when looping through existing elements in _mass_changed
	count = _nChange_matrixElemChange;
      } // end of if-block for count != 0

    } // end of for-block with iter looping through the atomic proposal

    // make up _matrixElemChange
    for (unsigned int m = 0; m<_nChange_matrixElemChange; ++m){
      _matrixElemChange.push_back(boost::make_tuple(_Row_changed[m],_Col_changed[m],_mass_changed[m]));
    }
 
  } // end of if-block for proposal size = 2 or above

} // end of extract_atomicProposal


// -----------------------------------------------------------------------------
// Extract from _new_atomicProposal
// The code is exactly the same as extract_atomicProposal, only using the 
// _new-quantities for the new proposal.

void GibbsSampler::extract_new_atomicProposal(char the_matrix_label){
  unsigned int bin, chRow, chCol;
  double chmass;
  map<unsigned long long, double>::const_iterator iter;

  _new_nChange_matrixElemChange = 0;
  _new_nChange_atomicProposal = _new_atomicProposal.size(); 

  if (_new_nChange_atomicProposal == 0) {    // atomic proposal size = 0
    _new_nChange_matrixElemChange = 0;
  } 

  else if (_new_nChange_atomicProposal == 1){  // atomic proposal size = 1 
    iter = _new_atomicProposal.begin();
    switch(the_matrix_label){
    case 'A':
      {
	bin = _AAtomicdomain.getBin(iter->first);
	chRow = getRow('A',bin);	  
	chCol = getCol('A',bin);
	break;}
    case 'P':
      {
	bin = _PAtomicdomain.getBin(iter->first);
	chRow = getRow('P',bin);
	chCol = getCol('P',bin);
	break;}
    } // end of switch-block for atomic proposal size = 1
    chmass = iter->second;
    _new_Row_changed.push_back(chRow);
    _new_Col_changed.push_back(chCol);
    _new_mass_changed.push_back(chmass);
    _new_nChange_matrixElemChange = 1;
    _new_matrixElemChange.push_back(boost::make_tuple(chRow,chCol,chmass));

  } // end of if-block for atomic proposal size = 1

  else {     // atomic proposal size = 2 or above
    unsigned int count = 0;

    for (iter = _new_atomicProposal.begin(); iter != _new_atomicProposal.end(); ++ iter) {

      switch (the_matrix_label){
      case 'A':
	{
	  bin = _AAtomicdomain.getBin(iter->first);
	  chRow = getRow('A',bin);
	  chCol = getCol('A',bin);
	  break;}
      case 'P':
	{
	  bin = _PAtomicdomain.getBin(iter->first);
	  chRow = getRow('P',bin);
	  chCol = getCol('P',bin);
	  break;}
      } // end of switch-block
      chmass = iter->second;

      if (count == 0){    // nothing to check for the first count
	_new_Row_changed.push_back(chRow);
	_new_Col_changed.push_back(chCol);
	_new_mass_changed.push_back(chmass);
	count += 1;
	_new_nChange_matrixElemChange += 1;
      } else {
	for (unsigned int m = 0; m < count; ++m){
	  if (chRow == _new_Row_changed[m] && chCol == _new_Col_changed[m]){
	    _new_mass_changed[m] += chmass;
	    if (_new_mass_changed[m] == 0) {	 
	      _new_nChange_matrixElemChange -= 1;
	      _new_Row_changed.erase(_new_Row_changed.begin()+m);
	      _new_Col_changed.erase(_new_Col_changed.begin()+m);
	      _new_mass_changed.erase(_new_mass_changed.begin()+m);
	    }
	  } else {
	    _new_Row_changed.push_back(chRow);
	    _new_Col_changed.push_back(chCol);
	    _new_mass_changed.push_back(chmass);
	    _new_nChange_matrixElemChange += 1;
	  } // end of if-block when chRow and chRol refers to new matrix elements
	} // end of for-block when looping through existing elements in _mass_changed
	count = _new_nChange_matrixElemChange;
      } // end of if-block for count != 0

    } // end of for-block with iter looping through the atomic proposal

    // make up _new_matrixElemChange
    for (unsigned int m = 0; m<_new_nChange_matrixElemChange; ++m){
      _new_matrixElemChange.push_back(boost::make_tuple(_new_Row_changed[m],
						      _new_Col_changed[m],_new_mass_changed[m]));
    }
 
  } // end of if-block for proposal size = 2 or above

} // end of extract_new_atomicProposal



// -----------------------------------------------------------------------------
// make proposal and update for the matrices
// It output "newMatrix", which is the final product of all the steps like 
// birth / death, or move / exchange. 
void GibbsSampler::update(char the_matrix_label){
  double rng = 0.1; // no use, filling up the list 
  bool updateIterCount = true;
  double ** D = _DMatrix.get_matrix();
  double ** S = _SMatrix.get_matrix();
  double ** AOrig = _AMatrix.get_matrix();
  double ** POrig = _PMatrix.get_matrix();
  vector<vector<double> > del_matrix;
  map<unsigned long long, double>::iterator iter;

  bool Q_update;
 
  switch(the_matrix_label){
  case 'A':
    {
      // ----------- making a proposal from atomic space A:
      _AAtomicdomain.makeProposal(rng);
      get_oper_type('A');
      _atomicProposal = _AAtomicdomain.getProposedAtoms();
      extract_atomicProposal('A');
      
      // ---- checking, display atomic proposal ----
      //cout << "Inside update(), matrix: A, oper_type = " << _oper_type << endl;    
      //map<unsigned long long, double>::iterator iter;
      //iter = _atomicProposal.begin();
      //cout << "Checking: _AAtomicdomain.getProposedAtoms() gives : " << endl;
      //while(iter != _atomicProposal.end()){
	//cout << "Bin(iter->first) = " << _AAtomicdomain.getBin(iter->first) 
	//   << " , iter->second = " << iter->second << endl;
	//iter++;   
      //}
      if (_nChange_atomicProposal == 1 && (_oper_type =='E' || _oper_type =='M'))
	{cout << "update inconsistency A1! _nChange_atomicProposal = " << _nChange_atomicProposal <<
	         ", _nChange_matrixElemChange = " << _nChange_matrixElemChange << 
	         ", _oper_type = " << _oper_type << endl;}         
      if (_nChange_atomicProposal == 2 && (_oper_type =='D' || _oper_type =='B'))
	{cout << "update inconsistency A2! _nChange_atomicProposal = " << _nChange_atomicProposal << 
 	         ", _nChange_matrixElemChange = " << _nChange_matrixElemChange << 
	         ", _oper_type = " << _oper_type << endl;}   

      // ----------------------------------
      
 

      // the proposal is translated into a proposal to matrix A:
      //vector<vector<double> > proposed_A = atomicProposal2Matrix('A',AOrig);

      // ----------- modifiy the proposal in a Gibbs way:      
      unsigned int iRow, iCol, iFactor;
      if ( _nChange_atomicProposal == 0){}
      if ( _nChange_atomicProposal> 2){
	throw logic_error("GibbsSampler: can't chnage more than two atoms!!");
      }
      if (_nChange_atomicProposal == 2){
	//del_matrix = move_exchange('A',D,S,AOrig,POrig);
        //_AMatrix.matrix_update(del_matrix);
        Q_update = move_exchange('A',D,S,AOrig,POrig);
	//cout << "Q_update = " << Q_update << endl;
	if (Q_update == true){
	   _AMatrix.matrix_Elem_update(_new_matrixElemChange,_oper_type,_new_nChange_matrixElemChange);
	}
      }
      if (_nChange_atomicProposal == 1){
	//del_matrix = birth_death('A',D,S,AOrig,POrig);
	//_AMatrix.matrix_update(del_matrix);	
	Q_update = birth_death('A',D,S,AOrig,POrig);
	//cout << "Q_update = " << Q_update << endl;
	if (Q_update == true){
	   _AMatrix.matrix_Elem_update(_new_matrixElemChange,_oper_type,_new_nChange_matrixElemChange);
	}
      }
 
      break;
    } // end of case 'A'
  case 'P':
    {
      // ----------- making a proposal from atomic space P:
      _PAtomicdomain.makeProposal(rng);
      get_oper_type('P');
      _atomicProposal = _PAtomicdomain.getProposedAtoms();
      extract_atomicProposal('P');    

     // ---- checking, display atomic proposal ----
     //cout << "Inside update(), matrix: P, oper_type = " << _oper_type << endl;     
      //map<unsigned long long, double>::iterator iter;
      //iter = _atomicProposal.begin();
      //cout << "Checking: _PAtomicdomain.getProposedAtoms() gives : " << endl;
      //while(iter != _atomicProposal.end()){
      //cout << "Bin(iter->first) = " << _PAtomicdomain.getBin(iter->first) 
      //   << " , iter->second = " << iter->second << endl;
      //iter++;
      //}
      if (_nChange_atomicProposal == 1 && (_oper_type =='E' || _oper_type =='M'))
	{cout << "update inconsistency P1! _nChange_atomicProposal = " << _nChange_atomicProposal <<
	         ", _nChange_matrixElemChange = " << _nChange_matrixElemChange << 
	         ", _oper_type = " << _oper_type << endl;}         
      if (_nChange_atomicProposal == 2 && (_oper_type =='D' || _oper_type =='B'))
	{cout << "update inconsistency P2! _nChange_atomicProposal = " << _nChange_atomicProposal << 
 	         ", _nChange_matrixElemChange = " << _nChange_matrixElemChange << 
	         ", _oper_type = " << _oper_type << endl;}   

      // ----------------------------------
     


      // the proposal is translated into a proposal to matrix P:
      //vector<vector<double> > proposed_P = atomicProposal2Matrix('P',POrig);

      // ----------- modifiy the proposal in a Gibbs way:     
      unsigned int iRow, iCol, iFactor;
      if (_nChange_atomicProposal== 0){}
      if (_nChange_atomicProposal > 2){
	throw logic_error("GibbsSampler: can't chnage more than two atoms!!");
      }
      if (_nChange_atomicProposal == 2){
	//del_matrix = move_exchange('P',D,S,AOrig,POrig);
	//_PMatrix.matrix_update(del_matrix);
        Q_update = move_exchange('P',D,S,AOrig,POrig);
	//cout << "Q_update = " << Q_update << endl;
	if (Q_update == true){
	   _PMatrix.matrix_Elem_update(_new_matrixElemChange,_oper_type,_new_nChange_matrixElemChange);
	}
      }
      if (_nChange_atomicProposal== 1){
	//del_matrix = birth_death('P',D,S,AOrig,POrig);
	//_PMatrix.matrix_update(del_matrix);
        Q_update = birth_death('P',D,S,AOrig,POrig);
	//cout << "Q_update = " << Q_update << endl;
	if (Q_update == true){
	   _PMatrix.matrix_Elem_update(_new_matrixElemChange,_oper_type,_new_nChange_matrixElemChange);
	}
      }

      break;
    } // end of case 'P'
  } // end of switch block

  /*
  //  Output everything:
  cout << "After all proposal making processes and acceptance, we have " <<
  "the atomic spaces read: " << endl;
  _AAtomicdomain.printAtomicInfo();
  _PAtomicdomain.printAtomicInfo();
  cout << endl;
  cout << "The updated matrices read: " << endl;
  _AMatrix.display_matrix();
  _PMatrix.display_matrix();
 
  cout << endl;
  */
  // clear Proposal for the next run
  clear_Proposal();
  clear_new_Proposal();

} // end of update()

// ----------------------------------------------------------------------------
void GibbsSampler::init_sysChi2(){
  _sysChi2 = 2.*cal_logLikelihood();

}

void GibbsSampler::update_sysChi2(double delsysChi2){
  _sysChi2 -= 2.*delsysChi2;
}

double GibbsSampler::get_sysChi2(){
  return _sysChi2;
}


// -----------------------------------------------------------------------------
void GibbsSampler::get_oper_type(char the_matrix_label){
  switch(the_matrix_label){
  case 'A':
    {
      _oper_type = _AAtomicdomain.get_oper_type();
      break;
    }
  case 'P':
    {
      _oper_type = _PAtomicdomain.get_oper_type();
      break;
    }
  }

}



// -----------------------------------------------------------------------------
bool GibbsSampler::birth_death(char the_matrix_label,
						  double const * const * D,
						  double const * const * S,
						  double ** AOrig,
						  double ** POrig)// in progress
{

  double rng = 0.1; // no use, just fill up the list
  //map<unsigned long long, double> newProposal;
  //map<unsigned long long, double> Proposal;
  //vector<vector<double> > newMatrix, nullMatrix;  // the thing to return in birth_death()
  //vector<vector<double> > FullnewMatrix;
  double newMass = 0;
  double attemptMass = 0;

  /*
  switch(the_matrix_label){
  case 'A':
    { Proposal = _AAtomicdomain.getProposedAtoms();
      newMatrix.resize(_nRow,vector<double>(_nFactor,0.0));
      nullMatrix.resize(_nRow,vector<double>(_nFactor,0.0));
      break;}
  case 'P':
    { Proposal = _PAtomicdomain.getProposedAtoms();
      newMatrix.resize(_nFactor,vector<double>(_nCol,0.0));
      nullMatrix.resize(_nFactor,vector<double>(_nCol,0.0));
      break;}
  }
  */

  // read in the original _atomicProposal made from the prior
  unsigned long long location = _atomicProposal.begin()->first;
  double origMass = _atomicProposal.begin()->second;
  unsigned int bin;
  unsigned int iRow = _Row_changed[0];
  unsigned int iCol = _Col_changed[0];
  double delLL = 0;
  double delLLnew = 0;

  // ---------- consistency check
  if ((origMass < 0 && _oper_type =='B') || (origMass >=0 && _oper_type == 'D')) {
    cout << "Birth-death inconsistency!! origMass = " << origMass << 
      ", _oper_type = " << _oper_type << endl;
  }
  // -----------


  // ------------------------- DEATH -------------------------------------------
  if (_oper_type == 'D'){
 
    // put in the changes to the atomic space, and compute the corresponding change in 
    // the loglikelihood. 
    switch(the_matrix_label){
    case 'A':
      { 
	//newMatrix = atomicProposal2Matrix('A',AOrig);
        delLL = computeDeltaLL('A',D,S,AOrig,POrig,_nChange_matrixElemChange,_matrixElemChange); 
        _AAtomicdomain.acceptProposal(false); // "false" only means not to update _iter!
	// -----
	update_sysChi2(delLL);  // update system Chi2
	// -----
	//_AMatrix.matrix_update(newMatrix);
	_AMatrix.matrix_Elem_update(_matrixElemChange,_oper_type,_nChange_matrixElemChange);
	//FullnewMatrix = atomicProposal2FullMatrix('A',AOrig);
	break;}
    case 'P':
      { 
	//newMatrix = atomicProposal2Matrix('P',POrig);
	delLL = computeDeltaLL('P',D,S,AOrig,POrig,_nChange_matrixElemChange,_matrixElemChange);  
        _PAtomicdomain.acceptProposal(false); // "false" only means not to update _iter!
	// -----
	update_sysChi2(delLL);  // update system Chi2
	// -----
	//_PMatrix.matrix_update(newMatrix);
	_PMatrix.matrix_Elem_update(_matrixElemChange,_oper_type,_nChange_matrixElemChange);
	//FullnewMatrix = atomicProposal2FullMatrix('P',POrig);
	break;}
    } // end of switch-block
     
    // an attempt to rebirth
    attemptMass = -origMass;
    switch(the_matrix_label){
    case 'A':
      {
        if (!performUpdateKill('A',iRow, iCol, POrig)) { 
	  newMass = attemptMass;
	} else {
	  // ----
          //FullnewMatrix = atomicProposal2FullMatrix('A',AOrig);
	  //newMass = getMass('A',attemptMass,iRow,iCol,POrig,FullnewMatrix,D,S,rng);
	  newMass = getMass('A',attemptMass,iRow,iCol,POrig,AOrig,D,S,rng);
	  // ------- Q: think about it
	  if (newMass <= epsilon) {
	    newMass = attemptMass;
	  }
	} // end of if-block 

	_new_atomicProposal.insert(pair<unsigned long long,double>(location,newMass));
	extract_new_atomicProposal('A');
	_AAtomicdomain.setProposedAtomMass(_new_atomicProposal,true);  
        delLLnew = computeDeltaLL('A',D,S,AOrig,POrig,_new_nChange_matrixElemChange,_new_matrixElemChange);
	break;
      } // end of switch-block for A
    case 'P':
      {
	if (!performUpdateKill('P',iRow, iCol, AOrig)) { 
	  newMass = attemptMass;
	} else {
	  // -----
          //FullnewMatrix = atomicProposal2FullMatrix('P',POrig);
	  //newMass = getMass('P',attemptMass,iRow,iCol,AOrig,FullnewMatrix,D,S,rng);
          newMass = getMass('P',attemptMass,iRow,iCol,AOrig,POrig,D,S,rng);
	  // ----- Q: think about it
	  if (newMass <= epsilon) {
	    newMass = attemptMass;
	  }
	} // end of if-block 

	_new_atomicProposal.insert(pair<unsigned long long,double>(location,newMass));
	extract_new_atomicProposal('P');
	_PAtomicdomain.setProposedAtomMass(_new_atomicProposal,true);
        delLLnew = computeDeltaLL('P',D,S,AOrig,POrig,_new_nChange_matrixElemChange,_new_matrixElemChange);
	break;
      } // end of switch-block for P

    } // end of switch-block
    
    // -------- checking
    /*
      double delLLnew2;
      switch(the_matrix_label){
      case 'A':
      {delLLnew2 = computeDeltaLL2('A',D,S,AOrig,POrig,_new_nChange_matrixElemChange,_new_matrixElemChange);
      break;}
      case 'P':
      {delLLnew2 = computeDeltaLL2('P',D,S,AOrig,POrig,_new_nChange_matrixElemChange,_new_matrixElemChange);
      break;}
      }
      cout << "matrix_label = " << the_matrix_label << endl;
 
      cout << "delLLnew = " << delLLnew << endl;
      cout << "delLLnew2 = " << delLLnew2 << endl;
      unsigned int chkdelLL;
      if (fabs(delLLnew-delLLnew2) < 1.e-5){
      chkdelLL = 1;
      }else{
      chkdelLL = 0;
      }
      cout << "chkdelLL = " << chkdelLL << endl;

      cout << endl;
    */
    // ---------------



    // M-H sampling
    if (delLLnew*_annealingTemperature  < log(randgen('U',0,0))) {

      switch(the_matrix_label){
      case 'A':
	{ _AAtomicdomain.rejectProposal(false);
	  //return nullMatrix;
	  return false;
	  break;}
      case 'P':
	{ _PAtomicdomain.rejectProposal(false);
	  //return nullMatrix;
	  return false;
	  break;}
      } // end of switch-block

    } else {
      //newProposal.insert(pair<unsigned long long, double>(location, newMass)); 
      switch(the_matrix_label){
      case 'A':
	{ //_AAtomicdomain.setProposedAtomMass(newProposal, false);
          //newMatrix = atomicProposal2Matrix('A',AOrig);
          _AAtomicdomain.acceptProposal(false); 
	  // -----
	  update_sysChi2(delLLnew);  // update system Chi2
	  // -----
	  //return newMatrix;
	  return true;
          break;}
      case 'P':
	{ //_PAtomicdomain.setProposedAtomMass(newProposal, false);
          //newMatrix = atomicProposal2Matrix('P',POrig);
          _PAtomicdomain.acceptProposal(false); 
	  // -----	  
	  update_sysChi2(delLLnew);  // update system Chi2
	  // -----
	  //return newMatrix;
	  return true;
          break;}
      } // end of switch-block
    } // else of if-block for M-H sampling
    

    //return newMatrix;
    return false;    

  } // end of if-block for origMass < 0 (end of death process)
 
  
  // -------------- BIRTH ------------------------------------------------------
  if (_oper_type == 'B'){

    switch(the_matrix_label){
    case 'A':
      {
	//newMatrix = atomicProposal2Matrix('A',AOrig);
	//FullnewMatrix = atomicProposal2FullMatrix('A',AOrig);

	// checking conditions for update
	if (iRow >= _nRow || iCol >= _nFactor) {
	  throw logic_error("Cannot update pattern out of range in A.");
	}
	if ( !performUpdate('A',origMass, iRow, iCol, AOrig, POrig)) {
	  _AAtomicdomain.acceptProposal(false);
	  // -----------
          delLL = computeDeltaLL('A',D,S,AOrig,POrig,_nChange_matrixElemChange,_matrixElemChange);
	  update_sysChi2(delLL);  // update system Chi2
	  // ----------
	  //return newMatrix;
	  _new_atomicProposal.insert(pair<unsigned long long,double>(location,origMass));
    	  extract_new_atomicProposal('A');
	  return true;
	  // ------------
	}
	// ------------------------
	//FullnewMatrix[iRow][iCol] -= origMass;
	//newMass = getMass('A',origMass,iRow,iCol,POrig,FullnewMatrix,D,S,rng);
        newMass = getMass('A',origMass,iRow,iCol,POrig,AOrig,D,S,rng);
	_new_atomicProposal.insert(pair<unsigned long long,double>(location,newMass));
  	extract_new_atomicProposal('A');
	_AAtomicdomain.setProposedAtomMass(_new_atomicProposal,false);
        delLLnew = computeDeltaLL('A',D,S,AOrig,POrig,_new_nChange_matrixElemChange,_new_matrixElemChange);
	// ------------------------
	break;
      } // end of case-A block
    case 'P':
      {
	//newMatrix = atomicProposal2Matrix('P',POrig);
 	//FullnewMatrix = atomicProposal2FullMatrix('P',POrig);

	// checking conditions for update
	if (iRow >= _nFactor || iCol >= _nCol) {
	  throw logic_error("Cannot update pattern out of range in P.");
	}
	if ( !performUpdate('P',origMass, iRow, iCol, AOrig, POrig)) {
	  _PAtomicdomain.acceptProposal(false);
	  // -------
	  delLL = computeDeltaLL('P',D,S,AOrig,POrig,_nChange_matrixElemChange,_matrixElemChange);  
  	  update_sysChi2(delLL);  // update system Chi2
	  // -----------------
	  //return newMatrix;
	  _new_atomicProposal.insert(pair<unsigned long long,double>(location,origMass));
    	  extract_new_atomicProposal('P');
	  return true;
	  // -----------------
	}
	// ------------------------
	//FullnewMatrix[iRow][iCol] -= origMass;
	//newMass = getMass('P',origMass,iRow,iCol,AOrig,FullnewMatrix,D,S,rng);
        newMass = getMass('P',origMass,iRow,iCol,AOrig,POrig,D,S,rng);
	_new_atomicProposal.insert(pair<unsigned long long,double>(location,newMass));
  	extract_new_atomicProposal('P');
	_PAtomicdomain.setProposedAtomMass(_new_atomicProposal,false);
        delLLnew = computeDeltaLL('P',D,S,AOrig,POrig,_new_nChange_matrixElemChange,_new_matrixElemChange);
	// ------------------------
	break;
      } // end of case-P block
    } // end of switch-block

    // -------- checking
    /*
      double delLLnew2;
      switch(the_matrix_label){
      case 'A':
      {delLLnew2 = computeDeltaLL2('A',D,S,AOrig,POrig,_new_nChange_matrixElemChange,_new_matrixElemChange);
      break;}
      case 'P':
      {delLLnew2 = computeDeltaLL2('P',D,S,AOrig,POrig,_new_nChange_matrixElemChange,_new_matrixElemChange);
      break;}
      }
      cout << "matrix_label = " << the_matrix_label << endl;
 
      cout << "delLLnew = " << delLLnew << endl;
      cout << "delLLnew2 = " << delLLnew2 << endl;
      unsigned int chkdelLL;
      if (fabs(delLLnew-delLLnew2) < 1.e-5){
      chkdelLL = 1;
      }else{
      chkdelLL = 0;
      }
      cout << "chkdelLL = " << chkdelLL << endl;

      cout << endl;
    */
    // ---------------



    //newProposal.insert(pair<unsigned long long, double>(location, newMass)); 
    // This incorporates the modified change only, if any.
    switch(the_matrix_label){
    case 'A':
      { 
	//_AAtomicdomain.setProposedAtomMass(newProposal, false);
	//newMatrix = atomicProposal2Matrix('A',AOrig);
	_AAtomicdomain.acceptProposal(false);
	// -----
	update_sysChi2(delLLnew);  // update system Chi2
	// -----
	//return newMatrix;
	return true;
	break;}
    case 'P':
      { 
	//_PAtomicdomain.setProposedAtomMass(newProposal, false);
	//newMatrix = atomicProposal2Matrix('P',POrig);
	_PAtomicdomain.acceptProposal(false);
	// -----
	update_sysChi2(delLLnew);  // update system Chi2
	// -----
	//return newMatrix;
	return true;
	break;}
    }

    //return newMatrix;
    return false;

  } // end of if-block for origMass >=0
 

  //return newMatrix;
  return false;
}  // end of method birth_death

bool GibbsSampler::move_exchange(char the_matrix_label,
						    double const * const * D,
						    double const * const * S,
						    double ** AOrig,
						    double ** POrig)// in progress
{

  /*
  map<unsigned long long, double> newProposal;
  map<unsigned long long, double> Proposal;
  vector<vector<double> > newMatrix, nullMatrix;  // the thing to return in birth_death()

  switch(the_matrix_label){
  case 'A':
    { newMatrix.resize(_nRow,vector<double>(_nFactor,0.0));
      nullMatrix.resize(_nRow,vector<double>(_nFactor,0.0));
      break;}
  case 'P':
    { newMatrix.resize(_nFactor,vector<double>(_nCol,0.0));
      nullMatrix.resize(_nFactor,vector<double>(_nCol,0.0));
      break;}
  }
  */
 
  map<unsigned long long, double>::const_iterator atom;
  double chmass1,chmass2;       
  unsigned long long loc1, loc2;
  unsigned int bin1, bin2;
  double mass1, mass2;
  double newMass1, newMass2;
  atom = _atomicProposal.begin();
  chmass1 = atom->second;
  atom++;
  chmass2 = atom->second;

  // extract location, bin #, mass and changed mass corresponding to the 
  // atomic proposal such that "1" refers to a positive mass change and 
  // "2" a negative one.  
  switch(the_matrix_label){
  case 'A':
    {
      if (chmass1 > chmass2) {
	atom--;
	loc1 = atom->first;
	bin1 = _AAtomicdomain.getBin(loc1);
	mass1 = _AAtomicdomain.getMass(loc1);
	newMass1 = atom->second + mass1;
	atom++;
	loc2 = atom->first;
	bin2 = _AAtomicdomain.getBin(loc2);
	mass2 = _AAtomicdomain.getMass(loc2);
	newMass2 = atom->second + mass2;
      } else {
	loc1 = atom->first;
	bin1 = _AAtomicdomain.getBin(loc1);
	mass1 = _AAtomicdomain.getMass(loc1);
	newMass1 = atom->second + mass1;
	atom--;
	loc2 = atom->first;
	bin2 = _AAtomicdomain.getBin(loc2);
	mass2 = _AAtomicdomain.getMass(loc2);
	newMass2 = atom->second + mass2;
      }  // end of if-block for comparing chmass1 and chmass2
      break;} // end of case 'A' block
  case 'P':
    {
      if (chmass1 > chmass2) {
	atom--;
	loc1 = atom->first;
	bin1 = _PAtomicdomain.getBin(loc1);
	mass1 = _PAtomicdomain.getMass(loc1);
	newMass1 = atom->second + mass1;
	atom++;
	loc2 = atom->first;
	bin2 = _PAtomicdomain.getBin(loc2);
	mass2 = _PAtomicdomain.getMass(loc2);
	newMass2 = atom->second + mass2;
      } else {
	loc1 = atom->first;
	bin1 = _PAtomicdomain.getBin(loc1);
	mass1 = _PAtomicdomain.getMass(loc1);
	newMass1 = atom->second + mass1;
	atom--;
	loc2 = atom->first;
	bin2 = _PAtomicdomain.getBin(loc2);
	mass2 = _PAtomicdomain.getMass(loc2);
	newMass2 = atom->second + mass2;
      }  // end of if-block for comparing chmass1 and chmass2
      break;}  // end of case 'P' block
  } // end of switch-block for extracting the atomic proposal info

  // return nullMatrix if bin1 == bin2
  if (bin1 == bin2){
    // cout << "Exchanges in the same bin!! So, no change!" << endl;
    //return nullMatrix;
    return false;
  }


  // preparing quantities for possible Gibbs computation later.
  bool exchange = false;
  double priorLL = 0.;

  unsigned int jGene, jSample, jPattern;
  bool anyNonzero = false;
  bool useGibbs = true;


  unsigned int iGene1, iPattern1, iGene2, iPattern2, iSample1, iSample2;
  switch(the_matrix_label){
  case 'A':
    {
      iGene1 = getRow('A',bin1);
      iPattern1 = getCol('A',bin1);
      iGene2 = getRow('A',bin2);
      iPattern2 = getCol('A',bin2);
      break;}
  case 'P':
    {
      iPattern1 = getRow('P',bin1);
      iSample1 = getCol('P',bin1);
      iPattern2 = getRow('P',bin2);
      iSample2 = getCol('P',bin2);
      break;}
  }

  // ---------------------------------------------
  switch(the_matrix_label){
  case 'A':
    {
      for (jSample = 0; jSample < _nCol; jSample++) {
	if (POrig[iPattern1][jSample] > epsilon) {
	  anyNonzero = true;
	  break;
	}
	if (POrig[iPattern2][jSample] > epsilon) {
	  anyNonzero = true;
	  break;
	}
      }  // end of for-block to determine the existence of corresponding 
      // non-zero elements in P
      if (!anyNonzero)  {  // cannot update in Gibbs way
	useGibbs = false;
      }
      break;}
  case 'P':
    {
      for (jGene = 0; jGene < _nRow; jGene++) {
	if (AOrig[jGene][iPattern1] > epsilon) {
	  anyNonzero = true;
	  break;
	}
	if (AOrig[jGene][iPattern2] > epsilon) {
	  anyNonzero = true;
	  break;
	}
      }  // end of for-block to determine the existence of corresponding 
      // non-zero elements in P
      if (!anyNonzero)  {  // cannot update in Gibbs way
	useGibbs = false;
      }
      break;}
  }  // end of switch-block

  // -------------------------------------------------------------------------
  // EXCHANGE ACTION when initial useGibbs = true and check if Gibbs is usable.

  double s = 0.0;
  double su = 0.0;
  double mock = 0.0;
  double mock1 = 0.0;
  double mock2 = 0.0;
  double mean = 0; 
  double sd = 0;
  double gibbsMass1, gibbsMass2;

  if (useGibbs == true){

    switch(the_matrix_label){
    case 'A':{
      // ---------- EXCHANGE ACTION WITH A ----------------------------
      if (_AAtomicdomain.inDomain(loc1) && _AAtomicdomain.inDomain(loc2)) 

	{
	  exchange = true;

	  double Aeff;
	  // compute the distribution parameters
	  if (iGene1 == iGene2) {
	    for (jSample = 0; jSample < _nCol; jSample++) {
	      // Calculate the mock term
	      mock = D[iGene1][jSample];
	      for (jPattern = 0; jPattern < _nFactor; jPattern++) {
		Aeff = AOrig[iGene1][jPattern];
		mock -= Aeff * POrig[jPattern][jSample];
	      } // end of for-block that make changes to elements in A

	      s += ( ( POrig[iPattern1][jSample]-POrig[iPattern2][jSample] ) *
		     ( POrig[iPattern1][jSample]-POrig[iPattern2][jSample] ) ) /
		pow(S[iGene1][jSample],2);
	      su += mock *( POrig[iPattern1][jSample]-POrig[iPattern2][jSample] ) /
		pow( S[iGene1][jSample],2);
	    }
	  }  // end of case iGene1 = iGene2
	  else  { // iGene1 != iGene2 
	    for (jSample = 0; jSample < _nCol; jSample++) {
	      mock1 = D[iGene1][jSample];
	      mock2 = D[iGene2][jSample]; 

	      for (jPattern = 0; jPattern < _nFactor; jPattern++) {
		Aeff = AOrig[iGene1][jPattern];
		mock1 -= Aeff * POrig[jPattern][jSample];
		Aeff = AOrig[iGene2][jPattern];
		mock2 -= Aeff * POrig[jPattern][jSample];
	      } // end of for-block for adding changes to A

	      s  += pow(POrig[iPattern1][jSample] / S[iGene1][jSample],2) +
		pow(POrig[iPattern2][jSample] / S[iGene2][jSample],2);
	      su += mock1*POrig[iPattern1][jSample]/pow(S[iGene1][jSample],2) -
		mock2*POrig[iPattern2][jSample]/pow(S[iGene2][jSample],2);
	    }
	  } // end of if-block for calculating s and su for case 'A'

	  s = s * _annealingTemperature;
	  su = su * _annealingTemperature; 
	  mean = su / s;
	  sd = 1./sqrt(s);
	  // end of compute distribution parameters for A	       

	}  // end of if-block for checking whether the changes are in domain (the exchange block)
	
      break;} // end of switch block for EXCHANGE ACTION with A
      // ---------- EXCHANGE ACTION WITH P ----------------------------
    case 'P': {
      if (_PAtomicdomain.inDomain(loc1) && _PAtomicdomain.inDomain(loc2)) 

	{
	  exchange = true;

	  double Peff;
	  // compute the distribution parameters
	  if (iSample1 == iSample2) {
	    for (jGene = 0; jGene < _nRow; jGene++) {
	      // Calculate the mock term
	      mock = D[jGene][iSample1];
	      for (jPattern = 0; jPattern < _nFactor; jPattern++) {
		Peff = POrig[jPattern][iSample1];
		mock -=  AOrig[jGene][jPattern]*Peff;
	      } // end of for-block that make changes to elements in P

	      s += pow( ((AOrig[jGene][iPattern1]-AOrig[jGene][iPattern2])/
			 S[jGene][iSample1]),2);
	      su += mock *( AOrig[jGene][iPattern1]-AOrig[jGene][iPattern2] ) /
		pow( S[jGene][iSample1],2);
	    }
	  }  // end of case iSample1 = iSample2
	  else  { // iSample1 != iSample2 
	    for (jGene = 0; jGene < _nRow; jGene++) {
	      mock1 = D[jGene][iSample1];
	      mock2 = D[jGene][iSample2]; 

	      for (jPattern = 0; jPattern < _nFactor; jPattern++) {
		Peff = POrig[jPattern][iSample1];
		mock1 -= AOrig[jGene][jPattern]*Peff;
		Peff = POrig[jPattern][iSample2];
		mock2 -= AOrig[jGene][jPattern]*Peff;
	      } // end of for-block for adding changes to P

	      s  += pow(AOrig[jGene][iPattern1] / S[jGene][iSample1],2) +
		pow(AOrig[jGene][iPattern2] / S[jGene][iSample2],2);
	      su += mock1*AOrig[jGene][iPattern1]/pow(S[jGene][iSample1],2) -
		mock2*AOrig[jGene][iPattern2]/pow(S[jGene][iSample2],2);
	    }
	  } // end of if-block for calculating s and su for case 'P'

	  s = s * _annealingTemperature;
	  su = su * _annealingTemperature; 
	  mean = su / s;
	  sd = 1./sqrt(s);
	  // end of compute distribution parameters for P	       

	}  // end of if-block for checking whether the changes are in domain (the exchange block)
	
      break;} // end of switch block for EXCHANGE ACTION with P
 
    } // end of switch block for EXCHANGE ACTION

  } // end of if-block for operations with possibly Gibbs sampling


  if (s == 0. && su == 0.){
    useGibbs = false;
    //cout << "Parameters aren't updated -> useGibbs = false, do M-H" << endl;
  }

  // -------------------------------------------------------------------------
 
  if(useGibbs == true){

    // set newMass1	 
    // need to retain exponential prior
    //double mean = (2.*_annealingTemperature*su - lambda) / 
    double plower = sub_func::pnorm(-mass1, mean, sd, DOUBLE_NEGINF, 0);
    double pupper = sub_func::pnorm(mass2, mean, sd, DOUBLE_NEGINF, 0);
    double u = plower + randgen('U',0,0)*(pupper - plower);
 
    // must sample from prior if the computed parameters are not good for Gibbs
    if (plower >  0.95 || 
	pupper < 0.05 ||
	s < epsilon || 
	newMass1 == DOUBLE_POSINF ||
	newMass1 == DOUBLE_NEGINF) {
      // do not make a change
      useGibbs = false;

    }

    if (useGibbs == true){

      gibbsMass1 = sub_func::qnorm(u, mean, sd, DOUBLE_NEGINF, 0);
      if (gibbsMass1 < -mass1) gibbsMass1 = -mass1;
      if (gibbsMass1 > mass2) gibbsMass1 = mass2;
      gibbsMass2 = - gibbsMass1;

      // update new masses
      double delLLnew;
      _new_nChange_matrixElemChange = 2;
      _new_atomicProposal.insert(pair<unsigned long long, double>(loc1,gibbsMass1));
      _new_atomicProposal.insert(pair<unsigned long long, double>(loc2,gibbsMass2));
      switch(the_matrix_label){
      case 'A':
	{
	  extract_new_atomicProposal('A');
	  delLLnew = computeDeltaLL('A',D,S,AOrig,POrig,_new_nChange_matrixElemChange,_new_matrixElemChange);
	  _AAtomicdomain.setProposedAtomMass(_new_atomicProposal, false);
	  //newMatrix = atomicProposal2Matrix('A',AOrig);
	  _AAtomicdomain.acceptProposal(false); 
  	  // -----
	  update_sysChi2(delLLnew);  // update system Chi2
	  // -----
	  break;}
      case 'P':
	{
	  extract_new_atomicProposal('P');
	  delLLnew = computeDeltaLL('P',D,S,AOrig,POrig,_new_nChange_matrixElemChange,_new_matrixElemChange);
	  _PAtomicdomain.setProposedAtomMass(_new_atomicProposal, false);
	  //newMatrix = atomicProposal2Matrix('P',POrig);
	  _PAtomicdomain.acceptProposal(false);
 	  // -----
	  update_sysChi2(delLLnew);  // update system Chi2
	  // -----
	  break;}
      }  // end of switch-block       
      // return newMatrix;
      /*
      // ---- checking Gibbs in exchange
      cout << " ---------------------------------- " << endl;
      cout << " Check Gibbs Exchange: " << endl;
      cout << "s = " << s << ", su = " << su << ", mean = " << mean << ", sd = " << sd << endl;
      cout << "Exchange with matrix: " << the_matrix_label << ", gibbsMass1 = " << gibbsMass1 
           << ", gibbsMass2 = " << gibbsMass2 
	   << ", delLLnew = " << delLLnew << endl;
      //if(fabs(delLLnew) > 500.0 ){
	for (unsigned int m = 0; m < 2; ++m ){
	  cout << "chRow[" << m << "] = " << _new_Row_changed[m] << ", chCol[" << m << "] = " 
               << _new_Col_changed[m] << ", mass changed = " << _new_mass_changed[m] << endl; 
	}
	display_matrix('A');
	display_matrix('P');
        if (fabs(delLLnew) > 500.0 ){
            cout << "delLLNew bigger than 500! " << endl;
        }

	//}
	*/
      // ----------------
      return true;


    } // end of inner if-block for final updating with Gibbs

  } // end of outer if-block for useGibbs == true && s > epsilon


  // ----------------------------
  // Metropolis-Hasting

  double pold = 0.;
  double pnew = 0.;
  double lambda;
  switch(the_matrix_label){
  case 'A':
    {lambda = _lambdaA;
      break;}
  case 'P':
    {lambda = _lambdaP;
      break;}
  }
  
  if (_oper_type == 'E'){
    if (mass1 > mass2) {
      pnew = sub_func::dgamma(newMass1, 2., 1./lambda, false);
      if (newMass1 > newMass2) {
	pold = sub_func::dgamma(mass1, 2., 1./lambda, false);
      } else {
	pold = sub_func::dgamma(mass2, 2., 1./lambda, false);
      }
    } else {
      pnew = sub_func::dgamma(newMass2, 2., 1./lambda, false);
      if (newMass1 > newMass2) {
	pold = sub_func::dgamma(mass1, 2., 1./lambda,  false);
      } else {
	pold = sub_func::dgamma(mass2, 2., 1./lambda, false);
      }
    }  
  } // end of if-block for Exchange

  /*
  if(_oper_type == 'M'){
    
      if (mass1 > mass2) {
      pnew = newMass1;
      if (newMass1 > newMass2) {
      pold = mass1;
      } else {
      pold = mass2;
      }
      } else {
      pnew = newMass2;
      if (newMass1 > newMass2) {
      pold = mass1;
      } else {
      pold = mass2;
      }
      } 
   
  } // end of if-block for move
  */

  if (pnew == 0. && pold == 0.) {
    priorLL = 0.0;
  } else if(pnew != 0. && pold == 0.) {
    priorLL = DOUBLE_POSINF;
  } else {
    priorLL = log(pnew / pold);
  }
	 
  double delLLnew;
  switch(the_matrix_label){
  case 'A':
    {delLLnew = computeDeltaLL('A',D,S,AOrig,POrig,_nChange_matrixElemChange,_matrixElemChange);
      break;}
  case 'P':
    {delLLnew = computeDeltaLL('P',D,S,AOrig,POrig,_nChange_matrixElemChange,_matrixElemChange);
      break;}
  }

  // -------- checking
  /*
    double delLLnew2;
    switch(the_matrix_label){
    case 'A':
    {delLLnew2 = computeDeltaLL2('A',D,S,AOrig,POrig,_nChange_matrixElemChange,_matrixElemChange);
    break;}
    case 'P':
    {delLLnew2 = computeDeltaLL2('P',D,S,AOrig,POrig,_nChange_matrixElemChange,_matrixElemChange);
    break;}
    }
    cout << "matrix_label = " << the_matrix_label << endl;
    cout << "_oper_type = " << _oper_type << endl;
    cout << "delLLnew = " << delLLnew << endl;
    cout << "delLLnew2 = " << delLLnew2 << endl;
    unsigned int chkdelLL;
    if (fabs(delLLnew-delLLnew2) < 1.e-5){
    chkdelLL = 1;
    }else{
    chkdelLL = 0;
    }
    cout << "chkdelLL = " << chkdelLL << endl;

    switch(the_matrix_label){
    case 'A':
    {
    cout << "iGene1 = " << iGene1 << endl;
    cout << "iGene2 = " << iGene2 << endl;
    break;
    }
    case 'P':
    {
    cout << "iSample1 = " << iSample1 << endl;
    cout << "iSample2 = " << iSample2 << endl;
    break;
    }
    }

    cout << endl;
  */
  // ---------------

  double totalLL = priorLL + delLLnew * _annealingTemperature;

  _new_nChange_matrixElemChange = 2;
  _new_atomicProposal.insert(pair<unsigned long long, double>(loc1,newMass1-mass1));
  _new_atomicProposal.insert(pair<unsigned long long, double>(loc2,newMass2-mass2));
  switch(the_matrix_label){
  case 'A': {
    //newMatrix = atomicProposal2Matrix('A',AOrig);
    extract_new_atomicProposal('A');
    break;
  }
  case 'P': {
    //newMatrix = atomicProposal2Matrix('P',POrig);
    extract_new_atomicProposal('P');
    break;
  }
  }

  double tmp;
  if (priorLL == DOUBLE_POSINF){
    //return newMatrix;
    return true;
  } else {
    tmp = priorLL + delLLnew*_annealingTemperature;
  }

  double rng = log(randgen('U',0,0));

  if (_oper_type == 'E'){
    if (tmp  < rng) {
      switch(the_matrix_label){
      case 'A':
	{ _AAtomicdomain.rejectProposal(false);
	  //return nullMatrix;
	  return false;
	  break;}
      case 'P':
	{ _PAtomicdomain.rejectProposal(false);
	  //return nullMatrix;
	  return false;
	  break;}
      } // end of switch-block
    } else {
 
      switch(the_matrix_label){
      case 'A':
	{ 	  
          _AAtomicdomain.acceptProposal(false); 
	  // -----
	  update_sysChi2(delLLnew);  // update system Chi2
	  // -----
	  //return newMatrix;
	  return true;
          break;
        }
      case 'P':
	{           
          _PAtomicdomain.acceptProposal(false); 
	   // -----
	   update_sysChi2(delLLnew);  // update system Chi2
	   // -----
	  //return newMatrix;
	  return true;
          break;
        }
      } // end of switch-block       
    }  
  } // end of the M-H determination block for _oper_type = E

 
  if (_oper_type == 'M'){
    if (tmp < rng) {
      switch(the_matrix_label){
      case 'A':
	{ _AAtomicdomain.rejectProposal(false);
	  //return nullMatrix;
	  return false;
	  break;}
      case 'P':
	{ _PAtomicdomain.rejectProposal(false);
	  //return nullMatrix;
	  return false;
	  break;}
      } // end of switch-block
    } else { 
      switch(the_matrix_label){
      case 'A':
	{ 	  
	  _AAtomicdomain.acceptProposal(false); 
	  // -----
	  update_sysChi2(delLLnew);  // update system Chi2
	  // -----
	  //return newMatrix;
	  return true;
	  break;
	}
      case 'P':
	{           
	  _PAtomicdomain.acceptProposal(false); 
	  // -----
	  update_sysChi2(delLLnew);  // update system Chi2
	  // -----
	  //return newMatrix;
	  return true;
	  break;
	}
      } // end of switch-block       
    }  
  } // end of the M-H determination block for _oper_type = M
 

  // end of M-H sampling
 
  // ----------

  //return nullMatrix;
  return false;

} // end of method move_exchange


// ************ METHODS FOR LOOPING AND CONTROL *****************************
void GibbsSampler::set_iter(unsigned long ext_iter){
  _iter = ext_iter;
} 

double GibbsSampler::get_AnnealingTemperature(){
  return _annealingTemperature;
}


// -----------------------------------------------------------------------------
// Calculate the _new annealingTemperature as a function of iteration step _iter.
// (Note: the annealingTemperature here is really the inverted temperature!)
void GibbsSampler::set_AnnealingTemperature() { 
  double SASteps = _nEquil;
  double SATemp = ( (double) _iter + 1. ) / (SASteps / 2.);

  if (SATemp > 1.) SATemp = 1;
  if (SATemp < 0) {
    throw logic_error("Invalid annealing temperature.");
  }

  _annealingTemperature = SATemp;
}


// -----------------------------------------------------------------------------
// 
void GibbsSampler::check_atomic_matrix_consistency(char the_matrix_label)
{
  double total_atom_mass;
  double total_matrix_mass;
  
  switch(the_matrix_label){
  case 'A': {
    total_atom_mass = _AAtomicdomain.get_atomicDomain_totalmass();
    total_matrix_mass = _AMatrix.cal_totalsum();
    break;
  } 
  case 'P': {
    total_atom_mass = _PAtomicdomain.get_atomicDomain_totalmass();
    total_matrix_mass = _PMatrix.cal_totalsum();
    break;
  }
  } // end of switch-block

  double diff_total_mass = fabs(total_atom_mass - total_matrix_mass);

  if(diff_total_mass > 1.e-5){
    cout << "Mass inconsistency!! Total mass difference = " << diff_total_mass << endl;
    cout << "total atom mass = " << total_atom_mass << endl;
    cout << "total matrix mass = " << total_matrix_mass << endl;
    cout << "Oper_type = " << _oper_type << endl;
    throw logic_error("Mass inconsistency between atomic domain and matrix!");
  } //else {
    //cout << "Things beautiful! Total mass difference = " << diff_total_mass << endl << endl;
    //} 
  
}

// -----------------------------------------------------------------------------
// Here we form the matrices _Amean, _Asd, _Pmean, _Psd. They are temporary 
// matrices to be modified in compute_statistics() to calculate the means and the variances
// of individual matrix entries. Note tha here we have cumulated to each matrix element all its
// realizations in time as:
// _Amean_{ij} = \sum_{t} A_{ij}*k_{j} 
// _Asd_{ij} = \sum_{t}  (A_{ij}*k_{j} )^2
// _Pmean_{ij} = \sum_{t} P_{ij}/k_{i} 
// _Psd_{ij} = \sum_{t}  (P_{ij}/k_{i} )^2

void GibbsSampler::compute_statistics_prepare_matrices(unsigned long statindx){

  double ** A = _AMatrix.get_matrix();
  double ** P = _PMatrix.get_matrix();
  vector<double> k(_nFactor);  // normalized vector

  // compute the normalization vector
  for (int m=0; m< _nFactor; ++m){
    k[m] = 0.;
    for (int n=0; n< _nCol; ++n){
      k[m] += P[m][n];
    }
    if (k[m] == 0){  // when the whole row of P is zero, then don't do anything
      k[m] = 1.0;
    }
  }


  // construct the mean and var matrices at statindx = 1
  if (statindx == 1){
        
    _Amean = new double * [_nRow];
    for (int m=0; m < _nRow ; ++m)
      {_Amean[m] = new double [_nFactor];}
    for (int m=0; m < _nRow; ++m){
      for (int n=0; n < _nFactor; ++n){
	_Amean[m][n] = A[m][n] * k[n];
      }
    }

    _Asd = new double * [_nRow];
    for (int m=0; m < _nRow ; ++m)
      {_Asd[m] = new double [_nFactor];}
    for (int m=0; m < _nRow; ++m){
      for (int n=0; n < _nFactor; ++n){
	_Asd[m][n] = pow(A[m][n]*k[n],2);
      }
    } 

    _Pmean = new double * [_nFactor];
    for (int m=0; m < _nFactor ; ++m)
      {_Pmean[m] = new double [_nCol];}
    for (int m=0; m < _nFactor; ++m){
      for (int n=0; n < _nCol; ++n){
	_Pmean[m][n] = P[m][n] / k[m];
      }
    } 

    _Psd = new double * [_nFactor];
    for (int m=0; m < _nFactor ; ++m)
      {_Psd[m] = new double [_nCol];}
    for (int m=0; m < _nFactor; ++m){
      for (int n=0; n < _nCol; ++n){
	_Psd[m][n] = pow(P[m][n] / k[m] , 2);
      }
    } 

  } // end of if-block for matrix construction statindx == 1

  // increment the mean and var matrices at statindx != 1
  if (statindx > 1){

    for (int m=0; m < _nRow; ++m){
      for (int n=0; n < _nFactor; ++n){
	_Amean[m][n] += A[m][n] * k[n];
      }
    }

    for (int m=0; m < _nRow; ++m){
      for (int n=0; n < _nFactor; ++n){
	_Asd[m][n] += pow(A[m][n]*k[n],2);
      }
    }

    for (int m=0; m < _nFactor; ++m){
      for (int n=0; n < _nCol; ++n){
	_Pmean[m][n] += P[m][n] / k[m];
      }
    }

    for (int m=0; m < _nFactor; ++m){
      for (int n=0; n < _nCol; ++n){
	_Psd[m][n] += pow(P[m][n] / k[m],2);
      }
    }

  } // end of if-block for matrix incrementation statindx > 1

  /*
    cout << "statinx = " << statindx << endl;
    cout << " A = " << endl;
    local_display_matrix2(A,_nRow,_nFactor);
    cout << " Amean = " << endl;
    local_display_matrix2(_Amean,_nRow,_nFactor);
    cout << " Asd = " << endl;
    local_display_matrix2(_Asd,_nRow,_nFactor);
    for (int m=0; m < _nFactor ; ++m){
    cout << "k" << m << " = " << k[m] << endl;
    }
    cout << " P = " << endl;
    local_display_matrix2(P,_nFactor,_nCol);
    cout << " Pmean = " << endl;
    local_display_matrix2(_Pmean,_nFactor,_nCol);
    cout << " Psd = " << endl;
    local_display_matrix2(_Psd,_nFactor,_nCol);
  */  

} // end of method compute_statistics_prepare_matrices

// -----------------------------------------------------------------------------
// Here, we use the matrices prepared in compute_statistics_prepare_matrices()
// to compute the means and variances of individual elements in A and P.
void GibbsSampler::compute_statistics(char outputFilename[],
                                      char outputAmean_Filename[],char outputAsd_Filename[],
                                      char outputPmean_Filename[],char outputPsd_Filename[],
				      char outputAPmean_Filename[],
                                      unsigned int Nstat){

  // compute statistics for A

  for (int m=0; m < _nRow ; ++m){
    for (int n=0; n<_nFactor; ++n){
      _Amean[m][n] = _Amean[m][n] / Nstat;
    }
  }

  for (int m=0; m < _nRow ; ++m){
    for (int n=0; n<_nFactor; ++n){
      _Asd[m][n] = sqrt((_Asd[m][n] - Nstat*pow(_Amean[m][n],2)) / (Nstat-1));
    }
  }

  // compute statistics for P
  for (int m=0; m < _nFactor ; ++m){
    for (int n=0; n<_nCol; ++n){
      _Pmean[m][n] = _Pmean[m][n] / Nstat;
    }
  }

  for (int m=0; m < _nFactor ; ++m){
    for (int n=0; n<_nCol; ++n){
      _Psd[m][n] = sqrt((_Psd[m][n] - Nstat*pow(_Pmean[m][n],2)) / (Nstat-1));
    }
  }


  double ** APmean;
  APmean = new double * [_nRow];
  for (int m=0; m < _nRow ; ++m)
    {APmean[m] = new double [_nCol];}
  for (unsigned int m=0; m < _nRow; ++m){
    for (unsigned int n=0; n < _nCol; ++n){
      APmean[m][n]=0.0;
      for (unsigned int k=0; k < _nFactor; ++k){
	APmean[m][n] += _Amean[m][k]*_Pmean[k][n];
      }
    }
  }



  ofstream outputFile;
  outputFile.open(outputFilename,ios::out|ios::app);

  outputFile << " ************************************************* " << endl;
  outputFile << " ---------- OUTPUT FINAL STATISTICS -------------- " << endl;
  outputFile << " ************************************************* " << endl;

  outputFile << " Number of samples for computing statistics = " << Nstat << endl;
  outputFile << " Amean = " << endl << endl;
  local_display_matrix2F(outputFile,_Amean,_nRow,_nFactor);
  outputFile << endl;
  outputFile << " Asd = " << endl << endl;
  local_display_matrix2F(outputFile,_Asd,_nRow,_nFactor);
  outputFile << endl;
  outputFile << " Pmean = " << endl << endl;
  local_display_matrix2F(outputFile,_Pmean,_nFactor,_nCol);
  outputFile << endl;
  outputFile << " Psd = " << endl << endl;
  local_display_matrix2F(outputFile,_Psd,_nFactor,_nCol);
  outputFile << endl;
  outputFile << "The product Amean*Pmean gives " << endl << endl;
  local_display_matrix2F(outputFile,APmean,_nRow,_nCol);
  outputFile << endl;

  outputFile.close();

  // independent output of Amean, Asd, Pmean, Psd, APmean
  ofstream output_Amean;
  output_Amean.open(outputAmean_Filename,ios::out);
  local_display_matrix2F(output_Amean,_Amean,_nRow,_nFactor);
  output_Amean.close();

  ofstream output_Asd;
  output_Asd.open(outputAsd_Filename,ios::out);
  local_display_matrix2F(output_Asd,_Asd,_nRow,_nFactor);
  output_Asd.close();

  ofstream output_Pmean;
  output_Pmean.open(outputPmean_Filename,ios::out);
  local_display_matrix2F(output_Pmean,_Pmean,_nFactor,_nCol);
  output_Pmean.close();

  ofstream output_Psd;
  output_Psd.open(outputPsd_Filename,ios::out);
  local_display_matrix2F(output_Psd,_Psd,_nFactor,_nCol);
  output_Psd.close();

  ofstream output_APmean;
  output_APmean.open(outputAPmean_Filename,ios::out);
  local_display_matrix2F(output_APmean,APmean,_nRow,_nCol);
  output_APmean.close();



}



// *****************************************************************************
// Adaptation from the original code:

// -----------------------------------------------------------------------------
bool GibbsSampler::performUpdate(char the_matrix_label, double origMass, 
				 unsigned int iRow, unsigned int iCol,
				 double const * const * A, double const * const * P) {

  unsigned int otherDim;
  bool nonZero = false;
  double llike;

  switch(the_matrix_label){

    // check that there is mass for this pattern in the P Matrix
  case 'A': 
    {
      for (otherDim = 0; otherDim < _nCol; otherDim++) {
	if (P[iCol][otherDim] > epsilon) {
	  nonZero = true;
	}
      }
      if (!nonZero) return false;
      if (origMass < 0) {
	throw logic_error("Should not be checking death during birth update.");
      }
      break;
    } // end of case 'A'

    // check that there is mass for this pattern in the A Matrix
  case 'P': 
    {
      for (otherDim = 0; otherDim < _nRow; otherDim++) {
	if (A[otherDim][iRow] > epsilon) {
	  nonZero = true;
	}
      }
      if (!nonZero) return false;      
      if (origMass < 0) {
	throw logic_error("Should not be checking death during birth update.");
      }
      break;
    } // end of case 'P'
  } // end of switch

  return true;
} // end of performUpdate

// -----------------------------------------------------------------------------
bool GibbsSampler::performUpdateKill(char the_matrix_label, unsigned int iRow, unsigned int iCol,
				     double const * const * otherMatrix){

  unsigned int otherDim;
  bool nonZero = false;

  switch(the_matrix_label){


    // check that there is mass for this pattern in the P Matrix
  case 'A': 
    {
      for (otherDim = 0; otherDim < _nCol; otherDim++) {
	if (otherMatrix[iCol][otherDim] > epsilon) {
	  nonZero = true;
	}
      }

      if (!nonZero) return false;
      break;
    }

    // check that there is mass for this pattern in the A Matrix
  case 'P': 
    {
      for (otherDim = 0; otherDim < _nRow; otherDim++) {
	if (otherMatrix[otherDim][iRow] > epsilon) {
	  nonZero = true;
	}
      }
      if (!nonZero) return false;

      break;
    }
  } // end of switch block

  return true; 
}

// -----------------------------------------------------------------------------
/*
double GibbsSampler::getMass(char the_matrix_label, double origMass,
			     unsigned int iRow,
			     unsigned int iCol,
			     double const * const * otherMatrix, 
			     const vector<vector<double> > & currentChainMatrix,
			     double const * const * D, double const * const * S,
			     double rng)
{
  double DOUBLE_POSINF = numeric_limits<double>::max();
  unsigned int iGene, iPattern, iSample, jPattern;
  double lambda;

  switch(the_matrix_label)
    {
    case 'A': 
      {
	iGene = iRow;
	iPattern = iCol;
	lambda = _AAtomicdomain.getLambda();
	break;
      }
    case 'P': 
      {
	iPattern = iRow;
	iSample = iCol;
	lambda = _PAtomicdomain.getLambda();
	break;
      }
    }

  bool anyNonzero = false;

  // determine the parameters for finding the mass
  double s  = 0.;
  double su = 0.;
  double mock;

  switch(the_matrix_label){

  case 'A': 
    {
      double Aeff;
      for (iSample = 0; iSample < _nCol; iSample++) 
	{
	  mock = D[iGene][iSample];
	  for (jPattern = 0; jPattern <_nFactor; jPattern++) 
	    {
	      Aeff = currentChainMatrix[iGene][jPattern];
	      if (jPattern == iPattern && origMass < 0) 
		{
		  Aeff -= origMass;
		}
	      mock -= Aeff * otherMatrix[jPattern][iSample];
	      if (otherMatrix[jPattern][iSample] > epsilon) 
		{
		  anyNonzero = true;
		}
	    }
	  s += _annealingTemperature * 
	    (otherMatrix[iPattern][iSample] * 
	     otherMatrix[iPattern][iSample]) /
	    (2*S[iGene][iSample]*S[iGene][iSample]);
	  su += _annealingTemperature * 
	    otherMatrix[iPattern][iSample]*mock /
	    (2*S[iGene][iSample]*S[iGene][iSample]);	    
	}
      break;
    }

  case 'P': 
    {
      double Peff;
      for (iGene = 0; iGene < _nRow; iGene++) 
	{
	  mock = D[iGene][iSample];
	  for (jPattern = 0; jPattern < _nFactor; jPattern++) 
	    {
	      Peff = currentChainMatrix[jPattern][iSample];
	      if (jPattern == iPattern && origMass < 0) 
		{
		  Peff -= origMass;
		}
	      mock -= otherMatrix[iGene][jPattern] * Peff;

	      if (otherMatrix[iGene][jPattern] > epsilon) 
		{
		  anyNonzero= true;
		}
	    }
	  s += _annealingTemperature * 
	    (otherMatrix[iGene][iPattern]*otherMatrix[iGene][iPattern]) / 
	    (2*S[iGene][iSample]*S[iGene][iSample]);
	  su += _annealingTemperature * otherMatrix[iGene][iPattern]*mock / 
	    (2*S[iGene][iSample]*S[iGene][iSample]);
	}
      break;
    }

  }

    
  double mean  = (2*su - lambda)/(2*s);
  double sd = 1./sqrt(2*s);
    
  // note: is bounded below by zero so have to use inverse sampling!
  // based upon algorithm in DistScalarRmath.cc (scalarRandomSample)

  double plower = sub_func::pnorm(0., mean, sd, DOUBLE_NEGINF, 0);
  double pupper = 1.;
  double u = plower + randgen('U',0,0) * (pupper - plower);
  // -------------------------------------------------------------
  // this line seems to be misplaced.
  // double newMass = sub_func::qnorm(u, mean, sd, DOUBLE_NEGINF, 0);
  double newMass = 0;
  // ------------------

  // if the likelihood is flat and nonzero, 
  // force to sample strictly from the prior
  if ( plower == 1 || s < 1.e-5 || newMass == DOUBLE_POSINF || pupper == 0) { 
    if (origMass < 0) {    // death case
      newMass = abs(origMass);
    } else {
      newMass = 0.;  // birth case
    }
  } // end of first comparison
  else if (plower >= 0.99) {
    double tmp1 = sub_func::dnorm(0, mean, sd, false); 
    double tmp2 = sub_func::dnorm(10*lambda, mean, sd, false);
    if ( (tmp1 > epsilon) && (fabs(tmp1-tmp2) < epsilon) )   {
      if (origMass < 0) {   // death case
	return 0.;
      }
      return origMass;   // birth case
    }  
  } // end of second comparison
  else {
    newMass = sub_func::qnorm(u, mean, sd, DOUBLE_NEGINF, 0);  // both death and birth
  }  // end of if-block for the remaining case


  // limit the mass range
    switch(the_matrix_label) {
    case 'A':
      {
	if (newMass > _max_gibbsmassA)
	  newMass = _max_gibbsmassA;
	break;
      }
    case 'P':
      {
	if (newMass > _max_gibbsmassP)
	  newMass = _max_gibbsmassP;
	break;
      }
    }

  if (newMass < 0) newMass = 0;   // due to the requirement that newMass > 0 

  return newMass;

}
*/


double GibbsSampler::getMass(char the_matrix_label, double origMass,
			     unsigned int iRow,
			     unsigned int iCol,
			     double const * const * otherMatrix, 
			     double const * const * currentChainMatrix,
			     double const * const * D, double const * const * S,
			     double rng)
{
  double DOUBLE_POSINF = numeric_limits<double>::max();
  unsigned int iGene, iPattern, iSample, jPattern;
  double lambda;

  // ---- check
  //cout << "Inside getMass, _oper_type = " << _oper_type << ", origMass = " << origMass << endl;
  // ------------

  switch(the_matrix_label)
    {
    case 'A': 
      {
	iGene = iRow;
	iPattern = iCol;
	lambda = _AAtomicdomain.getLambda();
	break;
      }
    case 'P': 
      {
	iPattern = iRow;
	iSample = iCol;
	lambda = _PAtomicdomain.getLambda();
	break;
      }
    }

  //bool anyNonzero = false;

  // determine the parameters for finding the mass
  double s  = 0.;
  double su = 0.;
  double mock;

  switch(the_matrix_label){

  case 'A': 
    {
      //double Aeff;
      for (iSample = 0; iSample < _nCol; iSample++) 
	{
	  mock = D[iGene][iSample];
	  for (jPattern = 0; jPattern <_nFactor; jPattern++) 
	    {
	      /*
	      Aeff = currentChainMatrix[iGene][jPattern];
	      if (jPattern == iPattern && origMass < 0.) 
	      //if (_oper_type == 'D' && jPattern == iPattern) 
		{
		  // --- checking ---
		  //if (_oper_type != 'D'){
		  //cout << "Inside getMass (A) re-adding mass, _oper_type = " << _oper_type << endl;
		     //}
		  // -------------
		  Aeff += origMass;
		}
	      mock -= Aeff * otherMatrix[jPattern][iSample];
	      */
              mock -= currentChainMatrix[iGene][jPattern]* otherMatrix[jPattern][iSample];


	    }
	  s += _annealingTemperature * pow(otherMatrix[iPattern][iSample],2) /
	    ( 2* pow(S[iGene][iSample],2) );
	  su += _annealingTemperature * otherMatrix[iPattern][iSample]*mock /
	    ( 2* pow(S[iGene][iSample],2) );	    
	}
      break;
    }

  case 'P': 
    {
      //double Peff;
      for (iGene = 0; iGene < _nRow; iGene++) 
	{
	  mock = D[iGene][iSample];
	  for (jPattern = 0; jPattern < _nFactor; jPattern++) 
	    {
	      /*
	      Peff = currentChainMatrix[jPattern][iSample];
	        if (jPattern == iPattern && origMass < 0.) 
		  //if (_oper_type == 'D' && jPattern == iPattern) 
		{
		  // --- checking ---
		  //if (_oper_type != 'D'){
		  //cout << "Inside getMass (P) re-adding mass, _oper_type = " << _oper_type << endl;
		     //}
		  // -------------
		  Peff += origMass;
		}
	      mock -= otherMatrix[iGene][jPattern] * Peff;
	      */
	      mock -= otherMatrix[iGene][jPattern] * currentChainMatrix[jPattern][iSample];

	    }
	  s += _annealingTemperature * pow(otherMatrix[iGene][iPattern],2) / 
	    ( 2* pow(S[iGene][iSample],2) );
	  su += _annealingTemperature * otherMatrix[iGene][iPattern]*mock / 
	    ( 2* pow(S[iGene][iSample],2) );
	}
      break;
    }

  }

    
  double mean  = (2*su - lambda)/(2*s);
  double sd = 1./sqrt(2*s);
    
  // note: is bounded below by zero so have to use inverse sampling!
  // based upon algorithm in DistScalarRmath.cc (scalarRandomSample)

  double plower = sub_func::pnorm(0., mean, sd, DOUBLE_NEGINF, 0);
  double pupper = 1.;
  double u = plower + randgen('U',0,0) * (pupper - plower);
  // -------------------------------------------------------------
  // this line seems to be misplaced.
  // double newMass = sub_func::qnorm(u, mean, sd, DOUBLE_NEGINF, 0);
  double newMass = 0;
  // ------------------

  // if the likelihood is flat and nonzero, 
  // force to sample strictly from the prior
  if ( plower == 1 || s < 1.e-5 || newMass == DOUBLE_POSINF || pupper == 0) { 
    if (origMass < 0) {    // death case
      newMass = abs(origMass);
    } else {
      newMass = 0.;  // birth case
    }
  } // end of first comparison
  else if (plower >= 0.99) {
    double tmp1 = sub_func::dnorm(0, mean, sd, false); 
    double tmp2 = sub_func::dnorm(10*lambda, mean, sd, false);
    if ( (tmp1 > epsilon) && (fabs(tmp1-tmp2) < epsilon) )   {
      if (origMass < 0) {   // death case
	return 0.;
      }
      return origMass;   // birth case
    }  
  } // end of second comparison
  else {
    newMass = sub_func::qnorm(u, mean, sd, DOUBLE_NEGINF, 0);  // both death and birth
  }  // end of if-block for the remaining case


  // limit the mass range
    switch(the_matrix_label) {
    case 'A':
      {
	if (newMass > _max_gibbsmassA)
	  newMass = _max_gibbsmassA;
	break;
      }
    case 'P':
      {
	if (newMass > _max_gibbsmassP)
	  newMass = _max_gibbsmassP;
	break;
      }
    }

  if (newMass < 0) newMass = 0;   // due to the requirement that newMass > 0 

  return newMass;

}





// -----------------------------------------------------------------------------
void GibbsSampler::detail_check(char outputchi2_Filename[]){

      double chi2 = 2.*cal_logLikelihood();
      cout << "oper_type: " << _oper_type <<
              " ,nA: " << getTotNumAtoms('A') <<
	      " ,nP: " << getTotNumAtoms('P') << 
	//" ,chi2 = " << chi2 << 
              " ,_sysChi2 = " << _sysChi2 << endl;


      ofstream outputchi2_File;
      outputchi2_File.open(outputchi2_Filename,ios::out|ios::app);
      outputchi2_File << chi2 << endl;
      outputchi2_File.close();

}
