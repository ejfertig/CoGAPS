#ifndef GAPSNORM_H_
#define GAPSNORM_H_

#include <iostream>
#include <vector>
#include <boost/tuple/tuple.hpp>

//using namespace gaps;
using namespace std;
//using std::vector;

namespace gaps
{

  class GAPSNorm
  {

  private:
    static void computeMock(double ** M, double const * const * A, double const * const * P,
			    unsigned int nRow, unsigned int nCol, unsigned int nFactor); 

  public:
    GAPSNorm();
    ~GAPSNorm();

    void local_display_matrix(double const * const * Mat, unsigned int n_row,
				     unsigned int n_col);
 
    static double calChi2(double const * const * D, double const * const * S,
			  double const * const * A, double const * const * P,
			  unsigned int nRow, unsigned int nCol,
			  unsigned int nFactor);

    static double calcDeltaLL1E(char matrix_label,
			      double const * const * D, double const * const * S, 
			      double const * const * A, double const * const * P, 
			      const vector<boost::tuple<unsigned int, unsigned int, double> > ElemChange, 
			      unsigned int nRow, unsigned int nCol, unsigned int nFactor); 

    static double calcDeltaLL2E(char matrix_label,
			      double const * const * D, double const * const * S, 
			      double const * const * A, double const * const * P, 
			      const vector<boost::tuple<unsigned int, unsigned int, double> > ElemChange, 
			      unsigned int nRow, unsigned int nCol, unsigned int nFactor);

    static double calcDeltaLLGen(char matrix_label,
			      double const * const * D, double const * const * S, 
			      double const * const * A, double const * const * P, 
			      const vector<boost::tuple<unsigned int, unsigned int, double> > ElemChange, 
			       unsigned int nRow, unsigned int nCol, unsigned int nFactor);


      /*
	static double logLikelihoodRatio(const double ** D, const double ** S,
	const double ** A, const double ** P,
	const vector<double> newMatrix, bool delA,
	unsigned int nRow, unsigned int nCol,
	unsigned int nFactor);

	static double logLikelihoodRatioAtomic(const double * D, const double * S,
	const double * A, const double * P,
	unsigned int iRowChanged[],
	unsigned int iColChanged[],
	double massesChanged[],
	unsigned int nChange,
	bool delA,
	unsigned int nRow, unsigned int nCol,
	unsigned int nFactor);

      */
      };
}
#endif 
