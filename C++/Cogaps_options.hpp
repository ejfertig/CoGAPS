#ifndef _COGAPS_OPTIONS_HPP_
#define _COGAPS_OPTIONS_HPP_


#include <iostream>
#include <fstream>
#include <iomanip>
#include <boost/program_options.hpp>
#include <limits>

namespace opt = boost::program_options;

using namespace std;

 

//exception for obligatory parameters if they are missing
struct nFactor_exception:exception { 
  virtual const char * what() const throw() {
    return "nFactor is undefined\n";
  } 
};
struct datafile_exception:exception { 
  virtual const char * what() const throw() {
    return "datafile is undefined\n";
  } 
};
struct variancefile_exception:exception { 
  virtual const char * what() const throw() {
    return "variancefile is undefined\n";
  } 
};
struct simulation_id_exception:exception { 
  virtual const char * what() const throw() {
    return "#Oh hell, no simulation I.D.!!!???### \n";
  } 
};


class Cogaps_options_class{

private:
  Cogaps_options_class(); //forbids default constructor
  Cogaps_options_class(const Cogaps_options_class & ); //forbids copy constructor
  Cogaps_options_class & operator=(const Cogaps_options_class & ); //forbids copy operator
  opt::variables_map vm;
  opt::options_description description;
  // other reference variables:
  unsigned long unsig_long_max; 

public:
  bool initialized,to_help;
  //obligatory values
  unsigned int nFactor;
  string datafile;
  string variancefile;
  string simulation_id;
  //all the values until here are obligatory
  unsigned long nEquil, nSample;
  unsigned long long max_atoms;
  double alphaA,alphaP,epsilon;
  double nMaxA,nMaxP;
  unsigned long nIterA, nIterP;
 	
  //	Cogaps_options_class(int ac, char* av[]) throw(exception,genes_exception):
  Cogaps_options_class(int ac, char* av[]):
    description("A template example for Cogaps/c++, config file or command line parser"),
    initialized(false),
    to_help(false) {
                
    unsig_long_max = std::numeric_limits<unsigned long>::max(); 
    description.add_options()
      ("help", "produce help message")
      // ---- OBLIGATORY arguments! Must be in!!!! -----------------------------
      ("common.nFactor,F", boost::program_options::value<unsigned int>(&nFactor), 
       "#of patterns")
      ("common.datafile,d", boost::program_options::value<string>(&datafile), 
       "name of the data file")
      ("common.variancefile,v", boost::program_options::value<string>(&variancefile), 
       "name of the sigma file")
      ("common.simulation_id", boost::program_options::value<string>(&simulation_id), 
       "unique id of the program run")
      // ------ NON-OBLIGATORY arguments! If not in, default will be used! -----
      ("common.nEquil,E", boost::program_options::value<unsigned long>(&nEquil)->default_value(500000000ul), 
       "number of iterations to equilibrium")
      ("common.nSample,S", boost::program_options::value<unsigned long>(&nSample)->default_value(unsig_long_max), 
       "maximal number of iterations")
      ("A.alphaA", boost::program_options::value<double>(&alphaA)->default_value(0.01), 
       "prior sparsity of A matrix")
      ("A.nMaxA", boost::program_options::value<double>(&nMaxA)->default_value(unsig_long_max),
       "Max Number of atoms in atomic space A")
      ("A.nIterA", boost::program_options::value<unsigned long>(&nIterA)->default_value(100000),
       "Number of iterations for updating A")
      ("P.alphaP", boost::program_options::value<double>(&alphaP)->default_value(0.01), 
       "prior sparsity of P matrix")
      ("P.nMaxP", boost::program_options::value<double>(&nMaxP)->default_value(unsig_long_max),
       "Max Number of atoms in atomic space P")
      ("P.nIterP", boost::program_options::value<unsigned long>(&nIterP)->default_value(100000),
       "Number of iterations for updating P")
      //("constants.epsilon,e", boost::program_options::value<double>(&epsilon)->default_value(1e-10), 
      //"the minimal difference for nonequal doubles")
      ("config-file", boost::program_options::value<string>(), 
       "configuration file name")
      // -------------------------------------------------------------------------------------------------------
      ;

    opt::positional_options_description pd; 
    pd.add("config-file", -1);
    //we need this to provide the config file name as positional parameter
		
    //here, the variables will be stored in their test repesentation after the parsing
    opt::store(opt::command_line_parser(ac, av).options(description).positional(pd).run(), vm);
    //here, the command-line parameters are read and stored
    //ac and av are command-line parameters passed to main
    opt::notify(vm);    
    //magic. no idea. do not forget to do it

    if (vm.count("help") || vm.count("?")) {
      to_help=true;
      return;
    }
    //'--help' was in command-line


    //we test whether we have --config-file switch or positional nonamed swich 
    //(first three lines of the try block was to see it)
    //if we have it, we parse the config file
    if (vm.count("config-file")) {
      ifstream conf(vm["config-file"].as<string>().c_str());
      opt::store(opt::parse_config_file(conf, description, true), vm);
    }

    opt::notify(vm);
    //magic. no idea. do not forget to do it


    // WARNING statements when the OBLIGATORY arguments are not given: ---------
    //		if (!vm.count("common.genes"))
    //{
    //	throw (*new genes_exception());
    //};
		
    //if (!vm.count("common.samples"))
    //{
    //	throw (*new samples_exception());
    //};

    if (!vm.count("common.nFactor"))
      {
	throw (*new nFactor_exception());
      };

    if (!vm.count("common.datafile"))
      {
	throw (*new datafile_exception());
      };

    if (!vm.count("common.variancefile"))
      {
	throw (*new variancefile_exception());
      };

    if (!vm.count("common.simulation_id"))
      {
	throw (*new simulation_id_exception());
      };
    // ---------------------------------------------------------------------------------------------------------

    initialized=true;
  };
	
  void const help(ostream & of)
  {
    of<<description;
  }

};



ostream & operator << (ostream & o, const Cogaps_options_class & options)
{
  if (!options.initialized) return o;
  o<<
    "[common]"<<endl<<


    "nFactor="<<options.nFactor<<endl<<
    "datafile="<<options.datafile<<endl<<
    "variancefile="<<options.variancefile<<endl<<
    "simulation_id="<<options.simulation_id<<endl<<
    "nEquil="<<options.nEquil<<endl<<
    "nSample="<<options.nSample<<endl<<


    "[A]"<<endl<<
    "alphaA="<<options.alphaA<<endl<<
    "nMaxA="<< options.nMaxA << endl <<
    "nIterA="<< options.nIterA << endl << endl <<

    "[P]"<<endl<<
    "alphaP="<<options.alphaP<<endl<<
    "nMaxP="<< options.nMaxP << endl <<
    // "nIterP="<< options.nIterP << endl << endl <<
    "nIterP="<< options.nIterP << endl << endl;

    // "[constants]"<<endl<<
    //"epsilon="<<options.epsilon<<endl;
  return o;
}
#endif // _COGAPS_OPTIONS_HPP_

