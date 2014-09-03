include ./boostdirs

.PHONY: all

#all: $(name)
	
cogapsmap: cogapsmap.cpp GibbsSamplerMap.o GibbsSampler.o GAPSNorm.o Matrix.o AtomicSupport.o sub_func.o randgen.o
	g++ -I $(boost_include) -o cogapsmap cogapsmap.cpp GibbsSamplerMap.o GibbsSampler.o GAPSNorm.o Matrix.o AtomicSupport.o sub_func.o randgen.o -L$(boost_lib) -lboost_program_options

GibbsSamplerMap.o: GibbsSamplerMap.cpp GibbsSamplerMap.h GibbsSampler.o
	g++ -c -I $(boost_include) GibbsSamplerMap.cpp GibbsSamplerMap.h GibbsSampler.o -L$(boost_lib)

GibbsSampler.o: GibbsSampler.cpp GibbsSampler.h GAPSNorm.o randgen.o sub_func.o Matrix.o AtomicSupport.o
	g++ -c -I $(boost_include) GibbsSampler.cpp GibbsSampler.h GAPSNorm.o randgen.o sub_func.o Matrix.o AtomicSupport.o -L$(boost_lib)

GAPSNorm.o: GAPSNorm.cpp GAPSNorm.h
	g++ -c -I $(boost_include) GAPSNorm.cpp GAPSNorm.h -L$(boost_lib)

Matrix.o: Matrix.cpp Matrix.h
	g++ -c -I $(boost_include) Matrix.cpp Matrix.h -L$(boost_lib)

AtomicSupport.o: AtomicSupport.cc AtomicSupport.h sub_func.o randgen.o
	g++ -c -I $(boost_include) AtomicSupport.cc AtomicSupport.h sub_func.o randgen.o -L$(boost_lib)

sub_func.o: sub_func.cpp sub_func.h randgen.o
	g++ -c -I $(boost_include) sub_func.cpp sub_func.h randgen.o -L$(boost_lib)

randgen.o: randgen.cpp randgen.h
	g++ -c -I $(boost_include) randgen.cpp randgen.h -L$(boost_lib)