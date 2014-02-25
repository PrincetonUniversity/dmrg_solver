CXX=g++ 
CXXFLAGS=-O3 -fopenmp -D_HAS_CBLAS -D_HAS_INTEL_MKL -std=c++11

BLASDIR=#/opt/intel/composer_xe_2013.0.079/mkl
BLASINC=#-I$(BLASDIR)/include
BLASLIB=-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lmkl_core

BOOSTDIR=/home/boxiao/usr
BOOSTINC=-I$(BOOSTDIR)/include
BOOSTLIB=-L$(BOOSTDIR)/lib -lboost_serialization

BTASDIR=/home/boxiao/mps/btas-master
BTASINC=-I$(BTASDIR)/include
BTASLIB=-L$(BTASDIR)/lib -lbtas

NEWMATLIB=-L./newmat10 -lnewmat

INCLUDEFLAGS=-I. $(BLASINC) $(BOOSTINC) $(BTASINC)
LIBRARYFLAGS=    $(BLASLIB) $(BOOSTLIB) $(NEWMATLIB)

SRC_SAMPLE = main.C dmrg.C driver.C btas_template_specialize.C mpogen.C mpogen2.C
OBJ_SAMPLE = $(SRC_SAMPLE:.C=.o)

.C.o	:
	$(CXX) -c $< -o $@ $(CXXFLAGS) $(INCLUDEFLAGS)

all	: dmrg.x

dmrg.x	: $(OBJ_SAMPLE)
	$(CXX) $(OBJ_SAMPLE) -o dmrg.x $(CXXFLAGS) $(LIBRARYFLAGS) $(BTASLIB)

clean	:
	rm *.o; rm *.x; rm *.a;
