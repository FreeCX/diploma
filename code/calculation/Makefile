PROG=vortex
ALGP=./include/alglib/src
HLBP=./include/HLBFGS
ALGLIB=$(ALGP)/alglibinternal.cpp $(ALGP)/alglibmisc.cpp \
	$(ALGP)/ap.cpp $(ALGP)/dataanalysis.cpp \
	$(ALGP)/diffequations.cpp $(ALGP)/fasttransforms.cpp \
	$(ALGP)/integration.cpp $(ALGP)/interpolation.cpp \
	$(ALGP)/linalg.cpp $(ALGP)/optimization.cpp \
	$(ALGP)/solvers.cpp $(ALGP)/specialfunctions.cpp \
	$(ALGP)/statistics.cpp
HLBFGS=$(HLBP)/HLBFGS_BLAS.cpp $(HLBP)/HLBFGS.cpp \
	$(HLBP)/ICFS.cpp $(HLBP)/LineSearch.cpp \
	$(HLBP)/Lite_Sparse_Matrix.cpp
MAIN=parser.cpp vortex.cpp

SOURCES=$(ALGLIB) $(HLBFGS) $(MAIN)
CFLAGS=-Wall -g -O2 -pipe -march=native -mfpmath=sse -omp -lpthread -I./include
LDFLAGS=-I./include

OBJS := $(patsubst %.cpp, %.o, $(SOURCES))
DEPS := $(patsubst %.o, %.d, $(OBJS))

$(PROG): $(OBJS) $(DEPS)
	$(CXX) $(OBJS) -o $@ $(LDFLAGS)

%.o: %.cpp %.d
	$(CXX) $(CFLAGS) -c $< -o $@

%.d: %.cpp
	@set -e; $(CXX) -M $< | \
		sed -e 's%\($*\)\.o[ :]*%\1.o $@ : %g' > $@; \
		[ -s $@ ] || rm -f $@
	@echo create $@

clean:
	$(RM) $(PROG) $(OBJS) $(DEPS) *~

-include $(DEPS)
