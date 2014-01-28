CC=g++
PROG=main
LIB=alglib
SOURCES=$(LIB)/ap.cpp $(LIB)/alglibinternal.cpp $(LIB)/specialfunctions.cpp main.cpp
CFLAGS=-Wall -O3
LDFLAGS=
ifeq ($(shell uname -s),Linux)
	CFLAGS+=-lGL -lGLU -lglut -lm
	EXT=
else
	CFLAGS+=-lopengl32 -lglu32 -lfreeglut -lm
	EXT=.exe
endif

OBJS := $(patsubst %.cpp, %.o, $(SOURCES))
DEPS := $(patsubst %.o, %.d, $(OBJS))

$(PROG): $(OBJS) $(DEPS)
	$(CC) $(OBJS) -o $@ $(CFLAGS)

%.o: %.cpp %.d
	$(CC) $(LDFLAGS) -c $< -o $@

%.d: %.cpp
	@set -e; $(CC) -M $< | \
		sed -e 's%\($*\)\.o[ :]*%\1.o $@ : %g' > $@; \
		[ -s $@ ] || rm -f $@
	@echo create $@

clean:
	$(RM) $(PROG) $(OBJS) $(DEPS)
	$(RM) *~

-include $(DEPS)
