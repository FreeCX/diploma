COMPILER=latexmk
FLAGS=-jobname=specification
SOURCE=00-specification.tex

default:
	$(COMPILER) -pdf $(FLAGS) $(SOURCE)
clean:
	$(COMPILER) -C
	$(RM) *.{aux,fls,aux,fdb_latexmk,log,ptc}
