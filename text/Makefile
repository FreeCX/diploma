COMPILER=latexmk
FLAGS=-jobname=diploma
SOURCE=00-diploma.tex

default:
	$(COMPILER) -pdf $(FLAGS) $(SOURCE)
clean:
	$(COMPILER) -C
	$(RM) *.{aux,fls,aux,fdb_latexmk,log,ptc}
