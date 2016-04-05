COMPILER=latexmk
FLAGS=-jobname=00-synopsis
SOURCE=00-synopsis.tex

default:
	$(COMPILER) -pdf $(FLAGS) $(SOURCE)
clean:
	$(COMPILER) -C
	$(RM) *.{aux,fls,aux,fdb_latexmk,log,ptc}
