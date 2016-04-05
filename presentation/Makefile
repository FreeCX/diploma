COMPILER=latexmk
SOURCE=presentation.tex

default:
	$(COMPILER) -pdf $(SOURCE)
clean:
	$(COMPILER) -C
	$(RM) *.{aux,fls,aux,fdb_latexmk,log,ptc,nav,snm}
