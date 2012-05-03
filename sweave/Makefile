all: tidy
	R CMD BATCH Build.R
	mv ../MCSMAC* .
	R CMD pdflatex MCSMAC.tex
	bibtex MCSMAC
	R CMD pdflatex MCSMAC.tex


tangle: tidy
	R CMD BATCH Tangle.R
	mv ../MCSMAC.R MCSMAC.R

tidy:
	rm -rf *Rout *~ *pdf *aux *log *tex MCSMAC.R .RData *.bbl *.bak
	rm -rf ../MCSMAC*

clean: tidy
	rm -rf cacheSweave*
