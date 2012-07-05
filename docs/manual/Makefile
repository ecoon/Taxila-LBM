#========================================================================
# UNKNOWN-file
#    author:              Ethan T. Coon
#    filename:            Makefile
#    version:             
#    created:             12 July 2011
#      on:                11:05:46 MDT
#    last modified:   12 July 2011
#      at:            11:05:59 MDT
#    URL:                 http://www.ldeo.columbia.edu/~ecoon/
#    email:               ecoon _at_ lanl.gov
#
#========================================================================
all: manual

manual:
	@echo "Making manual"
	pdflatex manual
	bibtex manual
	pdflatex manual
	pdflatex manual

fresh: clean manual clean

clean:
	rm -f *.dvi *.aux *.log *.bbl *.blg *.out *.snm *.nav *.toc *~;


