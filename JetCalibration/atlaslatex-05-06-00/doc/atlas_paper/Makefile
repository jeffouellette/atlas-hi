#------------------------------------------------------------------------------
# Makefile to make atlas_paper documentation
#------------------------------------------------------------------------------

BASENAME = atlas_paper
TEXLIVE = 2016

.PHONY: clean cleanpdf

atlas_paper:
	make -f ../Makefile BASENAME=$(BASENAME)

texlive:
	make -f ../Makefile texlive BASENAME=$(BASENAME) TEXLIVE=$(TEXLIVE)

clean:
	make -f ../Makefile clean
	
cleanpdf:
	-rm $(BASENAME).pdf

cleanall: clean cleanpdf
