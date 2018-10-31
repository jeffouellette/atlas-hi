#------------------------------------------------------------------------------
# Makefile to make atlas_tables documentation
#------------------------------------------------------------------------------

BASENAME = atlas_tables
TEXLIVE = 2016

.PHONY: clean cleanpdf

atlas_tables:
	if [ $(TEXLIVE) -ge 2011 ]; then \
		make -f ../Makefile BASENAME=$(BASENAME); \
	fi

texlive:
	make -f ../Makefile texlive BASENAME=$(BASENAME) TEXLIVE=$(TEXLIVE)

clean:
	make -f ../Makefile clean
	
cleanpdf:
	-rm $(BASENAME).pdf

cleanall: clean cleanpdf
	