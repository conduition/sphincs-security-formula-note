# Usage:
#   make            # build main.pdf
#   make clean      # remove build artifacts

TEX  = pdflatex
BIBC = bibtex
MAIN = spx_sec_formula
BIB  = references

BUILDTX  = $(MAIN).tex $(BIB).bib

.PHONY: all clean

all: $(MAIN).pdf

$(MAIN).pdf: $(BUILDTX)
	$(TEX) -halt-on-error $(MAIN).tex
	-$(BIBC) $(MAIN)
	$(TEX) -halt-on-error $(MAIN).tex
	$(TEX) -halt-on-error $(MAIN).tex
	$(TEX) -halt-on-error $(MAIN).tex

clean:
	rm -f $(MAIN).aux $(MAIN).log $(MAIN).out $(MAIN).toc $(MAIN).pdf \
	      $(MAIN).bbl $(MAIN).blg $(MAIN).fls $(MAIN).fdb_latexmk
