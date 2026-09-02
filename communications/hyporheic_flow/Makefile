.Phony: all

MAIN ?= slides
OUTDIR ?= _latex
BIB ?= References.bib
ENGINE ?= xelatex

main:
	mkdir $(OUTDIR) || true
	$(ENGINE) -output-directory=$(OUTDIR) --shell-escape $(MAIN).tex
	mv $(OUTDIR)/$(MAIN).pdf .

bib:
	cp $(BIB) $(OUTDIR)/$(BIB)
	cd $(OUTDIR) && bibtex $(MAIN)
# bibtex $(OUTDIR)/$(MAIN).aux

clean:
	rm -r $(OUTDIR) || true
	rm -r svg-inkscape || true

all:
	make clean
	mkdir _latex || true
	make main
	make bib
	make main
	make main

svg:
	make main
	pdf2svg $(MAIN).pdf $(MAIN).svg