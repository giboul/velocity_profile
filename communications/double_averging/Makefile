.Phony: all

MAIN ?= averaging_rules.tex
OUTDIR ?= _latex
BIB ?= references.bib
ENGINE ?= xelatex
README ?= README.md

MAIN_BASE = $(basename $(MAIN))
MAIN_STEM = $(notdir $(MAIN_BASE))
MAIN_LATEX = $(MAIN_BASE).tex

LATEX_FLAGS = 
LATEX_FLAGS += --halt-on-error
LATEX_FLAGS += -output-directory=$(OUTDIR)
LATEX_FLAGS += --shell-escape
LATEX_FLAGS += -jobname=$(MAIN_STEM)

main:
	mkdir -p $(OUTDIR)
	$(ENGINE) $(LATEX_FLAGS) $(MAIN_LATEX) $(FLAGS)

svg:
	make main
	pdf2svg $(MAIN_BASE).pdf $(MAIN_BASE).svg

bib:
	cp $(BIB) $(OUTDIR)/$(BIB)
	cd $(OUTDIR) && bibtex $(MAIN_STEM)

all:
	make main FLAGS=-no-pdf
	make bib
	make main FLAGS=-no-pdf
	make main
	mv $(OUTDIR)/$(MAIN_STEM).pdf .

clean:
	rm -r $(OUTDIR) || true
	rm -r svg-inkscape || true
	rm *.pdf