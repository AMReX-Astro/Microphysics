EPStoPDF = epstopdf

ALL: userguide.pdf

eps_source = $(wildcard *.eps) $(wildcard */*.eps)

pdf_source = $(eps_source:.eps=.pdf)

png_source := $(wildcard */*.png)

section_sources = userguide.tex \
                  basics/basics.tex \
                  data_structures/data_structures.tex \
                  EOS/eos.tex \
                  runtime_parameters/rp_intro.tex \
                  runtime_parameters/runtime_parameters.tex \
                  networks/networks.tex \
                  unit_tests/unit_tests.tex 

userguide.pdf: userguide.tex $(pdf_source) $(section_sources) $(pdf_source) $(png_source) refs.bib
	pdflatex userguide.tex < /dev/null
	bibtex userguide
	pdflatex userguide.tex < /dev/null
	pdflatex userguide.tex < /dev/null

pdf:	userguide.pdf

%.pdf: %.eps
	$(EPStoPDF) $<

clean:
	$(RM) $(pdf_source) userguide.dvi
	$(RM) userguide.blg userguide.log userguide.toc
	$(RM) userguide.aux userguide.ps userguide.bbl userguide.pdf
	$(RM) userguide.lot userguide.lof userguide.out
	$(RM) *~

.PHONY: clean
