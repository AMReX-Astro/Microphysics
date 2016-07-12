EPStoPDF = epstopdf

ALL: MicrophysicsUsersGuide.pdf

eps_source = $(wildcard *.eps) $(wildcard */*.eps)

pdf_source = $(eps_source:.eps=.pdf)

png_source := $(wildcard */*.png)

section_sources = MicrophysicsUsersGuide.tex \
                  basics/basics.tex \
                  data_structures/data_structures.tex \
                  EOS/eos.tex \
                  runtime_parameters/rp_intro.tex \
                  runtime_parameters/runtime_parameters.tex \
                  networks/networks.tex \
                  unit_tests/unit_tests.tex 

MicrophysicsUsersGuide.pdf: MicrophysicsUsersGuide.tex $(pdf_source) $(section_sources) $(pdf_source) $(png_source) refs.bib
	pdflatex MicrophysicsUsersGuide.tex < /dev/null
	bibtex MicrophysicsUsersGuide
	pdflatex MicrophysicsUsersGuide.tex < /dev/null
	pdflatex MicrophysicsUsersGuide.tex < /dev/null

pdf:	MicrophysicsUsersGuide.pdf

%.pdf: %.eps
	$(EPStoPDF) $<

clean:
	$(RM) $(pdf_source) MicrophysicsUsersGuide.dvi
	$(RM) MicrophysicsUsersGuide.blg MicrophysicsUsersGuide.log MicrophysicsUsersGuide.toc
	$(RM) MicrophysicsUsersGuide.aux MicrophysicsUsersGuide.ps MicrophysicsUsersGuide.bbl MicrophysicsUsersGuide.pdf
	$(RM) MicrophysicsUsersGuide.lot MicrophysicsUsersGuide.lof MicrophysicsUsersGuide.out
	$(RM) *~

.PHONY: clean
