EPStoPDF = epstopdf

ALL: userguide.pdf 

eps_source = $(wildcard *.eps) $(wildcard */*.eps)

pdf_source = $(eps_source:.eps=.pdf)

section_sources = userguide.tex

userguide.dvi: userguide.tex $(eps_source) $(section_sources) $(eps_source)
	latex userguide.tex < /dev/null
	bibtex userguide
	latex userguide.tex < /dev/null
	latex userguide.tex < /dev/null

userguide.pdf: userguide.tex $(pdf_source) $(section_sources) $(pdf_source)
	pdflatex userguide.tex < /dev/null
	bibtex userguide
	pdflatex userguide.tex < /dev/null
	pdflatex userguide.tex < /dev/null

pdf:	userguide.pdf 

%.ps: %.dvi
	dvips -t letter -o $@ $<

%.pdf: %.eps
	$(EPStoPDF) $<

clean:
	$(RM) $(pdf_source) userguide.dvi
	$(RM) userguide.blg userguide.log userguide.toc
	$(RM) userguide.aux userguide.ps userguide.bbl userguide.pdf
	$(RM) userguide.lot userguide.lof userguide.out
	$(RM) *~

.PHONY: clean
