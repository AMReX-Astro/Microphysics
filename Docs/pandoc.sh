#!/bin/sh

for i in `find . -name "*.tex"`
do
    echo "\\input{textsymbols}" > _temp.tex
    cat $i >> _temp.tex
    pandoc _temp.tex --mathjax --wrap=preserve -o ../docs_new/`basename $i .tex`.rst
    rm -f _temp.tex
done

