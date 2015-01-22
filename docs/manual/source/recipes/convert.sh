#!/bin/sh

for i in `ls *.ipynb`; do
	ipython nbconvert --to html $i;
	ipython nbconvert --to python $i;
done
