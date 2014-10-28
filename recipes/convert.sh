#!/bin/sh

for i in `ls recipe*.ipynb`; do
	echo $i;
	ipython nbconvert --template full --to python $i;
	ipython nbconvert --template full --to html $i;
done
