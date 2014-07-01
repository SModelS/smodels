#!/bin/sh

for i in `ls example*.ipynb`; do
	echo $i;
	ipython nbconvert --template full --to python $i;
	ipython nbconvert --template full --to html $i;
done
