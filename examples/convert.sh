#!/bin/sh

for i in `ls recipce*.ipynb`; do
	echo $i;
	ipython nbconvert --template full --to python $i;
	ipython nbconvert --template full --to html $i;
done
