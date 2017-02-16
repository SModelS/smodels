#!/bin/sh

# CMD=ipython
CMD=jupyter

for i in `ls *.ipynb`; do
	$CMD nbconvert --to html $i;
  $CMD nbconvert --to python $i;
done

mkdir -p ../../build/html/_downloads/
cp *.html ../../build/html/_downloads/
