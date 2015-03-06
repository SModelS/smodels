#!/bin/sh

for i in `ls *.ipynb`; do
	echo 
	echo
	echo $i;
	cmd="ipython nbconvert --to python $i"
	echo $cmd;
  `$cmd`;
	cmd="ipython nbconvert --to html $i"
	echo $cmd;
  `$cmd`;
done
