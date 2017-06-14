#!/bin/sh

CMD=`which jupyter`
[ "$CMD"=="" ] && CMD=`which ipython`;
echo $CMD
## CMD=jupyter

for i in `ls *.ipynb`; do
	$CMD nbconvert --to html $i || echo "\nERROR: cannot execute $CMD nbconvert. Maybe install $CMD-nbconvert?\n\n";
  $CMD nbconvert --to python $i;
done

mkdir -p ../../build/html/_downloads/
cp *.html ../../build/html/_downloads/
