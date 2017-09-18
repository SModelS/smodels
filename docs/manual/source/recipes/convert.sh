#!/bin/sh

CMD=`which jupyter-nbconvert`; ## first choice: jupyter-nbconvert
if [ "$CMD" = "" ]; then
	CMD=`which jupyter`; ## second choice: jupyter nbconvert
	if [ "$CMD" = "" ]; then
		CMD="`which ipython` nbconvert"; ## third choice: ipython nbconvert
	else
		CMD="`which jupyter` nbconvert";
	fi
fi
echo "CMD: $CMD"
## CMD=jupyter

for i in `ls *.ipynb`; do
	$CMD --to html $i || { echo "\nERROR: cannot execute $CMD nbconvert. Maybe install jupyter-nbconvert?\n\n"; exit; }
  $CMD --to python $i;
done

mkdir -p ../../build/html/_downloads/
cp *.html ../../build/html/_downloads/
