#!/bin/sh

for i in `ls */*/convert.py`; 
	do echo $i;
	DIR=`dirname $i`;
	BASE=`basename $i`;
	cd $DIR && $BASE;
	cd ../..;
done
