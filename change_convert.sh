#!/bin/zsh

for i in `ls **/convert.py`; do
    echo $i;
		cp $i ${i}.old
		cat $i | sed -e 's/prettyname/prettyName/' | sed -e 's/superseded_by/supersededBy/' | sed -e 's/\.condition/.conditionDescription/' | sed -e 's/fuzzycondition/condition/' | sed -e 's/implemented_by/implementedBy/' | tee $i
done
