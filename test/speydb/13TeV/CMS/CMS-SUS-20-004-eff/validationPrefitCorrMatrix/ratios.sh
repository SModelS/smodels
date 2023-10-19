#!/bin/sh

DIR=./
DIR=../../../../../smodels-utils/covariances/
DBPATH=$HOME/git/smodels-database/

$DIR/plotRatio.py -d $DBPATH -a1 CMS-SUS-20-004 -a2 CMS-SUS-20-004-eff -v1 TChiHH_2EqMassAx_EqMassBy.py -v2 TChiHH_2EqMassAx_EqMassBy_combined.py -l1 UL -l2 SLv2 -o "ratios_@a_@t_xy.png"
$DIR/plotRatio.py -d $DBPATH -a1 CMS-SUS-20-004 -a2 CMS-SUS-20-004-eff -v1 TChiHH_2EqMassAx_EqMassB1.0.py -v2 TChiHH_2EqMassAx_EqMassB1.0_combined.py -l1 UL -l2 SLv2
$DIR/plotRatio.py -d $DBPATH -a1 CMS-SUS-20-004 -a2 CMS-SUS-20-004-eff -v1 T5HH_2EqMassAx_EqMassBx-50_EqMassC1.0.py -v2 T5HH_2EqMassAx_EqMassBx-50_EqMassC1.0_combined.py -l1 UL -l2 SLv2
