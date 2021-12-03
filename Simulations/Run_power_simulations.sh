#!/bin/sh

DIR=$WK_DIR # working directory.
loopindex=1000 #simulation replicates.

for casetype in smooth sharp
do 

for CausalPercent in  `echo "0.2"; echo "0.5"`; do

for PositivePercent in  `echo "1.0"; echo "0.8"`; do

if [ $CausalPercent == "0.2" ]
then
 Constant=0.877 
else
 Constant=0.555
fi

#step 1 simulate phenotype
Rscript $DIR/Power_step1_simulate_phenotype.R $CausalPercent $PositivePercent $Constant $casetype $DIR $loopindex

#step 2 perform matching
#matching scheme should be set to 1 (1:1 matching), 3 (1:3 matching), or fullmatch.
for Match in 1 3 fullmatch
do
for ((i=1; i<=$loopindex;  i=$[i+1]))
do
Rscript $DIR/Power_step2_perform_matching.R $CausalPercent $PositivePercent $Constant $Match ${i} $casetype $DIR 
done
done
wait

#step 3 perform test
for Match in 1 3 fullmatch
do
Rscript $DIR/Power_step3_perform_test.R $CausalPercent $PositivePercent $Constant $Match $casetype $DIR $loopindex
done

done

done

done


















