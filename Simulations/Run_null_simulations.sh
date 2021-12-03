#!/bin/sh

loopindex=1000 #simulation replicates
DIR=$WK_DIR #working directory

for casetype in smooth sharp
do 

#step 1 simulate phenotype
Rscript $DIR/Null_step1_simulate_phenotype.R $casetype $DIR $loopindex

#step 2 perform matching 
#matching scheme should be set to 1 (1:1 matching), 3 (1:3 matching), or fullmatch.
for Match in 1 3 fullmatch
do
for ((i=1; i<=$loopindex; i=$[i+1]))
do 
Rscript $DIR/Null_step2_perform_matching.R $Match ${i} $casetype $DIR
done
done
wait

#step 3 perform test
for Match in 1 3 fullmatch
do
for ((i=1; i<=$loopindex; i=$[i+1]))
do
Rscript $DIR/Null_step3_perform_test.R $Match ${i} $casetype $DIR
done
done
wait

done