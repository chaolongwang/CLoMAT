#!/bin/sh

casetype=smooth #spatial risk distribution, should be set to "smooth" or "sharp".
Match=1 #matching scheme, should be set to 1 (1:1 matching), 3 (1:3 matching), or fullmatch.
loopindex=1000 #simulation replicates
DIR=$WK_DIR #working directory

#step 1 simulate phenotype
Rscript $DIR/Null_step1_simulate_phenotype.R $casetype $DIR $loopindex

#step 2 perform matching
for ((i=1; i<=$loopindex; i=$[i+1]))
do {
Rscript $DIR/Null_step2_perform_matching.R $Match ${i} $casetype $DIR
}
done
wait

#step 3 perform test
for ((i=1; i<=$loopindex; i=$[i+1]))
do {
Rscript $DIR/Null_step3_perform_test.R $Match ${i} $casetype $DIR
}
done
wait
