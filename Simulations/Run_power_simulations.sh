#!/bin/sh

CausalPercent=0.2 # casual percent, should be set to 0.2 or 0.5.
PositivePercent=0.8 # protective percant within causal variants, should be set to 0.8 or 1.0.
Constant=0.877 # set to 0.877 when CausalPercent=0.2; set to 0.555 when CausalPercent=0.5.
casetype=smooth # spatial risk distribution, should be set to "smooth" or "sharp".
Match=1 #matching scheme, should be set to 1 (1:1 matching), 3 (1:3 matching), or fullmatch.  
DIR=$WK_DIR # working directory.
loopindex=1000 #simulation replicates.

#step 1 simulate phenotype
Rscript $DIR/Power_step1_simulate_phenotype.R $CausalPercent $PositivePercent $Constant $casetype $DIR $loopindex

#step 2 perform matching
for ((i=1; i<=$loopindex;  i=$[i+1]))
do {
Rscript $DIR/Power_step2_perform_matching.R $Match ${i} $casetype $DIR 
} 
done
wait

#step 3 perform test
Rscript $DIR/Power_step3_perform_test.R $Match $casetype $DIR $loopindex



















