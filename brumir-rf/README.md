# Random Forest classifier of Brumir Predictions

## Building a random Forest from MirGenDB

### Step 1

Build the input matrix using the script ***compute_features_for_trainingRF.pl***

```
perl compute_features_for_trainingRF.pl ALL.fas > MirGenDB.rf.training.data
``` 

### Step 2 

use the MirGenDB.rf.training.data data to train and evaluate a random forest model with the RF_miRNAs-miRGenDB.Rmd 

##  Using a random forest to classify BRumiR predictions 

Use the RF model to classify Brumir predictions with the script ***classify_brumir_with_RF.pl***

```
perl classify_brumir_with_RF.pl -i SRR1734817-k14-BrumiR-rev-2-d50.candidate_miRNA.fasta -d ALL.fas -m MirGenDB_rf_model_v01.rds -s SRR1734817
```

To classify another BrumiR output is just necesary to rerun the script ***classify_brumir_with_RF.pl*** changing the -i and -s options (BrumiR result and sample name)

