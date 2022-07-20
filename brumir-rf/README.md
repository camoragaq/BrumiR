# Random Forest classifier for Brumir Predictions

## Building a random Forest from MirGeneDB



### Step 1


Gather the MirGeneDB data with:

```
wget https://mirgenedb.org/fasta/ALL?mat=1
```

Generate the ***ALL.fas*** file

### Step 2

Build the input matrix using the script ***compute\_features\_for\_trainingRF.pl***

```
perl compute_features_for_trainingRF.pl ALL.fas > MirGenDB.rf.training.data
``` 

### Step 3 

use the ***MirGenDB.rf.training.data*** data to train and evaluate a random forest model with the R code ***RF\_miRNAs-miRGenDB.Rmd*** 

The R code generate the random forest model called ***MirGenDB\_rf\_model\_v01.rds***

###  Step 4 

Use the RF model ***MirGenDB\_rf\_model\_v01.rds***  to classify Brumir predictions with the script ***classify\_brumir\_with\_RF.pl***

```
perl classify_brumir_with_RF.pl \\
-i SRR1734817-k14-BrumiR-rev-2-d50.candidate_miRNA.fasta \\
-d ALL.fas -m MirGenDB_rf_model_v01.rds -s SRR1734817
```

To classify another BrumiR output is just necessary to rerun the script ***classify\_brumir\_with\_RF.pl*** changing the -i and -s options (BrumiR result and sample name)


## Building a random forest from MiRbase (plants)

### Step 1 

Gather the miRbase data and select plant species with:

```
wget https://www.mirbase.org/ftp/CURRENT/mature.fa.gz
wget https://www.mirbase.org/ftp/CURRENT/organisms.txt.gz
zless organisms.txt.gz | grep "Viridiplantae" > plants.mirbase.txt
gzip -d mature.fa.gz
#get miRbase plant entries from mirbase mature.fa
python3 get_plants_entries_mirbase.py 
```

The above commands generate the file ***miRbase.mature.plants.fa*** with a total of 10414 miRNA mature plant sequences.

### Step 2 


Build the input matrix using the script ***compute\_features\_for\_trainingRF.pl***

```
perl compute_features_for_trainingRF.pl miRbase.mature.plants.fa > MirBase.plants.rf.training.data
``` 

### Step 3

use the ***MirBase.plants.rf.training.data*** data to train and evaluate a random forest model with the R code ***RF\_miRNAs-miRbase.Rmd***

The R code generate the random forest model called ***miRbase\_plants\_rf\_model\_v01.rds***

### Step 4 

Use the RF model  ***miRbase\_plants\_rf\_model\_v01.rds*** to classify plant Brumir predictions with the script ***classify\_brumir\_with\_RF.pl***

```
perl classify_brumir_with_RF.pl \\
-i  \\
-d miRbase.mature.plants.fa -m miRbase_plants_rf_model_v01.rds -s 
```

To classify another palnt BrumiR output is just necessary to rerun the script ***classify\_brumir\_with\_RF.pl*** changing the -i and -s options (BrumiR result and sample name)


