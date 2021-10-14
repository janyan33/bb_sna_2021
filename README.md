# Data analysis journal for the bedbug social network project
Last updated: 2021-10-14

## Data
I've combined all (day and night) live observation data into one spreadsheet named bbsna_raw_combined.csv  
This spreadsheet contains all information on matings and mountings  
The file bbsna_attributes_full.csv contains the attributes for each individual, some of which (strength, matings, mountings etc) were added after they were computed in R. 

## To do
Figure out how to properly do the permutations for each contrast in the strength analysis (one permutation test for M vs. F per treatment).  
Also write a script that analyzes assortativity (for both association based on touch and shared shelters). 

## Strength analysis
Using the bbsna_aggregations file, I've created aggregation based-social networks for each replicate and exported every individual's strength value into the master attribute file. 
I then used a glmm to look at strength as a function of sex and treatment (and the interaction) with block and size as random effects. The diagnostic plots (especially the QQ plot) looks wonky but I don't think this matters if I significance test via permutations? 

As for the permutation test for each of the two treatments, I used emmeans to extract the z-score for the specific contrast and used that as the test statistic. 
I'm still unsure about whether this makes sense. 

## Mating and mounting analysis
Using mate_mount_data.csv which comes from the raw_data_combined.csv file which is cleaned up in the data_cleaning.R script,  
I calculated how many matings and mountings each individual received and inputted those values into the bbsna_attributes_full.csv  
file. I then ran linear mixed models to look at mating and mounting rates as a function of treatment with block as a random factor.  

## Assortment analysis
To be completed. 