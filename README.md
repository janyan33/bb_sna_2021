# Data analysis journal for the bedbug social network project
Last updated: 2021-09-20

## Data
I've combined all (day and night) live observation data into one spreadsheet named bbsna_raw_combined.csv  
This spreadsheet contains all information on matings and mountings. 


## To do
Check if "attempted mount" should be treated as mounts. If so, add to patch table. Need to ask students what they meant. Same with "tried to mount". 

Re-run mating and mounting stats now that we have the light-time data added. 

Ask Reuven about the strength stats (whether the GLMM makes sense, whether the dianostic plots matter, how to permute the contrasts etc.)


## Strength analysis
Using the bbsna_aggregations file, I've created aggregation based-social networks for each replicate and exported every individual's strength value into the master attribute file. 
I then used a glmm to look at strength as a function of sex and treatment (and the interaction) with block and size as random effects. The diagnostic plots (especially the QQ plot) looks wonky but I don't think this matters if I significance test via permutations? 

As for the permutation test for each of the two treatments, I used emmeans to extract the z-score for the specific contrast and used that as the test statistic. 
I'm still unsure about whether this makes sense. 



