## /sed_sam ## 
## /sed_sam_zoe (replaced with final_sam_20200630 after re-running dada2 ##

contains modified sam files after running dada2 (workflow described in ../scripts/dada2_all.R)

analyze_dadasam function (../scripts/fn_dadasam.R) will read these files and put variants and counts data into a dataframe.


=====================================================

## files (output from ../scripts/dada2_all.R) ##

dada2result_yesol.csv & dada2result_zoe.csv:
Row = individual codon mutation
Columns:
    sample = indicates time point and replicate
    var = label from dada2 for each denoised variant
    nmut = number of codon mutations
           (increments by 1 if multiple rows (var with multiple muts)
    count = number of read counts per variant
    site = AA site with start codon as 1
    AA_change
    cod_change
    in_target = mut matches JT's primers list

dada2result_yesol2.csv & dada2result_zoe2.csv:
Row = individual variant
Columns:
    sample = same as previous
    var = same as previous
    AAsub = AA change information with residue number
            Multiple muts are separated by comma
    count = same as previous
    nmut_tot = total # of codon changes in the variant
    nmut_offtarget = # of off-target codon changes in the variant

Gmut_fitness_Lu.csv:
    Fitness assay results from some G mutants.
    Original file, "prot stab data_BEACON(1).xlsx" received from Lu on 10/01/20
    --> took out score and codon information on 6 mutants and saved separately to compare with my data.



