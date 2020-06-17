# GSE118959
Analysis of GSE118959 circRNAs in LNCaP cell lines
***

# Analysis

### Differentially expressed circRNAs
Run `DE_circRNA.Rmd` script to output results for DE circRNAs in Clone1 vs control & Clone9 vs control. 

To get the overlapping circRNAs present in upregulated/downregulated sets, use the following bash command:

```bash
awk -F'\t' 'NR==FNR{c[$2$2]++;next};c[$2$2] > 0' clone1_upreg.txt clone9_upreg.txt > up_regulated_overlap.txt

awk -F'\t' 'NR==FNR{c[$2$2]++;next};c[$2$2] > 0' clone1_downreg.txt clone9_downreg.txt > down_regulated_overlap.txt
```

The above code returns matching lines using column 2 as the key. As this analysis will be an unweighted graph, there is no need to worry about which Fold Change value to take, I am simply interested in the circRNA ID and its genomic coordiantes. 

### Predicting miRNA binding sites
Firstly, I will use targetscan (7.0) to perform this analysis on a test "UTR" file. Note that targetscan requires RNAfold and 2 perl libraries in order to work `Statistics::Lite` and `Bio::TreeIO`. `Bio::TreeIO` was a nasty install, I had to install `cpanminus` (`cpanm`) which will install missing dependencies for the library you are trying to install. It failed, so the `--force` option had to be used along with `sudo` :

```bash
curl -L https://cpanmin.us | perl - --sudo App::cpanminus
sudo cpanm Bio::TreeIO -f -n 
sudo cpanm Statistics::Lite -f -n ##-n no test 
```

```
For fresh install of perl (useful for lugh)  
delete /home/perl5 dir
follow instructions here: https://gist.github.com/ckandoth/1f01d8f3692bb8de7f2929f259a4035f
after install, run curl -L https://cpanmin.us | perl - App::cpanminus
installs automatically to /home/perl5
then install final 2 libraries via cpanm. 
have not tested install without cpanm part. 
```

I am still testing the software. Once it is running, I will have to develop a pipeline to extract DNA sequences, convert to RNA and then format them for targetscan input. (use samtools faidx here, also `sed 's/T/U/g'` for DNA to RNA conversion). 

### TargetScan Pipeline
For your benefit and others, detailed below is how to execute the TargetScan pipeline using perl scripts downloadable from http://www.targetscan.org/cgi-bin/targetscan/data_download.vert72.cgi (bottom of page has executable scripts). 

**Step 1: Run targetscan_70.pl**

The script takes two input files:
1.  A tab-delimited file that lists the miRNA seed sequences and the species in which they are present.
2.  A tab-delimited multiple sequence alignment of the 3' UTRs of genes from the desired species.

To generate file 1, download the miR Family info file from TargetScan, and subset it for human miRNA and remove duplicate miRNA entries: 
``` bash
wget http://www.targetscan.org/vert_72/vert_72_data_download/miR_Family_Info.txt.zip && unzip miR_Family_Info.txt.zip && rm miR_Family_Info.txt.zip

grep "hsa" miR_Family_Info.txt | awk '{print $1, $2, $3}' | sort -u -k1 | tr ' ' '\t' > hsa_miR_Family_Info.txt
```

For now, the UTR file has been manually created as a tab separated file:
``` bash
hsa_circ_0022393        9606    GGACCUGGCCUGGGCCGUCAGCUACUACAUCCGGUUCUUCAUCACCUACAUCCCUUUCUACGGCAUCCUGGGAGCCCUCCUUUUCCUCAACUUCAUCAGGUGCCUGGGCUUUGCUGAUUCAUCCCUGGGCCCUCCACUGGCACUGAUGAUGGCUUUGGCAUGAGAAGAGGCUCGGACGGCUCCACUUUGCCUGGGGACCUCCCAUCCUGGCCCCUUGGCAGGGGCUGUCCUGGGCUGGUUGGGAGAGUGGAUCUGUUCCCACUGUGGCUGGGCCCCCGGAGGCCCCUGCGCUGAGCUGUGUUCUCUUGCAGGUUCCUGGAGAGCCACUGGUUUGUGUGGGUCACACAGAUGAAUCACAUCGUCAUGGAGAUUGACCAGGAGGCCUACCGUGACUGGUUCAGUAGCCAGGUAGGGAAGUCAGGGCCGGUCACCAGAGCCUGGUCCCAAUGUCCAUGUCCUGGCCCCAGAUGAUCUGCACCUGCUAACACAGGGCUGUGCUGGGCCUCCCGGGGCUGCCUGUGUCCUUUUGCUGGGAGCUUCAGGGCUGAGCAAGAAAGGGCAGCAGGAAUCCCCCAGAGGGACAGGUGGGGCAUCUGGGUGGGCUGAGGCCAUCAGGCAGGACGGUAUGAUGUGGACGGGUGGCCUGGAGACCCUGGGUAAGGCUGGGCCCCCUGGGAAUGAGGCCGGGCCCUUGGGCCUUCCUUUGUUCCCUGACACUCAUCCCCUCCACUGACAGCUGACAGCCACCUGCAACGUGGAGCAGUCCUUCUUCAACGACUGGUUCAGUGGACACCUUAACUUCCAGAUUGAGCACCA
```

Run the command: 
```bash 
../../scripts/targetscan_exe/targetscan_70.pl ../../data/TargetScan/hsa_miR_Family_Info.txt ../../data/TargetScan/hsa_circ_0022392.txt targetscan_70_output.txt
```

**Step 2: Run targetscan_70_BL_bins.pl & targetscan_70_BL_PCT.pl**

Two Perl scripts calculates conserved branch length (BL) and probability of conserved targeting (PCT).

Input files for `targetscan_70_BL_bins.pl`: 
1. a UTR file = a tab-delimited multiple sequence alignment of the 3' UTRs of genes from the desired species
		which is the same as the input file for `targetscan_70.pl`.
2. The script also requires a directory called "PCT_parameters" containing the following file: `Tree.generic.txt`.
This is a tree file describing the phylogeny of all TargetScanHuman 7 species as defined by the conservation of their UTRs.

The script must be run in the same directory as `PCT_parameters`, you can write the output files to a results dir. 
```bash
../../scripts/targetscan_exe/targetscan_70_BL_bins.pl hsa_circ_0022392.txt > ../../results/circRNA_TargetScan/UTRs_median_BLs_bins.my_output.txt
```

Input files for `targetscan_70_BL_PCT.pl`:
1.  miRNA file, the same miRNA file used by `targetscan_70.pl`.
2.  Predicted targets, the output file generated by `targetscan_70.pl` (`targetscan_70_output.txt`)
3.  UTR bin info, the output file generated by `targetscan_70_BL_bins.pl` (`UTRs_median_BLs_bins.my_output.txt`)

The script must be run in the same directory as `PCT_parameters`, you can write the output files to a results dir. 
```bash
../../scripts/targetscan_exe/targetscan_70_BL_PCT.pl hsa_miR_Family_Info.txt ../../results/circRNA_TargetScan/targetscan_70_output.txt ../../results/circRNA_TargetScan/UTRs_median_BLs_bins.my_output.txt > ../../results/circRNA_TargetScan/targetscan_70_output.BL_PCT.my_output.txt
```

**Step 3: run targetscan_70_contect_scores.pl**

Perl script `targetscan_70_context_scores.pl` calculates context++ scores for a set of miRNA targets predicted by `targetscan_70.pl` and further processed by `targetscan_70_BL_PCT.pl`.  Both of these steps must precede context++ score calculation.

The script takes several input files:
*specified in command*:
1.  a miRNA file: a tab-delimited text file of mature miRNA information.  This is different from the file required by `targetscan_70.pl`, where we collapsed replicate miRNA entries in field 1. This contains full miRNA sequences which have variation so do not collapse the information. Generate the file using `miR_Family_Info.txt` that we downloaded in step 1:
```bash
grep "hsa" miR_Family_Info.txt | cut -f1,3,4,5 | tr ' ' '\t' > miR_for_context_scores.txt
```
2.  a UTR file: a tab-delimited multiple sequence alignment of the 3' UTRs of genes from the desired species which is the same as the input file for `targetscan_70.pl`.
3. a predicted targets file with BLSs and PCTs: output from `targetscan_70_BL_PCT.pl` (`targetscan_70_output.BL_PCT.my_output.txt`)
4. ORF lengths file (for ORFs matching 3' UTRs in UTR_file): contains the length of each ORF corresponding to aligned 3' UTRs. Is outputted in the step below.
```bash
hsa_circ_0022393	9606	816
```
5. ORF 8mer counts file: contains the number of 8mer sites in ORFs of file (4). You generate this file using a perl script `targetscan_count_8mers.pl`:
```bash
./targetscan_count_8mers.pl hsa_miR_Family_Info.txt hsa_circ_0022392.txt > hsa_circ_0022392_8mer_counts.txt
```
Gives 8mer counts and sequence length file. 

*with names hard-coded in the analysis script*:
6. "TA_SPS_by_seed_region.txt": contains TA and SPS parameters for each seed region
7. "Agarwal_2015_parameters.txt": with model parameters to calculate context++ score contributions
8. "All_cell_lines.AIRs.txt": AIRs for each region of each 3' UTR. Run the following code to download the full file and grab the first 4 columns:
```bash
wget http://www.targetscan.org/vert_70/vert_70_data_download/3Pseq_tags_AIRs.zip && unzip 3Pseq_tags_AIRs.zip && awk '{print $1,$2,$3,$4}' All_cell_lines.AIRs.txt | tr ' ' '\t' > All_cell_lines.AIRs_tmp.txt && rm All_cell_lines.AIRs.txt && mv All_cell_lines.AIRs_tmp.txt All_cell_lines.AIRs.txt
```

Run the command
```bash
../../scripts/targetscan_exe/targetscan_70_context_scores.pl miR_for_context_scores.txt hsa_circ_0022392.txt ../../results/circRNA_TargetScan/targetscan_70_output.BL_PCT.my_output.txt hsa_circ_0022392.lengths.txt hsa_circ_0022392_8mer_counts.txt ../../results/circRNA_TargetScan/targetscan_70_context_scores_output.txt
```



