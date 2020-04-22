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
sudo cpanm Bio::TreeIO --force
sudo cpanm Statistics::Lite
```



I am still testing the software. Once it is running, I will have to develop a pipeline to extract DNA sequences, convert to RNA and then format them for targetscan input. (use samtools faidx here, also `sed 's/T/U/g'` for DNA to RNA conversion). 
