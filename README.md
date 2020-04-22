# GSE118959
Analysis of GSE118959 circRNAs in LNCaP cell lines
***

# Analysis

#### Differentially expressed circRNAs
Run `DE_circRNA.Rmd` script to output results for DE circRNAs in Clone1 vs control & Clone9 vs control. 

To get the overlapping circRNAs present in upregulated/downregulated sets, use the following bash command:

```bash
awk -F'\t' 'NR==FNR{c[$2$2]++;next};c[$2$2] > 0' clone1_upreg.txt clone9_upreg.txt > up_regulated_overlap.txt

awk -F'\t' 'NR==FNR{c[$2$2]++;next};c[$2$2] > 0' clone1_downreg.txt clone9_downreg.txt > down_regulated_overlap.txt
```

The above code returns matching lines using column 2 as the key. As this analysis will be an unweighted graph, there is no need to worry about which Fold Change value to take, I am simply interested in the circRNA ID and its genomic coordiantes. 
