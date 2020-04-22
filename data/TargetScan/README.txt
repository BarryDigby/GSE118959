All 3P-seq tags (score = raw counts) per cluster (in tags_by_cell_line/)
	HEK293.all_clusters.bed
	HeLa.all_clusters.bed
	Huh7.all_clusters.bed
	IMR90.all_clusters.bed

3P-seq tags (score = normalized counts) per cluster that overlap 3' UTRs (in tags_by_cell_line/)
	HeLa.3UTR_overlap_clusters.bed
	Huh7.3UTR_overlap_clusters.bed
	IMR90.3UTR_overlap_clusters.bed
	HEK293.3UTR_overlap_clusters.bed

All_cell_lines.tags.bed: sum of normalized 3P-seq tags across all cell lines
All_cell_lines.tags+pseudocounts.bed: combined normalized 3P-seq tags and pseudocounts 
                                      added at end of annotated UTR (5) and at end of extended 3' UTR (0.01)
All_cell_lines.AIRs.txt: AIR by UTR and genome region with the following fields:
	UTR ID
	UTR region start for this AIR
	UTR region end for this AIR
	AIR (range of 0 to 100)
	total normalized 3P-seq tags and pseudocounts for this UTR
	strand
	hg19 genome chromosome for this AIR
	hg19 genome region start for this AIR
	hg19 genome region end for this AIR
