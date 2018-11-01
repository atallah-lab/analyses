# analyses

#### gene_counts_all
* Directory with R data files storing counts of L1 target loci throughout each protein-coding Ensembl gene.

#### gene_counts_exon
* Directory with R data files storing counts of L1 target loci in exonic ranges of each protein-coding Ensembl gene.

#### sim_out_analyses
* Contains analyses of simulation results.

#### cancer_type_site_counts.ipynb
* Notebook to estimate the probability of driver, passenger, and null (no effect) L1 insertions across lung, colon, and brain cancer.

#### check_target_site_map.ipynb
* Notebook for checking that the L1 target site annotation data is valid.

#### compare_gene_pd_tables.ipynb
* Notebook for comparing the probability of L1 insertion in each Ensembl gene when considering only exonic insertions, or insertions across the whole gene.

#### count_all_sites_for_each_gene.ipynb
* Counts the number of L1 target loci of each snap-velcro category across the whole range of each Ensembl gene.

#### count_cds_sites_for_each_chromosome.ipynb
* Notebook that uses the hg38 gff3 file to compute the number of L1 target sites of each SV category in known coding regions and non-coding regions of the genome (separately), and also the coding fraction of each chromosome.

#### count_exon_sites_for_each_chromosome.ipynb
* Equivalent to count_cds_sites_for_each_chromosome.ipynb, except for exonic regions rather than just coding regions.

#### count_exon_sites_for_each_gene.ipynb
* Uses the table exann.rda (output by sim-develop/src/sim_exon_annotation.ipynb) to determine the number of sites in the exonic regions of each gene in Ensembl v86.

#### count_sites_exonicvsnon.ipynb
* Notebook to count the distribution of L1 target loci in exonic vs. non-exonic regions of the human genome (hg38).

#### count_sites_genevsnon.ipynb
* Notebook to count the distribution of L1 target loci in gene vs. non-gene regions of the human genome (hg38).

#### make_gene_pd_table_from_counts.ipynb
* Notebook for computing probability of insertion for each Ensembl gene, using the counts computed by count_all_sites_for_each_gene.ipynb and count_exon_sites_for_each_gene.ipynb.
