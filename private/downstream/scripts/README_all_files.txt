The all_lines analysis groups samples per cell line, cat the fastq files and processes them together with Nanocompore.

Results are given as:

*_modified_sites.xls <- this file contains:
				all_significant_sites Sheet <- all sites found to pass LOR and GMM pval thresholds (sites on IVT or ORF junctions are NOT excluded)
				all_canonical_Us_significant_sites Sheet <- all sites found to pass LOR and GMM pval thresholds with a central U in the reference kmer and chosen in canonical transcripts (sites on IVT or ORF jucntions are excluded)
				5p_significant_sites Sheet <- sites with genomic position <100 found to pass LOR and GMM pval thresholds (sites on IVT or ORF junctions are NOT excluded)
				5p_canonical_Us_significant_sites Sheet <- sites with genomic position <100 found to pass LOR and GMM pval thresholds  with a central U in the reference kmer and chosen in canonical transcripts (sites on IVT or ORF junctions are NOT excluded)
				Burrows_redundant_ALL_sites Sheet <- redundant Burrows sites significant and non-significant (sites on IVT or ORF junctions are NOT excluded)
				Burrows_redundant_significant_sites <- redundant Burrows significant sites (sites on IVT or ORF junctions are NOT excluded)
				Burrows_non_redundant_CandNC_significant_sites <- non redundant Burrows significant sites in canonical and non canonical transcripts (sites on IVT or ORF junctions are NOT excluded)
				Burrows_non_redundant_CandNC_all_sites <- non redundant Burrows ALL sites in canonical and non canonical transcripts (sites on IVT or ORF junctions are NOT excluded)

*_plots_per_transcript.pdf -> this file contains the sharkfin plots for each transcript model of the assembly in which every point represent a position. Points labeled with their letter sequence are those which have a LOR and a pvalue over the set threshold in all datasets of the cell line and that DO NOT fall in an IVT or ORF junction. Colors represent RAPID fragments.

*_plots_per_transcript_5p.pdf -> same as the *_plots_per_transcript.pdf but focused ONLY on sites with genomic Position <= 100. 


