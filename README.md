# Lynch-PoolSeq-estimator
Maximum likelihood estimator of population allele frequencies from Pool-seq data using the method of Lynch et al. GBE.2014 

This python script applies to the allele frequency estimator for pool seq data from Lynch et al. 2014 GBE.

Pysam and Numpy must be installed in the version of python used. Bam files must be indexed.

To guard against mis-mappings and CNVs, the program only outputs sites with coverage between a third and twice the mean at the set of SNPs considered

The SNP file must have the tab seperated fields in the following order: chromosome, position (one-based), reference allele, alternate allele
If you want to change things like mapping and base quality threshold, edit the python code under the section "Input arguments"

Usage is ./LynchPool_called.py bamfile filenameout target_SNP_file reference_genome

