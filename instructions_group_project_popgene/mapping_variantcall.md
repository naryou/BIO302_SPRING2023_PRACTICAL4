# Mapping to a reference genome


To map our sequencing reads (.fastq.gz) to the reference genome we will use BWA. This program and Bowtie are the most common mapping programs to use. In Galaxy, this tools is called BWA-MEM. Outside galaxy, one usually creates an index of the reference genome with the function 'bwa index ref_genome.fasta' that creates a file called 'ref_genome.fasta.fai'. Fortunately, in Galaxy the index of the reference genome is done in parallel of the mapping.

In the tool parameters of BWA-MEM we will indicate that we will 'Use a genome from history and build and index' and select the reference genome that we uploaded (pve_haplotypeT_chromosome_1.fasta). To construct the index we will let BWA decide which is the best algorithm. We will also indicate that our sequencing reads are paired-end (*_R1.fastq.gz and *_R2.fastq.gz). We will also specify the analysis mode indicating that we are using ILlumina reads. We will indicate that we want the results to be sorted by chromosomal coordinates. The resulting file of the mapping is a '*.bam' file, but we need to indicate how do we want it to be sorted. You get one '*.bam' file per individual.

Remember to change the names of the resulting files by clicking on the 'Edit attributes' and change the name to something more meaning full. For example, instead of 'Map with BWA-MEM on data XX, data XX, and data XX (mapped reads in BAM format)' change it to 'popXX_indXX.bam' so that you can keep track which individual comes from which population.

We will repeat the same pipeline per individuals. To make it easier to keep track we will create a history of mapping per population. For instance, we can call the first History 'mapping to reference pop TR' and map all six individuals of the same population. Then you can copy the same pipeline and use it for populations (SK, CH, and EN).


# Variant calling

For the variant calling we will use a program called bcftools because its simple pipeline (and because is in Galaxy). Other common programs to make the variant calling are GATK and FreeBayes. In Galaxy, we will use the tool 'bcftools mpileup' first. This tool will xxxxxxx.

The variant calling can be done one individual at at time with the option 'Single BAM/CRAM' or 'Multiple BAM/CRAMs' in the Alignment Inputs part. In this case, let's run  'bcftools mpileup' per population (i.e., add all six individuals of a single population). We will leave all the option in default. In the 'Output options' we will select DP, AD, and SP (What does these parameter means?) and indicate that we want a compressed BCF (What is the difference of BCF and VCF?).

This will create a binary compressed version of a VCF file (format bcf). Remember to rename this file into something meaninging full such as 'popXX.bcf'.

The step above (piling up the sequecing reads) is only the first step to the variant calling. Now we are going to make the actual calling of variants (SNPs, insertions, and deletions) withe the tool 'bcftools call'. Again, there are many options but today we are going to leave eveything as default. For this part we will specify that the output is in 'uncompressed VCF' format.

We will do this for every population, meaning that in the end you should end up with four different VCF files with six individuals each (ideally with different names each).

Now we are going to use the tool 'bcftools merge' to merge all four VCF files into a single one in order to proceeed with the quality filtering of the VCF. Again, we will indicate that we want an 'uncompressed VCF' as output.
