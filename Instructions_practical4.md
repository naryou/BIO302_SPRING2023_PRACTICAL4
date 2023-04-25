# VCF filtering

In your second practical session you did a variant calling. Just as a reminder you usually take some sequencing reads (in .fastq format) and map them into a reference genome (.fasta) using BWA, bowtie or any other alignment program. The resulting file after the mapping is a file in either SAM or BAM format. With these files you can do the variant calling with GATK, FreeBayes, or BCFTOOLS. The resulting file after the variant calling is a Variant Calling File (VCF). In this file you have all the information about invariable and variable sites and its associated quality values.

Often after the variant calling we end up with a very messy and long VCF file that is full of sequencing errors, genotyping errors during variant calling, some invariant sites that tell nothing about diversity and so on. So we need to do some cleaning of the VCF file to keep high quality sequencing sites so that we can perform some analyses.

For this paer of the practical you will use the file called 'BIO302_practical_filtering.vcf.gz' that I sent you in the email. You can download the file [here](https://wetransfer.com/downloads/f95a2d784ff8c6eb2b43c8fcccb3ac9020230424070021/c59faa7fd8ec72778f2551a7543b619120230424070041/1febb0).

Upload this file into your Galaxy account.

After you had upload your vcf file, take a look on how many sites we have in the VCF file so that we can compare with what we end up.

If you inspect your VCF file all the way down where the genetic information starts, we can see that there are several columns. In the INFO column you will see a very long and confusing set of arguments each separated by a semicolon (;). These is the infromarion of the quality of each site. In this INFO column AC, AF, AN, BaseQRankSum, etc says something meaningful of each site. If you want to know what does each of these arguments mean, you can go to the header of the VCF file to the part that says '##INFO= ...' the meaning of the all arguments in the INFO column can be found in here. The information on this INFO column will allow us to make the filtering.

For the filtering today we will use 'bcftools filter'. We could also use GATK and vcftools that work in a very similar way.

In the 'bcftools filter' we will indicate the VCF file that we want to filter. Then we click on the 'Restrict to' option and move to where the Include and Exclude options are. Now we have two ways to set up the filter. We can tell bcftools to include only sites that have certain threshold (Include) or we can do the opposite (i.e., exclude sites that are below certain threshold) with the Exclude option.

The filters that we are going apply, using the 'Include' option, are the following:

```
QD > 2.0
FS < 60.0 
MQ > 40.0 
MQRankSum > -12.5 
ReadPosRankSum > -8.0 
FORMAT/DP > 10 
FORMAT/DP < 400 QUAL > 100
```

If you want to know more about what do these and other argument means, you can look [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants).

**Hint:** If you want to use the 'Exclude' option you will just need to flip the '> or <' around.

Now, in the command line version of 'bcftools filter' one could specify all the thresholds in a single line command. However, I don't know how to do it in Galaxy (duh!!). So you will have to run the tools for each filtered argument, one by one (emoji of 'I am sorry' face). In the output specify the 'uncompressed VCF' so that you can take a quick look at the file after each filter. Remember to change use the output VCF file after each filter for the next one. Also change the name of the file everytime you run a filter so that you know which filter you used in each run. For example:

```
original_unfiltered.vcf

After the 'QD > 2.0' filter you can call it

original_unfiltered_QD.vcf

After the 'FS < 60.0' filter you can call it

original_unfiltered_QD_FS.vcf
```

And so on and so forth. After you apply all the filters, you will end up a VCF file with super long name. It is a bit annoying but it helps you keep track on what you are doing. If you see the header of the VCF file, you can see that everytime you add a new filter it gets registered. This allows us to see what filter and which commands have been runned on the VCF files.

After using the 'bcftools filter' now we are going to move and use the 'bcftools view' tool. It looks very similar but has way much more options. We are going to click in the 'Filter' options and specify that we want a maximum allele per site of '2' in the 'Max Alleles' option (Why is this?). We are also going to click on 'Select Types' box and only keep 'snps' (Hint: alternatively, you could also 'Exclude Types' box and eliminate indels, mnps, and other). Remember to rename that last VCF file with all the information of the filters.

Again we specify 'uncompressed VCF' as output.

Now take a look at the number of lines of the in the very last VCF (the super amazingly filtered one) and the one in the very beginning. How many sites have been lost on the way?

If you finished all the filters you can now rename the VCF file with something more short such as 'original_filtered.vcf'.

With this VCF file we can now do a lot of population genetic analyses.

**Note:** Filtering is one of the most important parts because bad filtering can lead to bad results. Sometimes is hard to know how hard do we need to be in our filtering. IF we filter TOO much, we might end up with almost no sites for our population genetic analyses. However, very light filtering could give us very misleading results.

# Analzying population structure with: STRUCTURE

**Background**
To investigate the genetic structure of a population from genotype data we use one of the most widely used software called [Structure](https://web.stanford.edu/group/pritchardlab/structure.html). Structure uses a Bayesian algorithm to group the individuals based on similarities/differences in their allele frequencies. The user sets the number of *assumed* groups or clusters (K) as a prior and Structure tries to group the individuals in those K clusters while maximizing the likelihood of shared genetic variation in each group. It's common to run the software for a range of pre-selected clusters (for instance K from 1 to 5) and then compare the results. There are also available methods to compare the different clsuters and estimate the **Best K** which is the most likely number of genetic clusters in our data.  
You can read more about Structure in the overview article [here](https://www.frontiersin.org/articles/10.3389/fgene.2013.00098/full). *You can find this article also on OLAT.*

## Running Structure on Galaxy 

Download the structure file from [this link](URL by Emiliano) and upload it to your Galaxy history. Change the format to **tsv** (click on the pencil icon, convert, set **New Type** as **tsv** and save).
From Tools panel, choose [**Structure** using multi-locus genotype data to investigate population structure](https://usegalaxy.eu/root?tool_id=toolshed.g2.bx.psu.edu/repos/iuc/structure/structure/2.3.4+galaxy1)

Set the following parameters:


**Genotype data**: SNPs_mac_TR_SK_CH_pve001_var_832_3pop.tsv

**Number of runs**: 2

**Number of populations assumed** or [K]: 1

**Length of burnin period**: 50000

**Number of MCMC reps after burnin**: 200000

**Number of diploid individuals in data file**: 30

**Number of loci in data file**: 832

**Input file contains individual labels**: Yes

**Input file contains a user-defined population-of-origin for each individual**: Yes

**Data file contains row of marker names**: No

In `extraparams` section, set:

**Use POPDATA for location information:** No


and click on `Run Tool`

This generates four outputs:
- run_K_1.out
- run_K_1.log
- run_K_1.mainparams
- run_K_1.extraparams

When running Structure, we need to assume a range of populations or clusters (**K**) and later find out which number of clusters is best explaining the genetic variation in our data. Therefore, you need to repeat the Structure run again, each time for a different K value (*Number of populations assumed*): K=2, 3, 4 and 5. All other parameters are kept the same as for K=1.

## Visualizing the Structure results
Download the run_K_x.out files to your computer. If they are zipped, unzip them. You should have the following files on your computer:
```
K_1_run_1_f
K_1_run_2_f
K_2_run_1_f
K_2_run_2_f
K_3_run_1_f
K_3_run_2_f
K_4_run_1_f
K_4_run_2_f
K_5_run_1_f
K_5_run_2_f
```

Now zip all the files together (e.g., select all files and choose 'add to run.zip'). Go to the [CLUMPAK webpage](http://clumpak.tau.ac.il/) and upload the zipped file in **Upload zip file containing STRUCTURE/ADMIXTURE runs**. Give your email address and click on **Submit Form**. This will take a few minutes to run and will give you the results in a PDF file and zip summary. Download them to your computer. Inspect the plots and discuss them among yourselves. There are 5 plots, representing K 1-5. Each horizental bar is one individual and each color is one genetic cluster. What do you understand from the plots? 

## Estimating the best K

Now we go back to Clumpak webpage and choose [**Best K**](http://clumpak.tau.ac.il/bestK.html) from the top banner.
Upload the zip file again, give your email address and click on 'Submit Form'. It will take a few minutes to run and will give the results in PNG files and a summary file. Download them to your computer.

What is the best K estimated based on Evanno method? Does it differ from the estimated best K from the Ln Prob. method?

# Population genomic analyses with Rstudio

Galaxy has (apparently) many tools to analyze genetic structure and estimate genetic diveristy. However, some of these tools have been outdated so we cannot use them. For that reason we will use R to make some of the basic analyses. In reality we do these analyses in the command line using (mostly) UNIX or python languages. Fortunately, you will not need have to do any UNIX stuff. However, if you are curious on how this looks like you can ask me or any of the other TAs so that you see how each tool looks like in the command line.

# Analzying population structure with: Principal Components Analyses (PCA)

Now that you know how the STRUCTURE results looks like, we will compare with a the PCA. Both methods are clustering analyses. The main difference is that one (STRUCTURE) is a model-based method that works using some assumption about allelic frequencies of populations in an evolutionary framework (e.g., using Hardy-Weinberg equilibrium), whereas the other one (PCA) is a model-free analysis that has no underlying assumption. Usually both analyses give similar results. In papers of population genetics you usually see both of them next to each other.

For the PCA we will use the package 'vcfR' in Rstudio and a shortened version of the VCF that you used for the filtering (**merged_allpops_biallelic_mac_TR_SK_CH_pve001_var_1000.vcf**). You can download this file to your computer [here](https://github.com/EmilianoMora/BIO302_SPRING2023_PRACTICAL4/blob/main/merged_allpops_biallelic_mac_TR_SK_CH_pve001_var_1000.vcf).

For the analyses in Rstudio, we first will set up the directory where we have all the data in your computer by clicking in 'Session > Set Working Directory > Choose Directory ...'. Alternatively, you can use the following command:

```
setwd('~/where/ever/directory/your/data/is/')
```

Also, we will load the libraries that we will need by running the commands:

```
library(adegenet)
library(ape)
library(pegas)
library(seqinr)
library(ggplot2)
library(vcfR)
library(phytools)
```

I suggest you create an Rscript (Ctrl+Shift+N) so that you keep track of everything that you are doing today.

Now that we have the packages and the directory, we can start by reading the VCF file in R with the following command:

```
df <- read.vcfR('./merged_allpops_biallelic_mac_TR_SK_CH_pve001_var_300.vcf', limit = 1e+07, nrows = -1, skip = 0, cols = NULL,
convertNA = TRUE, verbose = TRUE)
```

Take a quick look on the file with

```
head(df)
```

You can take a look on the amount of missing data with the following commands:

```
gt_only <- extract.gt(df)
sum(is.na(gt_only[,1]))/length(gt_only[,1])*100  # gives the percentage of missing data in each column
```

Take a look at the gt_only. It is a more digested way to see at the results of a VCF. In this format the genotypes are shown as 0/0, 0/1, and 1/1 that mean homozygous for the reference, heterozygous, and homozygous for the alternative allele, respectively. We could also have a look at the average missing data per loci and per individual.

```
res_loci <- round(apply (gt_only, 1, function (i) sum (is.na(i))/length(i)), 2) # missing data for each loci
res_ind <- round(apply (gt_only, 2, function (i) sum (is.na(i))/length(i))*100, 2) # missing data for each individual
```

To do the PCA we need to change the format of the VCF file

```
df_genlight <- vcfR2genlight(df, n.cores = 7)

POPS <-c("RS170","RS170","RS170","RS170","RS170","RS170","RS170","RS170","RS170","RS170",
         "RS180","RS180","RS180","RS180","RS180","RS180","RS180","RS180","RS180","RS180",
         "RSBK01","RSBK01","RSBK01","RSBK01","RSBK01","RSBK01","RSBK01","RSBK01","RSBK01","RSBK01") # Here we create a vector with the population information
         
df_genlight@pop <- as.factor(POPS) #Here we add the population information (POPS) into the genlight file (df_genlight)
```

Now we can proceed to do the PCA with the following command:

```
pca1 <- glPca(df_genlight)

myCol <- s.class(pca1$scores,pop(df_genlight), col=c("blue", "red", "green"))
add.scatter.eig(pca1$eig,posi = "bottomright", 2,2,1)
```

You will need to specify the number of axes (eigenvalues) to use. One plot was created automatically, you can see that the very first eigenvalues explain a lot of the variance. The more axes we keep the more resulting we get in the PCA. In this case we indicate that we want to keep the first '10' eigenvalues.

**Question:** How do the results looks like? Do the populations overlap on the PCA? Do the results align with the results of the STRUCTURE analysis?


# Analyzing genetic diversity with PopGenome

There are many programs and packages to estimate genetic diversity of a population. Some programs/packages are better or more efficient than others. In R one can work with VCF files. However, VCF files cannot be use right away and they usually need to be trasnformed into a more efficient format.

Today we are going to use the package 'PopGenome'. It stopped been updated a couple of years ago but we can still get the package in the [archive](https://cran.r-project.org/src/contrib/Archive/PopGenome/). If you have not installed the package in Rstudio, you just need to download the file called 'PopGenome_2.7.5.tar.gz' and run the following command:

```
install.packages("~/DIRECTORY_YOU_HAVE_THE_FILE/PopGenome_2.7.5.tar.gz", repos=NULL, type="source") library(PopGenome_2.7.5.tar.gz) #To check if it installed correctly
```

Now we just need to input our VCF. For this part, we will use another VCF file that is very similar the one you have been working on but with some filters that I've done myself. This VCF file has 30 individuals from three populations (RS170=Turkey, RS180=Slovakia, and RSBK01=Swizterland). This file only include SNP data from the first 10mb of the first chromosome. In order to input the file you just need to run the following command:

```
GENOME.class <- readVCF("~/DIRECTORY_YOU_HAVE_THE_FILE/merged_allpops_biallelic_mac_TR_SK_CH_pve001_10mb.vcf.gz",numcols=10, tid="pve_haplotypeT_001", from=1, to=10000000, approx=FALSE, out="", parallel=FALSE, gffpath=FALSE)
```

Note that the program will require a binary-compressed VCF file (*vcf.gz) and its index (*vcf.gz.tbi) and that one needs to specify the chrosome (tid="pve_haplotypeT_001").

If you get the following error, you might need to reduce the 'numcols=100' option:

```
Error in 1:numusedcols : NA/NaN argument In addition: Warning message: In readVCF("merged_allpops_biallelic_mac_TR_SK_CH_pve001_10mb.vcf.gz", : NAs introduced by coercion
```

If you had no problems you should see the following:

```
| : | : | 100 % |===================================================| ;-)

Concatenate regions | : | : | 100 % |===================================================| ;-)
```

You can get a breif summary of the data with the following command:

```
get.sum.data(GENOME.class)
```
You can also take a look at specific information such as the number of sites (n.sites) or the number of biallelic sites (n.biallelic.sites) wiht the following commands:

```
GENOME.class@n.sitess#number of sites GENOME.class@n.biallelic.sites#number of biallelic sites
```

One could estimate the genetic diversity with the following command:

```
GENOME.class <- diversity.stats(GENOME.class) get.diversity(GENOME.class)[[1]] #with this command you see all different estimates of genetic diversity that can be estimated GENOME.class@nuc.diversity.within #with this command you only get the estimate of nucleotide diversity
```

Apparently the way nucleotide diversity is estimated in PopGenome still needs to be controlled by the number of sites analyzed. So in this case we need to divide the estimated 'nuc.diversity.within' by 10000000 (the number of sites used in our data).

```
GENOME.class@nuc.diversity.within/10000000
```

The resulting estimate is more realistic, usually nucleotide diverstiy ranges from (0.01-0.00001). Everything above or below is way too big (or too little).

Now, we just did the analysis as if all 30 individuals were part of the same population, but if you remember we have three populations (RS170=Turkey, RS180=Slovakia, and RSBK01=Swizterland). If you remember, the STRUCTURE and PCA analyses told us that we are working with three genetically different populations. So analyzing all of them as a single population can bias our results.

**Question:** How do you think different levels of genetic diversity among population can bias the results?

To asign individuals to population in our file (GENOME.class), we can use the following command:

```
GENOME.class <- set.populations(GENOME.class, list(c("RS170_WD04","RS170_WD06","RS170_WE03","RS170_WE05","RS170_WE07", "RS170_WE12","RS170_WG02","RS170_WG10","RS170_WG11","RS170_WH09"), c("RS180_WE07","RS180_WE08","RS180_WF01","RS180_WF02","RS180_WF03", "RS180_WF04","RS180_WF05","RS180_WF06","RS180_WF10","RS180_WF11"), c("RSBK01_WA03","RSBK01_WB03","RSBK01_WC07","RSBK01_WD05","RSBK01_WD07", "RSBK01_WD08","RSBK01_WD10","RSBK01_WD11","RSBK01_WF02","RSBK01_WH02")), diploid=TRUE) #basically, we asign the individuals to three different R vector
```

**Question:** Lets assume that RS180 and RSBK01 are from the same population, how would you set up the above command to have two populations instead of three?

Anyway, now that we asigned our individuals to the populations, we can estimate genetic diversity of each population with the following commands:

```
GENOME.class <- diversity.stats(GENOME.class) get.diversity(GENOME.class) [[1]] #shows results of genetic diversity in Population 1 get.diversity(GENOME.class) [[2]] #shows results of linkage diversity in Population 2 get.diversity(GENOME.class) [[3]] #shows results of linkage diversity in Population 3
```

Again, we need to account for the number of sites used in our VCF files so ...

```
GENOME.class@nuc.diversity.within/10000000 #divide by the number of sites included
```

**Question:** Are there differences of genetic diversity among populations? Which population has the highest genetic diversity?

We can also run another set of population genetic tests that allows us to gain some information about the evolutionary history of the population. You can run the following commands:

```
GENOME.class <- neutrality.stats(GENOME.class, FAST=TRUE) get.neutrality(GENOME.class)[[1]] #shows results of neutrality tests in Population 1 get.neutrality(GENOME.class)[[2]] #shows results of neutrality tests in Population 2 get.neutrality(GENOME.class)[[3]] #shows results of neutrality tests in Population 3
```

Again, take a look at which population has the higher number of segregating (i.e., polymorphic sites). What does this mean?

Another test I didnt mentioned in the class was Tajima's D. The program also gives you a result for this test. Do you see any difference among populations? What does these results mean? Go on wikipedia and read about how to interpret [Tajima's D](https://en.wikipedia.org/wiki/Tajima%27s_D) (this will take you 5 min). Based on the results of these population, can you say something about them?

We can also estimate the extent of linkage disequilibrium (LD). Based on the previous results and what we saw on class about genetic bottlenecks, can you make a prediction of the extent of LD that we expect on these three populations?

```
GENOME.class <- linkage.stats(GENOME.class) get.linkage(GENOME.class)[[1]] #shows results of recombination tests in Population 1 get.linkage(GENOME.class)[[2]] #shows results of recombination tests in Population 2 get.linkage(GENOME.class)[[3]] #shows results of recombination tests in Population 3
```

There is a chance that this not run if you don't have enough memory on your computer. What we can do is to reduce the dataset to help our poor computers to run this commands.

```
GENOME.class_small <- readVCF("merged_allpops_biallelic_mac_TR_SK_CH_pve001_10mb.vcf.gz",numcols=10, tid="pve_haplotypeT_001", from=1, to=5000000, approx=FALSE, out="", parallel=FALSE, gffpath=FALSE) #Here we are just reducing the amount of site used from 10 million to 1 million.

GENOME.class_small <- set.populations(GENOME.class_small, list(c("RS170_WD04","RS170_WD06","RS170_WE03","RS170_WE05","RS170_WE07", "RS170_WE12","RS170_WG02","RS170_WG10","RS170_WG11","RS170_WH09"), c("RS180_WE07","RS180_WE08","RS180_WF01","RS180_WF02","RS180_WF03", "RS180_WF04","RS180_WF05","RS180_WF06","RS180_WF10","RS180_WF11"), c("RSBK01_WA03","RSBK01_WB03","RSBK01_WC07","RSBK01_WD05","RSBK01_WD07", "RSBK01_WD08","RSBK01_WD10","RSBK01_WD11","RSBK01_WF02","RSBK01_WH02")), diploid=TRUE) # we need to asign individuals again

GENOME.class_small <- linkage.stats(GENOME.class_small) get.linkage(GENOME.class_small)[[1]] #shows results of LD tests in Population 1 get.linkage(GENOME.class_small)[[2]] #shows results of LD tests in Population 2 get.linkage(GENOME.class_small)[[3]] #shows results of LD tests in Population 3
```

Now it should have runned. If it didn't, your computer might be too weak for the job.

Anyways, we can also estimate the genetic differentiation (Fst) of the three populations with the following commands:

```
GENOME.class <- F_ST.stats(GENOME.class) get.F_ST(GENOME.class) #Show all estimates of Fst
```

As you can see, the estimates of genetic differentiation can be quite different from each other. This is because the way we measure genetic differentiation can be biased by the genetic diveristy of the population. For instance Fst can be inflated by low (or high) genetic diveristy of the population, whereas Gst controls for the genetic diversity of each population. See the following results:

```
GENOME.class@nucleotide.F_ST #if you want to keep one specific measurement of Fst GENOME.class@Nei.G_ST #if you want to keep one specific measurement of Gst
```

Another factor that can bias our over-all estimates of genetic diversity is how genetic differentiation changes along the genome. As we saw in class, different genomic regions can have differences in mutarion rate, recombination rate, selection, etc. Also Fst might change between coding and non-coding regions (intronic and intergenic) All this might lead to regions with higher and lower Fst. The message is that collapsing everything together might give an erroneous estimate of Fst. So now we are going to estimate Fst along the genome using sliding windows.

To do these analysis, we will tell PopGenome to split the data in windows of 1000bp (1kb)

```
GENOME.class.slide_1kb <- sliding.window.transform(object = GENOME.class, width = 1000, jump = 1000, type=2)
```

Now we can tell the package to estimate Fst on each of these windows with the following commands:

```
slide_1kb <- diversity.stats(GENOME.class.slide_1kb) #First we need some estimates of genetic diveristy before doing the Fst slide_1kb <- F_ST.stats(slide_1kb, mode="nucleotide") #Now here we tell R to estimate Fst on every window slide_1kb@nuc.F_ST.pairwise
```

The results should look like a pairwise table among the three populations and per window. However, sisnce we have too many windows this looks like a mess. So we run the following command to make it look prettier:

```
pairwise.FST_1kb <- t(slide_1kb@nuc.F_ST.pairwise) pairwise.FST_1kb #Now take a look each row is a window and each column is the comparison each population against another one. #The 't(xxx)' in t(slide_1kb@nuc.F_ST.pairwise) is to transpose the table
```

Since we input 10 million bpfrom our VCF, and we have windows of 1000, then we should have 1000 windows (10000000 sites/1000 bp window size). We can check this with the following command:

```
length(GENOME.class.slide_1kb@region.names) #if you change the window size in 'GENOME.class.slide_1kb <- sliding.window.transform(object = GENOME.class, width = 1000, jump = 1000, type=2)' then thiscomamnd tell you how many windows you have.
```

Now we are going to plot the results with the following commands:

```
ids <- 1:10000 #adjust for the number of windows loess.fst_pop1_pop2_1kb <- loess(pairwise.FST_1kb[,1] ~ ids, span=0.05) loess.fst_pop1_pop3_1kb <- loess(pairwise.FST_1kb[,2] ~ ids, span=0.05) loess.fst_pop2_pop3_1kb <- loess(pairwise.FST_1kb[,3] ~ ids, span=0.05) plot(predict(loess.fst_pop1_pop2_1kb), type = "l", xaxt="n", xlab="Position (Mb)", ylab="Fst", main = "Chromosome 1 (1kb windows)", ylim=c(0,1)) lines(predict(loess.fst_pop1_pop2_1kb), col="blue") lines(predict(loess.fst_pop1_pop3_1kb), col="red") lines(predict(loess.fst_pop2_pop3_1kb), col="green") axis(1,c(1,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000), c("0","1","2","3","4","5","6","7","8","9","10"))
```

**Question:** Do you know what each of these commands is doing? Can you change the colors of the populations? Based on the results of genetic diversity and global Fst, what do you expect in this result? Which population is the most genetically differentiated?

Now that you know how to run a Fst analysis along the genome, do the same analysis to estimate genetic diveristy (nucleotide diversity) using windows.

# Simulating Genetic Drift and Natural selection with learnPopGen

Great! you reached the final part of the practical! I have no idea how long is this going to take you ideally you would ask stuff that you dont understand and so.

Now, in class we saw some of the evolutionary forces that shape the genome. If you remember, I made emphasis in two evolutionary forces, Natural selection and Genetic drift. Some population genetisis consider these two as the most important shaping the complexity of the genome. It sometimes is not easy to understand their effects when they are both working at the same time. A good way to wrap our head around these evolutionary forces is to visualize it. Now close your eyes and imagine Genetic Drift ... haha just kidding you cannot imagine Genetic Drift. What we are going to do is run some simulation using the R package learnPopGen. This package can be downloaded in R and runned through the command line in Rstudio. However, this package also has a shiny interface on the internet that makes everything more easy, so follow the link:

http://www.phytools.org/PopGen/

Now **click on drift-selection**. You will se the shiny interface. The results of the simulations looks very similar to some plots I showed earlier in class. In p0 we have the frequency of a biallelic allele at generation '0'. If p0=0.5 it means that both alleles segregate at equal frequencies.If you click on 'new plot' below you can see that everytime you run it, the results looks different, this is because genetic drift is a random process that generated different results, but do you see a pattern?

Yeah!! in all cases some alleles go into fixation (f(A) reaches 1.0) and sometimes it gets eliminated from the population (f(A) reaches 0.0), whereas in some populations both alleles keep segregating and different frequencies.

Now, what happens if we reduce the initial frequency of an allele? For example to p0=0.3 or p0=0.3. What is the pattern that you see? Remember that a new mutation in a population (a SNP, indel, duplication, etc) occurs initially in a single individual in the population ... what do you think its p0 frequency would be? Based on what you saw juts before, what is the most common fate of a new mutation in a population?

In class we mentioned that the strength of Genetic Drift is determined by the effective population size (Ne). Play around with the Ne parameter by increasing it, what do you see? Is there any pattern? What happen to the alleles when we increase or reduce Ne?

Now, play with Ne and p0, what happen there?

Now let's add Natural Seletion into the game. Here we have the fitness among all three genotypes of a biallelic loci. When we see 1.00 in each it means that all of them have the same fitness. If we change one fo the to 1.10 it means that individuals with this genotype have a 10% increase in fitness. IT can be reproductive fitness or increase change in survival. For now just play with the W(AA). What happens if you increase the fitness of this genotype?

Now go back and put W(AA) back to 1.00 and play around with the initial frequnecy p0=0.3. What is the pattern?

Now we said earlier that Genetic Drift will determine the efficacy of Natural Selection. The stronger the genetic drift the less selection is. To see this, play around with fitness of W(AA) and reduce Ne. What happens? You can also play around with p0. What would happen to a new mutation that increase fitness but is in a population with a small Ne? For instance, what happens in a mutation causes an increase of 80-90% in fitness?

With this we ahve reached the end of our journey into the exciting field of population genomics. You can now go home and have some rest.


