# Group project for BIO302
## Population genomic consequences of the range expansion and the transition to self-fertilization in of Primula vulgaris

Primula vulgaris aka the common primrose is a very charismatic plant that flowers around this time of the year (March-May). This plant has also played 
a key role as a model for ecological and evolutionary work since the 1940's. It is also quite pretty. This plant has a wide distribution from Asia (the Caucacus region) to
all the British Isles. You can find it in almost all continental Europe and in many mediterranean Islands. It has been suggested that P. vulagris
and its most closely related species originated in Asia and expanded its distribution from to the west during the Quaternary. This expansion
could have lead to the accumulation of genomic consequences.

Primroses have also been a model to study the consequences of the transition from outcrossing to self-fertilization. This mating system transition
is very common in plants and can also lead to consequences at the genomic level (Wright et al., 2007; 2013). In P. vulgaris, this transition is represented by the shift from
heterostyly to homostyly (Li et al., 2016; Huu et al., 2016). Heterostyly is a floral polymorphism that promotes outcrossing that has been reported in 
at least 28 families of plants. Conversely, homostyly is associated with self-fertilization. This transistion has occured in many species of Primula. 
In this species, the transition to homostyly has been reported in England (Mora-Carrera et al., 2021).

Aims of the project:
- Determine what are the genomic consequences of the recent range expansion of Primula vulgaris
- Determine what are the genomic consequences of the shift from heterotyly to homostyly in P. vulgaris

What do you have:
```
-Paired-end sequences (*.fastq.gz) from 24 individuals from four populations of Primula vulgaris
  -Turkey (TR), Swizterland (CH), England 1 (EN1), and England 2 (EN2).
  -All populations are heterostylous except from EN2.
-A reference genome (*.fasta) from Primula veris 
```

What are you expected to do in this project:
```
- Mapping to the reference genome with BWA
	-Create an index
	-Map to reference
- Make the variant calling with bcftools
- Filter the VCF file with bcftools
- Do a STRUCTURE and PCA analyses
- Perform analyses the genetic diversity of the populations

If there is time (and you are interested):
- Determine what is the genetic basis of the transition from heterostyly to homostyly

Ask me questions about the analyses and the literature.

```

For the mapping and variant calling you can look into this tutorial (https://github.com/EmilianoMora/BIO302_SPRING2023_PRACTICAL4/blob/main/instructions_group_project_popgene/mapping_variantcall.md). For the VCF filtering and population genetic analyses you can look in the tutorial of practical session 4 (https://github.com/EmilianoMora/BIO302_SPRING2023_PRACTICAL4/blob/main/Instructions_practical4.md).

Key papers to read:
```
- Li et al. 2016. Genetic architecture and evolution of the S locus supergene in Primula vulgaris. Nature Plants 12:1-7
- Wright et al. 2013. Evolutionary consequences of self-fertilization in plants. Proceedings of the Royal Society B: Biological Sciences. 280: 20130133
- Wright et al. 2008. Genomic Consequences of Outcrossing and Selfing in Plants. 
- Laenen et al. 2018. Demography and mating system shape the genome-wide impact of purifying selection in Arabis alpina. 
- Bachmann et al. 2019. Genetic basis and timing of a major mating system shift in Capsella. New Phytologist, 224: 505-517 
- Sicard A, Lenhard M. 2011. The selfing syndrome: a model for studying the genetic and evolutionary basis of morphological adaptation in plants. Annals of Botany. 107:1433–43.
- González-Martínez et al. 2017. Range Expansion Compromises Adaptive Evolution in an Outcrossing Plant. Current Biology
- Huu et al. 2016. Presence versus absence of CYP734A50 underlies the style-length dimorphism in primroses. eLife
- Stubbs et al. 2022. Whole-genome analyses disentangle reticulate evolution of primroses in a biodiversity hotspot. New Phytologist.
- Mora-Carrera et al. 2021. Different molecular changes underlie the same phenotypic transition: Origins and consequences of independent shifts to homostyly within species. Molecular Ecology
```
