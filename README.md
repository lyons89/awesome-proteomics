# awesome-proteomics
An awesome list of proteomics tools and resources

I've been analyzing large scale proteomics data sets for over 5 years now. I've recently stumbled upon these awesome lists and wanted to make an in-depth one for proteomics.


## Proteomics 

### Table of Contents 
1. Learning Resources
2. Databases
3. Raw data search software/algorithms
4. Assorted Pipeline Tools
5. R Packages for Raw Data Analysis
6. R Packages for Stastical Analysis
7. R packages for Protein Pathway Enrichment
8. R Packages for Kinase Motif/Activity Analysis



## 1. Leaning Resources

[Ben Orsburn](https://proteomicsnews.blogspot.com/) has hands down the best protomeics blog I've ever seen. Ben is very knowledgable with a great sense of humor.

[Phil Wilmarth](https://github.com/pwilmart/Start_Here) has alot of goodlearning resources. He was python scripts for raw data analysis as well as blog posts detailing TMT data analysis techniques. 

[Review article](https://www.nature.com/articles/nrm1468) great proteomics review article

[Tutorial videos](https://www.youtube.com/channel/UC0v4sjdXLMa-OWR7IYeoFoA/videos) videos from NCQBCS, a project lead by the Coon lab. Contains lots of information regarding Experimental design, ionization, quantitative proteomics, analysis, post-translational modifications and more.

[ASMS video mass spec channel](https://vimeo.com/channels/asms) contains a lot of videos from leading researchers in the protoemics field. 

[Videos from Nikolai](https://www.youtube.com/c/NikolaiSlavovResearch/videos) most of his videos focus on single cell proteomics, and DIA

[Videos from Matthew Padula](https://www.youtube.com/c/MatthewPadula/videos) lots of great videos on the basis of mass spectrometry and proteomics.

## 2. Databases

Commonly used protoemics data repoitories. You can find raw data here, or publish your down for article publication requirements. 

[ProteomeXchange](http://www.proteomexchange.org/) is a global repository that contains links to all major databases, including MassIVE, Pride, iProX and more. Probably the best place to start.

[massIVE](https://massive.ucsd.edu/ProteoSAFe/static/massive.jsp)

[Pride database](https://www.ebi.ac.uk/pride/archive/)


This massive master list of databses from Pastel BioScience

http://www.pastelbioscience.co.uk/resources/databases.html - this list contains so much useful information that I'm sure the rest of my awesome list will be redundant

[Uniprot](https://www.uniprot.org/) - Complete protein database

[Biogrid](https://thebiogrid.org/) - protein-protein interaction database

[KEGG](https://www.genome.jp/kegg/) - biological pathway database

[Reactome](https://reactome.org/) - nicer looking biological pathway database


## 3. Raw data search software/algorithms

Software with a pretty GUI and a peptide search engine for raw file searches.


[Fragpipe](https://github.com/Nesvilab/FragPipe) currently my favorite raw file search software. It's much faster than maxquant and in my opinion has sleeker GUI. In a recent update it also allows searches to be set up on a linux server for even faster results. The software is very modular, it consists of [MSfragger](https://msfragger.nesvilab.org/) the database search algorithm, [Philosopher](https://philosopher.nesvilab.org/) that analyzes the database results, as well as others for PTM and TMT integration. 

[MaxQuant](https://www.maxquant.org/) is probably the most used and well known DDA software. Developed by Jurgon Cox, this completely free software is user friendly and is always being updated with new and original features. There is even a [youtube](https://www.youtube.com/c/MaxQuantChannel) that has tons of videos on how to use the software. 

[Peptide-shaker](https://compomics.github.io/projects/peptide-shaker) is like the swiss army knife of search tools. You can search data with multiple search engines inclduing, comet, tide, andromeda, mascot, X!Tandem and more that I've never heard of.  


## 4. Assorted pipeline Tools

2021 - [Monocle](https://github.com/gygilab/Monocle) a program built, in part by Deven Schweppe, for monoisotopic peak and accurate precursor m/z detection in shotgun proteomics experiments. I have only spent a little time tinkering around with this, but if I ever need to make a raw data analysis pipeline this would be included. - [paper](https://pubs.acs.org/doi/10.1021/acs.jproteome.0c00563)


2021 - [RawBeans](https://bitbucket.org/incpm/prot-qc/src/master/) is a upgraded program of RawMeat. It's a raw data quaility control tool that help identify insturment issues relating to spray instability, problems with fragmentation or unequal loading. This program can be used on a stand alone PC or included in a pipeline. - [paper](https://pubs.acs.org/doi/10.1021/acs.jproteome.0c00956?goto=supporting-info)


2021 - [mokapot](https://github.com/wfondrie/mokapot) - python -  Semisupervised Learning for Peptide Detection - [paper](https://pubs.acs.org/doi/10.1021/acs.jproteome.0c01010)

2021 - [qcloud2](https://github.com/proteomicsunitcrg/qcloud2-pipeline) - cloud based quality control pipeline, can be integrated with nextflow and openMS - [paper](https://pubs.acs.org/doi/10.1021/acs.jproteome.0c00853)

## 5. Raw Data Analysis

R is a great language to learn for proteomics data analysis. I highly recommend leanring R or Python. From what I've read R is more popular in Academia, while Python is more popular in Industry. Probably because python is great for making piplines and can do more than just statistcal analysis. While R has so many more novel packages for data analysis, that's probably not used a lot in Industry. 

I don't recommend it for large raw files but it's possible to read raw MS files into R and analyze the PSMs. Can be helpful if you want to plot or quantify specific PTMs or peptides for say a PRM experiment.

2018 - [rawDiag](https://github.com/fgcz/rawDiag) - R - Package that can be used in conjustion with rawrr - [paper](https://pubs.acs.org/doi/10.1021/acs.jproteome.8b00173)

2021 - [rawrr](https://github.com/fgcz/rawrr) - R - A great package that can read in raw thermo files! Thats great to me, because I always find it tedious to convert a raw file into a mzML or mzXML file - [paper](https://pubs.acs.org/doi/10.1021/acs.jproteome.0c00866)

2020 - [COSS](https://github.com/compomics/COSS) - java - user-friendly spectral library search tool - [paper](https://pubs.acs.org/doi/10.1021/acs.jproteome.9b00743)

2018 - [pyteomics](https://pypi.org/project/pyteomics/) - python - proteomics framework tools - [paper](https://pubs.acs.org/doi/10.1021/acs.jproteome.8b00576)

## 6. Stastical Analysis

2014 - [MSstats](https://github.com/Vitek-Lab/MSstats) - R - DDA/shotgun, bottom-up, SRM, DIA - [paper](https://academic.oup.com/bioinformatics/article/30/17/2524/2748156?login=false)

2020 - [MSstatsTMT](https://github.com/Vitek-Lab/MSstatsTMT) - R - TMT shotgun proteomics - [paper](https://www.mcponline.org/article/S1535-9476(20)35114-8/fulltext)

2020 - [proteiNorm](https://github.com/ByrumLab/proteiNorm) - R - TMT and unlabeled, has multiple options for normalization and statistical analysis - [paper](https://pubs.acs.org/doi/10.1021/acsomega.0c02564)

2020 - [DEqMS](https://github.com/yafeng/DEqMS) - R - Developed ontop of limma, but takes into account variability in PSMs. Works on both labelled and unlabelled samples - [paper](https://www.mcponline.org/article/S1535-9476(20)34997-5/fulltext)

2018 - [PaDuA](https://github.com/mfitzp/padua) - python - proteomics and phosphoproteomics data analysis - [paper](https://pubs.acs.org/doi/10.1021/acs.jproteome.8b00576)

2021 - [MSstatsPTM[(https://github.com/tsunghengtsai/MSstatsPTM) - labeled and unlabeled PTM data analysis

## 7. R packages for Protein Pathway Enrichment

2019 - [fgsea](https://github.com/ctlab/fgsea) - R - fast gene set enrichment analysis - [paper](https://www.biorxiv.org/content/10.1101/060012v2.full)

2021 - [phosphoRWHN](https://github.com/JoWatson2011/phosphoRWHN) - R - pathway enrichment for phosphoproteomics data - [paper](https://pubs.acs.org/doi/10.1021/acs.jproteome.1c00150)

2021 - [leapR](https://github.com/PNNL-CompBio/leapR) - R - package for multiple pathway analysis - [paper](https://pubs.acs.org/doi/10.1021/acs.jproteome.0c00963)

2019 - [pathfindR](https://github.com/egeulgen/pathfindr) - R - active subnetwork oriented pathway enrichment analyses that uses protein-protein ineteraction networks to enchance the standard pathway analysis method - [paper](https://www.frontiersin.org/articles/10.3389/fgene.2019.00858/full)




2020 - [lipidR](https://www.lipidr.org/) - R - lipidomics data analysis - [paper](https://pubs.acs.org/doi/10.1021/acs.jproteome.0c00082)


## Top down data analysis

2021 - [ClipsMS](https://github.com/loolab2020/ClipsMS) - python -  analysis of terminal and internal fragments in top-down mass spectrometry data - [paper](https://pubs.acs.org/doi/10.1021/acs.jproteome.0c00952)

## Miscellaneous

2021 - [ProteaseGuru](https://github.com/smith-chem-wisc/ProteaseGuru) - C# - tool for In Silico Database Digestion, optimize bottom up experiments - [paper](https://pubs.acs.org/doi/10.1021/acs.jproteome.0c00954)

2020 - [PeCorA](https://github.com/jessegmeyerlab/PeCorA) - R - peptide correlation analsis - [paper](https://pubs.acs.org/doi/10.1021/acs.jproteome.0c00602)

