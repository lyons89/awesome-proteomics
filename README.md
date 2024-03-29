
[![Awesome](https://awesome.re/badge.svg)](https://awesome.re)
[![Awesome](https://awesome.re/badge-flat.svg)](https://awesome.re)
[![Awesome](https://awesome.re/badge-flat2.svg)](https://awesome.re)

# awesome-proteomics
An awesome list of proteomics tools and resources.

I've been analyzing large scale proteomics data sets for over 5 years now. I've recently stumbled upon these awesome lists and wanted to make an one for proteomics.


## Proteomics 

### Table of Contents 
1A. Learning Resources - Proteomics
1B. Learning Resources - Programming
2. Databases
3. Raw data search software/algorithms
4. Assorted Pipeline Tools
5. Raw Data Analysis
6. Stastical Analysis
7. Protein Pathway Enrichment
8. Kinase Motif/Activity Analysis
9. Top down data analysis
10. Multi-Omics data analysis




## 1A. Leaning Resources - Proteomics

[Ben Orsburn](https://proteomicsnews.blogspot.com/) has hands down the best protomeics blog I've ever seen. Ben is very knowledgable with a great sense of humor.

[Phil Wilmarth](https://github.com/pwilmart/Start_Here) has alot of goodlearning resources. He was python scripts for raw data analysis as well as blog posts detailing TMT data analysis techniques. 

[biostars](https://www.biostars.org/) is like stackoverflow but only for bioinformatics, has forums for questions, job postings and tutorials.

[Review article](https://www.nature.com/articles/nrm1468) about protoemics.

[Tutorial videos](https://www.youtube.com/channel/UC0v4sjdXLMa-OWR7IYeoFoA/videos) from NCQBCS, a project led by the Coon lab. Contains lots of information regarding experimental design, ionization, quantitative proteomics, analysis, post-translational modifications and more.

[ASMS video mass spec channel](https://vimeo.com/channels/asms) contains a lot of videos from leading researchers in the protoemics field. 

[Videos from Nikolai](https://www.youtube.com/c/NikolaiSlavovResearch/videos) most of his videos focus on single cell proteomics, and DIA.

[Videos from Matthew Padula](https://www.youtube.com/c/MatthewPadula/videos) lots of great videos on the basis of mass spectrometry and proteomics.

[MayInstuite](https://github.com/MayInstitute) Computational proteomics short courses organized by Olga Vitek. They also have ALOT of videos on [youtube](https://www.youtube.com/@MayInstituteNEU).

## 1B. Learning Resources - Programming

[R books](https://bookdown.org/) large abundance of ebooks for learning R, from basic R, to advance R, Shiny and more!!

[Are-we-learning-yet](https://github.com/anowell/are-we-learning-yet?tab=readme-ov-file) a resource to learn machine learning in rust. I've noticed rust is starting to become a more popular language, not only in the wild but also in proteomics (see Sage below).

[conda/bioconda](https://github.com/MonashBioinformaticsPlatform/bioconda-tutorial) anaconda is a popular bioinformatics tool used mostly for python, and for some R, programming. Useful for creating reproducible enrivonments.

[Intro Math](https://github.com/erikaduan/introductory_maths) if you're like me and it's been a few years since you've had to use math check out this repo to brush up on it and learn some R, Python and Julia.

[python data science tips](https://github.com/erikaduan/python_data_science_tips/tree/master/notebooks) another resource to learn python/pandas with.

[Mass Spec Coding Club](https://github.com/michaelmarty/MassSpecCodingClub) great resource to learn python and then apply that knowledge to mass spectrometry.


## 2. Databases

[ProteomeXchange](http://www.proteomexchange.org/) is a global repository for raw MS data that contains links to all major databases, including MassIVE, Pride, iProX and more. Probably the best place to start.

[Pastel BioScience](http://www.pastelbioscience.co.uk/resources/databases.html) has a database that contains staggering amounts information that I'm sure the rest of my awesome list will be redundant.

[Uniprot](https://www.uniprot.org/) - Has all the information you will ever need to know for individual proteins and the go to for protein FASTA databases.

[Biogrid](https://thebiogrid.org/) - protein-protein interaction database

[KEGG](https://www.genome.jp/kegg/) - biological pathway database

[Reactome](https://reactome.org/) - nicer looking biological pathway database

2021 - [CPTAC](https://github.com/PayneLab/cptac) - python/R - API interaface to publically available cancer datasets - [paper](https://pubs.acs.org/doi/10.1021/acs.jproteome.0c00919)

2021 - [ppx](https://github.com/wfondrie/ppx) - python - Python interface to proteomics data repositories - [paper](https://pubs.acs.org/doi/10.1021/acs.jproteome.1c00454)

## 3. Raw data search software/algorithms

2017 - [Fragpipe](https://github.com/Nesvilab/FragPipe) - Java - It's a very fast search engine with a nice GUI. The software is modular, it consists of [MSfragger](https://msfragger.nesvilab.org/) the database search algorithm, [Philosopher](https://philosopher.nesvilab.org/) that analyzes the database results, as well as others for PTM and TMT integration. 

2008 - [MaxQuant](https://www.maxquant.org/) is probably the most used and well known DDA software. Developed by Jurgon Cox, this completely free software is user friendly and is always being updated with new and original features. There is even a [youtube](https://www.youtube.com/c/MaxQuantChannel) that has tons of videos on how to use the software. - [paper](https://www.nature.com/articles/nbt.1511)

2015 - [Peptide-shaker](https://compomics.github.io/projects/peptide-shaker) is like the swiss army knife of search tools. You can search data with multiple search engines inclduing, comet, tide, andromeda, mascot, X!Tandem and more that I've never heard of. [paper](https://www.nature.com/articles/nbt.3109)

2010 - [skyline](https://skyline.ms/project/home/software/skyline/begin.view) - software for targeted proteomics - [paper](https://pubmed.ncbi.nlm.nih.gov/20147306/)

2012 - [Comet](https://uwpr.github.io/Comet/) - C++ - Free and open-source search engine, lately it's had several  - [paper](https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/abs/10.1002/pmic.201200439)

2020 - [DIA-NN](https://github.com/vdemichev/DiaNN) - C/C++ - free and open source search tool for DIA data that uses neural networks, works using either a library or a FASTA database. - [paper](https://www.nature.com/articles/s41592-019-0638-x)

2023 - [Sage](https://github.com/lazear/sage) - Rust - most likely the current fasest search engine, it's completely terminal based but if you learn to use it, it will be worth it - [paper](https://pubs.acs.org/doi/10.1021/acs.jproteome.3c00486?ref=PDF)

## 4. Assorted pipeline Tools

2019 - [MaxQuant Live](https://maxquant.org/mqlive/) (not sure where to put this) for real time monitoring of MS data and acquistion. 

2009 - [PAW_pipeline](https://github.com/pwilmart/PAW_pipeline) - python - a pretty much stock python raw file protoemics pipeline tool. It includes functions, to convert files, run comet, produce histograms. Can also do TMT - [paper](https://pubmed.ncbi.nlm.nih.gov/20157357/)

2015 - [Ursgal](https://github.com/ursgal/ursgal) - python - combines multiple search engine algorithms, postprocessing algorithms, and stastis on the output from multiple search engines - [paper1](https://pubs.acs.org/doi/10.1021/acs.jproteome.5b00860) [paper2](https://pubs.acs.org/doi/10.1021/acs.jproteome.0c00799)

2019 - [DIAlignR](https://github.com/Roestlab/DIAlignR) - R - DIA retention time alignment of targetd MS data, including DIA and SWATH-MS - [paper](https://www.mcponline.org/article/S1535-9476(20)31843-0/fulltext)

2021 - [Monocle](https://github.com/gygilab/Monocle) - C# - for monoisotopic peak and accurate precursor m/z detection in shotgun proteomics experiments. - [paper](https://pubs.acs.org/doi/10.1021/acs.jproteome.0c00563)

2021 - [RawBeans](https://bitbucket.org/incpm/prot-qc/src/master/) is a upgraded program of RawMeat. It's a raw data quaility control tool that help identify insturment issues relating to spray instability, problems with fragmentation or unequal loading. This program can be used on a stand alone PC or included in a pipeline. - [paper](https://pubs.acs.org/doi/10.1021/acs.jproteome.0c00956?goto=supporting-info)

2021 - [mokapot](https://github.com/wfondrie/mokapot) - python -  Semisupervised Learning for Peptide Detection - [paper](https://pubs.acs.org/doi/10.1021/acs.jproteome.0c01010)

2021 - [qcloud2](https://github.com/proteomicsunitcrg/qcloud2-pipeline) - cloud based quality control pipeline, can be integrated with nextflow and openMS - [paper](https://pubs.acs.org/doi/10.1021/acs.jproteome.0c00853)

2021 - [DIAproteomics](https://www.openms.de/comp/diaproteomics/) - python - a module that can be added to a openMS workflow for the analysis of DIA data - [paper](https://pubs.acs.org/doi/10.1021/acs.jproteome.1c00123)



## 5. Raw Data Analysis

2012/2020 - [MSnbase](https://github.com/lgatto/MSnbase) - R - provides MS data structures, allows you to process, quantify, visualize raw data - [paper](https://academic.oup.com/bioinformatics/article/28/2/288/199094?login=false)

2015 - [MaRaCluster](https://github.com/statisticalbiotechnology/maracluster) - C++ - clustering technique to identify fragment spectra stemming from the same peptide species - [paper](https://pubs.acs.org/doi/10.1021/acs.jproteome.5b00749)

2015 - [pyproteome](https://github.com/white-lab/pyproteome) - python - analyzes proteomics data, can filter, normalize, perform motif and pathway enrichment. Currently only supports ProteomeDiscoverer .msf search files - [paper](https://pyproteome.readthedocs.io/en/latest/)

2018 - [pyteomics](https://pypi.org/project/pyteomics/) - python - proteomics framework tools - [paper](https://pubs.acs.org/doi/10.1021/acs.jproteome.8b00576)

2018 - [RawTools](https://github.com/kevinkovalchik/RawTools) - C# - quality control checking of raw files, can assist in method development and insturment quality control - [paper](https://pubs.acs.org/doi/full/10.1021/acs.jproteome.8b00721)

2018 - [MSstatsQC](https://www.bioconductor.org/packages/release/bioc/html/MSstatsQC.html) - R - provides methods for multiple peptide monitoring using raw MS files, works for DDA and DIA data - [paper](https://pubs.acs.org/doi/full/10.1021/acs.jproteome.8b00732)

2018 - [rawDiag](https://github.com/fgcz/rawDiag) - R - Package that can be used in conjustion with rawrr - [paper](https://pubs.acs.org/doi/10.1021/acs.jproteome.8b00173)

2020 - [COSS](https://github.com/compomics/COSS) - java - user-friendly spectral library search tool - [paper](https://pubs.acs.org/doi/10.1021/acs.jproteome.9b00743)

2021 - [rawrr](https://github.com/fgcz/rawrr) - R - A great package that can read in raw thermo files! Thats great to me, because I always find it tedious to convert a raw file into a mzML or mzXML file - [paper](https://pubs.acs.org/doi/10.1021/acs.jproteome.0c00866)

2021 - [PSpecteR](https://github.com/EMSL-Computing/PSpecteR) - R - User Friendly and Interactive for Visualizing the quality of Top-Down and Bottom-Up Proteomics - [paper](https://pubs.acs.org/doi/10.1021/acs.jproteome.0c00857)

2022 - [RforMassSpectrometry](https://www.rformassspectrometry.org/) - R - a massive project that contains multiple helpful packages including RforMassSpectrometry, MsExperiment, Spectra, QFeatures, PSMatch, Chromatograms, MsCoreUtils, and MetaboCoreUtils. 

2023 - [mpwR](https://github.com/OKdll/mpwR) - R - package that allows you to directly compare the output of raw search engines such as MQ, DIANN, spectronaut and I think PD. It's also helpful if you're testing out different settings within your search engine and you want to quickly see how each performs. - [paper](https://pubmed.ncbi.nlm.nih.gov/37267150/)

## 6. Stastical Analysis

2014 - [MSstats](https://github.com/Vitek-Lab/MSstats) - R - DDA/shotgun, bottom-up, SRM, DIA - [paper](https://academic.oup.com/bioinformatics/article/30/17/2524/2748156?login=false)

2018 - [PaDuA](https://github.com/mfitzp/padua) - python - proteomics and phosphoproteomics data analysis - [paper](https://pubs.acs.org/doi/10.1021/acs.jproteome.8b00576)

2020 - [MSstatsTMT](https://github.com/Vitek-Lab/MSstatsTMT) - R - TMT shotgun proteomics - [paper](https://www.mcponline.org/article/S1535-9476(20)35114-8/fulltext)

2020 - [proteiNorm](https://github.com/ByrumLab/proteiNorm) - R - TMT and unlabeled, has multiple options for normalization and statistical analysis - [paper](https://pubs.acs.org/doi/10.1021/acsomega.0c02564)

2020 - [DEqMS](https://github.com/yafeng/DEqMS) - R - Developed ontop of limma, but takes into account variability in PSMs. Works on both labelled and unlabelled samples - [paper](https://www.mcponline.org/article/S1535-9476(20)34997-5/fulltext)

2021 - [MSstatsPTM](https://github.com/tsunghengtsai/MSstatsPTM) - labeled and unlabeled PTM data analysis - [paper](https://www.mcponline.org/article/S1535-9476(22)00285-7/fulltext)

[PermFDP](https://github.com/steven-shuken/permFDP) - R - Package to perform multiple hypothesis correction using permutation based FDP. One of the better performing methods for multiple test corrections. - [paper](https://pubs.acs.org/doi/full/10.1021/acs.analchem.2c03719?casa_token=4CgZMMnAmjgAAAAA%3A-8SyKwz2Hs3L-yRXXQDkq45ZBPW8nexcpqaYzDB5Ok-Kgp0C_W9KPscLE-zUfN2nZUv8uiNYZZcCxy-7) the paper isn't on the tool, it's just a paper that uses it and compares it to other methods. 



## 7. Protein Pathway Enrichment

2019 - [fgsea](https://github.com/ctlab/fgsea) - R - fast gene set enrichment analysis - [paper](https://www.biorxiv.org/content/10.1101/060012v2.full)

2019 - [pathfindR](https://github.com/egeulgen/pathfindr) - R - active subnetwork oriented pathway enrichment analyses that uses protein-protein ineteraction networks to enchance the standard pathway analysis method - [paper](https://www.frontiersin.org/articles/10.3389/fgene.2019.00858/full)

2020 - [lipidR](https://www.lipidr.org/) - R - lipidomics data analysis - [paper](https://pubs.acs.org/doi/10.1021/acs.jproteome.0c00082)

2021 - [phosphoRWHN](https://github.com/JoWatson2011/phosphoRWHN) - R - pathway enrichment for phosphoproteomics data - [paper](https://pubs.acs.org/doi/10.1021/acs.jproteome.1c00150)

2021 - [leapR](https://github.com/PNNL-CompBio/leapR) - R - package for multiple pathway analysis - [paper](https://pubs.acs.org/doi/10.1021/acs.jproteome.0c00963)


## 8. Kinase Motif/Activity Analysis

2017 - [KSEAapp](https://github.com/casecpb/KSEA/) - R - Kinase substrate enrichment analysis. I would recommend using with a freshly downloaded kinase-substarte database from phosphositeplus - [paper](https://pubmed.ncbi.nlm.nih.gov/28655153/)

2015 - [rnotifx](https://github.com/omarwagih/rmotifx) - R - motif enrichment analyssis of PTMs on proteins, probably mostly used for phosphorylation - [paper](https://pubmed.ncbi.nlm.nih.gov/26572964/)


## 9. Top down data analysis

2021 - [ClipsMS](https://github.com/loolab2020/ClipsMS) - python -  analysis of terminal and internal fragments in top-down mass spectrometry data - [paper](https://pubs.acs.org/doi/10.1021/acs.jproteome.0c00952)



## 10. Multi-Omics data analysis

2015 - [moCluster](https://www.bioconductor.org/packages/release/bioc/html/mogsa.html) - R - Integration of multiple omics datasets to identify patterns - [paper](https://pubs.acs.org/doi/10.1021/acs.jproteome.5b00824)

2019 - [MOGSA](https://github.com/kerschke/mogsa) - R - Multiple omics data integrative clustering and gene set analysis - [paper](https://www.mcponline.org/article/S1535-9476(20)32768-7/fulltext)



## Miscellaneous

[cytoscape](https://cytoscape.org/) - visualizing protein-protein interaction netweorks

2012 - [ProteoWizard](https://proteowizard.sourceforge.io/publications.html) - Great software for converting one MS file type to another. I mostly use it ot convert thermo .raw files to mzML - [paper](https://www.nature.com/articles/nbt.2377)

2019 - [IPSC](https://github.com/coongroup/IPSA) - Interactive Peptide Spectrum Annotator, web based utility for shotgun mass spectrum annotation - [paper](https://www.mcponline.org/article/S1535-9476(20)32771-7/fulltext)

2020 - [PeCorA](https://github.com/jessegmeyerlab/PeCorA) - R - peptide correlation analysis - [paper](https://pubs.acs.org/doi/10.1021/acs.jproteome.0c00602)

2021 - [ProteaseGuru](https://github.com/smith-chem-wisc/ProteaseGuru) - C# - tool for In Silico Database Digestion, optimize bottom up experiments - [paper](https://pubs.acs.org/doi/10.1021/acs.jproteome.0c00954)

2021 - [DeepLC](https://github.com/compomics/DeepLC) - python - predicts retention times for peptides that have unseen modifications - [paper](https://www.nature.com/articles/s41592-021-01301-5)
