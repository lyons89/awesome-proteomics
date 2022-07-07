# awesome-proteomics
An awesome list of proteomics tools and resources

I've been analyzing large scale proteomics data sets for over 5 years now. I've recently stumbled upon these awesome lists and wanted to make an in-depth one for proteomics.


## Proteomics 

### Table of Contents 
1. Learning Resources
2. Databases
3. Raw data search software/algorithms
4. Pipeline tools
5. R Packages for Raw Data Analysis
6. R Packages for Stastical Analysis
7. R packages for Protein Pathway Enrichment
8. R Packages for Kinase Motif/Activity Analysis



## 1. Leaning Resources

[Proteomics blog](https://proteomicsnews.blogspot.com/) - hands down the best protomeics blog I've seen. Ben Orsburn is very knowledgable with a great sense of humor.

[Phil Wilmarth](https://github.com/pwilmart/Start_Here) has alot of goodlearning resources. He was python scripts for raw data analysis as well as blog posts detailing TMT data analysis techniques. 

[Review article](https://www.nature.com/articles/nrm1468) - great proteomics review article

[Tutorial videos](https://www.youtube.com/channel/UC0v4sjdXLMa-OWR7IYeoFoA/videos) - videos from NCQBCS, a project lead by the Coon lab. Contains lots of information regarding Experimental design, ionization, quantitative proteomics, analysis, post-translational modifications and more.

[ASMS video mass spec channel](https://vimeo.com/channels/asms) contains a lot of videos from leading researchers in the protoemics field. 

[Videos from Nikolai](https://www.youtube.com/c/NikolaiSlavovResearch/videos) most of his videos focus on single cell proteomics, and DIA

[Videos from Matthew Padula](https://www.youtube.com/c/MatthewPadula/videos) lots of great videos on the basis of mass spectrometry and proteomics.
## 2. Databases

Commonly used protoemics data repoitories. You can find raw data here, or publish your down for article publication requirements. 

[ProteomeXchange](http://www.proteomexchange.org/) - proteomeXchange is a global repository that contains links to all major databases, including MassIVE, Pride, iProX and more. Probably the best place to start.

[massIVE](https://massive.ucsd.edu/ProteoSAFe/static/massive.jsp)

[Pride database](https://www.ebi.ac.uk/pride/archive/)


This massive master list of databses from Pastel BioScience

http://www.pastelbioscience.co.uk/resources/databases.html - this list contains so much useful information that I'm sure the rest of my awesome list will be redundant

[Uniprot](https://www.uniprot.org/)

## 3. Raw data search software/algorithms

Maybe you're doing a collaboration with a mass spec lab or working with the core. They sent you back raw data files and you have no idea what to do with them. Alternatively, you're just starting out in proteomics and want to find a better software to use to search your raw MS files.


[Fragpipe](https://github.com/Nesvilab/FragPipe) currently my favorite raw file search software. It's much faster than maxquant and in my opinion has sleeker GUI. In a recent update it also allows searches to be set up on a linux server for even faster results. The software is very modular, it consists of [MSfragger](https://msfragger.nesvilab.org/) the database search algorithm, [Philosopher](https://philosopher.nesvilab.org/) that analyzes the database results, as well as others for PTM and TMT integration. 

[MaxQuant](https://www.maxquant.org/) is probably the most used and well known DDA software. Developed by Jurgon Cox, this completely free software is user friendly and is always being updated with new and original features. There is even a [youtube](https://www.youtube.com/c/MaxQuantChannel) that has tons of videos on how to use the software. 



