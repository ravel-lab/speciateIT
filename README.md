![speciateit_logo](https://user-images.githubusercontent.com/17168205/40029457-bf249c04-57b2-11e8-9d2e-85e4ea6f3d0c.png)

### SUMmmary:
SpeciateIT is an algorithm capable of fast, accurate individual sequence taxonomic classification. vSpeciateDB are models built from custom sets of reference sequences for classifying vaginal microbiota. Using a model guide tree and 7th order Markov Chain models to represent bacterial species trained on taxonomy-adjusted amplicon specific regions sequences, speciateIT requires little computational resources, and can quickly process large sequence datasets. SpeciateIT models for the vaginal microbiota include training sets for 16S rRNA gene V1-V3, V3-V4, and V4 regions sequences. "Cat Maps" are provided for each region to indicate which species are indistinguishable at the targeted variable regions.

### vSpeciateDB:
Holm, Johanna (2024). speciateIT: vSpeciateDB Models. figshare. Dataset. https://doi.org/10.6084/m9.figshare.25254229

### Requirements: 
- **RAM**: Minimum 1 GB for classification. 
- **Storage**: Minimum 7.5 GB (each set of models is 2.5 GB, uncompressed). 

### To install:
This software runs on **Linux** or **MacOSX**. See specific executables within bin. 

1. Clone repository.
2. Download vSpeciateDB (https://doi.org/10.6084/m9.figshare.25254229) and place into "vSpeciateDB_models".
3. Unzip vSpeciateDB directories.  

   cd /path/to/speciateIT/vSpeciateDB_models

   unzip vSpeciateIT_V1V3  
   unzip vSpeciateIT_V3V4  
   unzip vSpeciateIT_V4V4

4. Add classify to PATH or copy correct executable to your bin:

   export PATH="/path/to/speciateIT/bin/OS/:$PATH"

### To use: 

1. Classification:  
   classify -d < vSpeciateDB dir > -i < fasta file > -o < outDir >  
&nbsp;&nbsp;&nbsp;&nbsp;Example: classify -d vSpeciateDB_models/vSpeciateIT_V3V4 -i MyProject_ASV.fasta -o MyProject   
&nbsp;&nbsp;&nbsp;&nbsp;Output file name is always "MC_order7_results.txt". Format: "Sequence ID" "Classification" "posterior probability" "number of Decisions".  

&nbsp;&nbsp;&nbsp;&nbsp;To force species-level annotations (i.e. ignore error thresholds):   
&nbsp;&nbsp;&nbsp;&nbsp;classify -d < vSpeciateDB dir > -i < fasta file > -o < outDir > --skip-err-thld

2. Make sample x taxon count table using classifications:   
&nbsp;&nbsp;&nbsp;&nbsp;count_table.py -s < MC_order7_results.txt > -c < sample x ASV count table >

### To test:

classify -d vSpeciateDB_models/vSpeciateIT_V3V4 -i test.fasta -o test

more test/MC_order7_results.txt
<pre>
Seq	Classification		pp		nDecisions
ASV1	Lactobacillus_iners	0.970448	50  
ASV2	Lactobacillus_crispatus	0.974671	50  
ASV3	Lactobacillus_mulieris	0.976813	50  
ASV4	Ca_Lachnocurva_vaginae	0.969937	187  
ASV5	Gardnerella_vaginalis	0.974747	27  
ASV6	Lactobacillus_crispatus	0.965138	50  
ASV7	Lactobacillus_iners	0.973368	50  
ASV8	Fannyhessea_vaginae	0.984829	24  
ASV9	Leptotrichia_shahii	0.916945	26  
ASV10	Megasphaera_lornae	0.973400	27 

TO COMPILE FROM SOURCE: From src directories (either macosx or linux) of cloned repository.

   make all

   This will compile all classification scripts.

   NOTE1: May require "sudo" to access /usr/local/bin
   sudo make all
   NOTE2: All Makefile's are configured for either linux or macosx.

</pre>
