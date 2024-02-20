![speciateit_logo](https://user-images.githubusercontent.com/17168205/40029457-bf249c04-57b2-11e8-9d2e-85e4ea6f3d0c.png)

### Summary:
SpeciateIT is an algorithm capable of fast, accurate individual sequence taxonomic classification. vSpeciateDB are models built from custom sets of reference sequences for classifying vaginal microbiota. Using a model guide tree and 7th order Markov Chain models to represent bacterial species trained on taxonomy-adjusted amplicon specific regions sequences, speciateIT requires little computational resources, and can quickly process large sequence datasets. SpeciateIT models for the vaginal microbiota include training sets for 16S rRNA gene V1-V3, V3-V4, and V4 regions sequences. "Cat Maps" are provided for each region to indicate which species are indistinguishable at the targeted variable regions.

### vSpeciateDB:
Holm, Johanna (2024). speciateIT: vSpeciateDB Models. figshare. Dataset. https://doi.org/10.6084/m9.figshare.25254229.v1

### To install:
1. Clone repository.
2. Download vSpeciateDB and place into "vSpeciateDB_models". 
3. From parent directory
   
   make all
   
   This will compile all classification scripts, add them to /usr/local/bin, and unarchive vSpeciateDB.

   NOTE: May require "sudo" to access /usr/local/bin
   sudo make all

### To use: 

classify -d < vSpeciateDB dir > -i < fasta file > -o < outDir >

  Example: classify -d vSpeciateDB_models/vSpeciateIT_V3V4 -i MyProject_ASV.fasta -o MyProject 

   Classifier will produce output directory if needed.

   Output file name always "MC_order7_results.txt"

   Output file contains 4 columns: "Sequence ID" "Classification" "posterior probability" "number of Decisions"

To force species-level annotations (i.e. ignore error thresholds): 
classify -d < vSpeciateDB dir > -i < fasta file > -o < outDir > --skip-err-thld

NOTE: All Makefile's are configured for macosx. If you are on linux machine,
change CC, CXX and LINK to gcc, g++ and g++, respectively. You may also need to
modify LDFLAGS.
