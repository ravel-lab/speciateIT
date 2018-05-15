# PECAN
Per sEquenCe tAxonomic assigNer

Compile classify script:
  From top directory: 
    cd src
    make -f Makefile
    
Unarchive models: 
  From top directory: 
    cd PECAN_models
    tar -xvzf V1V3.tar.gz
    OR
    tar -xvzf V3V4.tar.gz
    OR
    tar -xvzf V4.tar.gz
   
  Warning ==> Each model set requires ~30Gb storage space
  
Classify sequences: 
  ~/bin/classify -d <model-directory> -i <input-fasta-file> -o <output-directory>
  Classifier will produce output directory if needed. 
  Output file name always "MC_order7_results.txt"
  
  Ex. ~/bin/classify -d V4 -i my_V4_sample_sequences.fa -o my_V4_sample_sequences_classification
      cd my_V4_sample_sequences_classification
      ls 
        MC_order7_results.txt
        
    
