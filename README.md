![speciateit_logo](https://user-images.githubusercontent.com/17168205/40029457-bf249c04-57b2-11e8-9d2e-85e4ea6f3d0c.png)

### Summary:
speciateIT is an algorithm capable of fast, accurate individual sequence taxonomic classification. Using 7th order Markov Chain models to represent microbial species and a model guide tree, speciateIT requires little computational resources, and can quickly process large sequencing datasets. Currently, speciateIT models have been built from full length 16S rRNA gene sequences for the popular V1-V3, V3-V4, and V4 amplicon regions. speciateIT models correctly classified 98, 97, and 94% of known sequences from the V1-V3, V3-V4, and V4 databases. 

### Compile classify script:  
 #### From top directory:   
  cd src   
  make -f Makefile     
    
### Unarchive models:   
  #### From top directory:   
    cd PECAN_models  
    tar -xvzf V1V3.tar.gz  
    OR  
    tar -xvzf V3V4.tar.gz  
    OR  
    tar -xvzf V4.tar.gz  
     
  Warning ==> Each model set requires ~30Gb storage space    
    
### Classify sequences:   
  ~/bin/classify -d <model-directory> -i <input-fasta-file> -o <output-directory>  
  Classifier will produce output directory if needed.   
  Output file name always "MC_order7_results.txt"  
    
  Ex. ~/bin/classify -d V4 -i my_V4_sample_sequences.fa -o my_V4_sample_sequences_classification  
      cd my_V4_sample_sequences_classification  
      ls   
        MC_order7_results.txt  
        
    
