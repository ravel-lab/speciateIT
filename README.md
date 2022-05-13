![speciateit_logo](https://user-images.githubusercontent.com/17168205/40029457-bf249c04-57b2-11e8-9d2e-85e4ea6f3d0c.png)

### Summary:
speciateIT is an algorithm capable of fast, accurate individual sequence taxonomic classification. Using 7th order Markov Chain models to represent microbial species and a model guide tree, speciateIT requires little computational resources, and can quickly process large sequencing datasets. Currently, speciateIT models have been built from full length 16S rRNA gene sequences for the popular V1-V3, V3-V4, and V4 amplicon regions. speciateIT models correctly classified 98, 97, and 94% of known sequences from the V1-V3, V3-V4, and V4 databases. 

### Compile classification script:  
 #### From top directory:   
  cd src   
  make -f Makefile     
    
  Warning ==> Each model set requires ~30Gb storage space    
    
### Build models:
To build models & estimate error thresholds with V3V4, for example:

  /bin/buildModelTree -l post_merge_reference_files/V3V4_spp_new.lineage -i  reference_sequences/V3V4_trimmed_noEuks_nr_Complete.fa -t post_merge_reference_files/V3V4_spp_new.tx -o V3V4

  /bin/buildMC -t mcDir/spp_paths.txt -k 8 -d V3V4

  /bin/est_error_thlds -d V3V4 -c 0.9 
  
### Classify sequences:   
  /bin/classify -d <model-directory> -i <input-fasta-file> -o <output-directory> --skip-err-thld 
  Classifier will produce output directory if needed.   
  Output file name always "MC_order7_results.txt"
  Output file contains 4 columns: "Sequence ID" "Classification" "posterior probability" "number of Decisions"
    
  Ex. /bin/classify -d V3V4 -i my_V3V4_sample_sequences.fa -o my_V3V4_sample_sequences_classification  
      cd my_V4_sample_sequences_classification  
      ls   
        MC_order7_results.txt  
 
        
    
