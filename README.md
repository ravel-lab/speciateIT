![speciateit_logo](https://user-images.githubusercontent.com/17168205/40029457-bf249c04-57b2-11e8-9d2e-85e4ea6f3d0c.png)

### Summary:
speciateIT is an algorithm capable of fast, accurate individual sequence taxonomic classification. Using 7th order Markov Chain models to represent microbial species and a model guide tree, speciateIT requires little computational resources, and can quickly process large sequencing datasets. Currently, speciateIT models have been built from full length 16S rRNA gene sequences for the popular V1-V3, V3-V4, and V4 amplicon regions. speciateIT models correctly classified 98, 97, and 94% of known sequences from the V1-V3, V3-V4, and V4 databases.

### Compile classification script:

#### From top directory execute:
  make Makefile

  and then

  make install

Make sure /usr/local/bin is in your PATH

  Warning ==> Each model set requires ~30Gb storage space

NOTE: All Makefile's are configured for macosx. If you are on linux machine,
change CC, CXX and LINK to gcc, g++ and g++, respectively. You may also need to
modify LDFLAGS.

### Build models:
To build models & estimate error thresholds with V3V4, for example:

  buildModelTree -l post_merge_reference_files/V3V4_spp_new.lineage -i  reference_sequences/V3V4_trimmed_noEuks_nr_Complete.fa -t post_merge_reference_files/V3V4_spp_new.tx -o V3V4

  buildMC -v -t V3V4/spp_paths.txt -k 8 -d V3V4

  est_error_thlds -v -d V3V4 -c 0.9

### Classify sequences:

  classify -d <model-directory> -i <input-fasta-file> -o <output-directory> --skip-err-thld

  Classifier will produce output directory if needed.
  Output file name always "MC_order7_results.txt"
  Output file contains 4 columns: "Sequence ID" "Classification" "posterior probability" "number of Decisions"

  Example

  classify -v -i test10k_V3V4.fa -d V3V4_vag_mcDir -o txclass_test10k_V3V4

  classify -v -i test10k_V3V4.fa -d V3V4_vaginal_Mar2021/V3V4_vag_mcDir -o txclass_test10k_V3V4




  The classification results are in txclass_test10k_V3V4/MC_order7_results.txt
