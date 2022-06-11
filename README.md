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

% cd ~/devel/speciateIT/data
% buildModelTree -l ../post_merge_reference_files/V3V4_spp_new.lineage -i  ../reference_sequences/V3V4_trimmed_noEuks_nr_Complete.fa -t ../post_merge_reference_files/V3V4_spp_new.tx -o sIT_models_V3V4

% buildMC -v -d sIT_models_V3V4

% est_error_thlds -v -d sIT_models_V3V4 -c 0.9

% classify --skip-err-thld -v -i ~/projects/ASV_files/data/CONTRA_ASVs_nr.fa -d sIT_models_V3V4 -o CONTRA_ASVs_sIT


% classify -v -i ~/projects/ASV_files/data/CONTRA_ASVs_nr.fa -d sIT_models_V3V4 -o CONTRA_ASVs_sIT

% classify -v -i ~/projects/ASV_files/data/DOUCHING_ASVs_nr.fa -d sIT_models_V3V4 -o DOUCHING_ASVs_sIT

% classify -v -i ~/projects/ASV_files/data/GALE_ASVs_nr.fa -d sIT_models_V3V4 -o GALE_ASVs_sIT

% classify -v -i ~/projects/ASV_files/data/HMP_ASVs_nr.fa -d sIT_models_V3V4 -o HMP_ASVs_sIT

% classify -v -i ~/projects/ASV_files/data/LSVF_ASVs_nr.fa -d sIT_models_V3V4 -o LSVF_ASVs_sIT

% classify -v -i ~/projects/ASV_files/data/V400_ASVs_nr.fa -d sIT_models_V3V4 -o V400_ASVs_sIT


### Building vaginal sIT models

% cd ~/devel/speciateIT/data

% buildModelTree -l V3V4_db/V3V4_trimmed_noEuks_nr_Complete_vaginal.lineage -i  V3V4_db/V3V4_trimmed_noEuks_nr_Complete_vaginal.fa -t V3V4_db/V3V4_trimmed_noEuks_nr_Complete_vaginal.tx -o vag_sIT_models_V3V4

% buildMC -v -d vag_sIT_models_V3V4

% cd ~/devel/speciateIT/src_tests

% readTable_test /Users/pgajer/devel/speciateIT/data/vag_sIT_models_V3V4/rLpps/Abiotrophia_defectiva.csv


### Gnerating lpp's of model random samples


### Classify sequences:

  classify -d <model-directory> -i <input-fasta-file> -o <output-directory> --skip-err-thld

  Classifier will produce output directory if needed.
  Output file name always "MC_order7_results.txt"
  Output file contains 4 columns: "Sequence ID" "Classification" "posterior probability" "number of Decisions"

  Example

  classify -v -i test10k_V3V4.fa -d V3V4_vag_mcDir -o txclass_test10k_V3V4

  classify -v -i test10k_V3V4.fa -d V3V4_vaginal_Mar2021/V3V4_vag_mcDir -o txclass_test10k_V3V4


### renaming leaves of a tree

/Users/pgajer/opt/anaconda3/bin/nw_rename

### Debugging

cd ~/devel/speciateIT/data/V3V4_vaginal_Mar2021

lldb -- sp_model_seq_lpps -v -d V3V4_vag_mcDir -s 1000


### ToDo's

* Encapsulate the timing routine in el_time()

* Remove from inPar_t parameters not used in the given program

* Check input parameters within inPar_t::check_parameters()

* Organize the code and simplify
