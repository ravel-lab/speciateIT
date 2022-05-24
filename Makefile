#############################################################################
#
#               Makefile for building all executables
#
#############################################################################

# Set this to your installation directory.
BIN_DIR = /usr/local/bin

first: make_default
MAKEFILE        = Makefile
DEL_FILE        = rm -f
CHK_DIR_EXISTS  = test -d
MKDIR           = mkdir -p
COPY            = cp -f
COPY_FILE       = cp -f
COPY_DIR        = cp -f -R
INSTALL_FILE    = $(COPY_FILE)
INSTALL_PROGRAM = $(COPY_FILE)
INSTALL_DIR     = $(COPY_DIR)
DEL_FILE        = rm -f
SYMLINK         = ln -sf
DEL_DIR         = rmdir
MOVE            = mv -f
CHK_DIR_EXISTS  = test -d
MKDIR           = mkdir -p

BUILDDIR= .build

src-make_default: src/Makefile_buildMC src/Makefile_buildModelTree src/Makefile_clError src/Makefile_classify src/Makefile_cut_tree src/Makefile_est_error_thlds  src/Makefile_fp_fn_rates src/Makefile_pp_embedding src/Makefile_pp_ref_sib_wr_ref_models src/Makefile_pp_wr_selected_models src/Makefile_rk_stats src/Makefile_tree_stitcher src/Makefile_vicut
	cd src && $(MAKE) -f Makefile_buildMC
	cd src && $(MAKE) -f Makefile_buildMC
	cd src && $(MAKE) -f Makefile_buildModelTree
	cd src && $(MAKE) -f Makefile_clError
	cd src && $(MAKE) -f Makefile_classify
	cd src && $(MAKE) -f Makefile_cut_tree
	cd src && $(MAKE) -f Makefile_est_error_thlds
	cd src && $(MAKE) -f Makefile_fp_fn_rates
	cd src && $(MAKE) -f Makefile_pp_embedding
	cd src && $(MAKE) -f Makefile_pp_ref_sib_wr_ref_models
	cd src && $(MAKE) -f Makefile_pp_wr_selected_models
	cd src && $(MAKE) -f Makefile_rk_stats
	cd src && $(MAKE) -f Makefile_vicut

src-clean:
	$(DEL_FILE) -f $(BUILDDIR)/*.o

install: bin/buildMC bin/buildModelTree bin/clError bin/classify bin/cut_tree bin/est_error_thlds bin/fp_fn_rates bin/pp_embedding bin/pp_ref_sib_wr_ref_models bin/pp_wr_selected_models bin/rk_stats bin/tree_stitcher bin/vicut
	$(shell $(CHK_DIR_EXISTS) $(BIN_DIR) || $(MKDIR) $(BIN_DIR))
	$(COPY) bin/buildMC $(BIN_DIR)
	$(COPY) bin/buildModelTree $(BIN_DIR)
	$(COPY) bin/clError $(BIN_DIR)
	$(COPY) bin/classify $(BIN_DIR)
	$(COPY) bin/cut_tree $(BIN_DIR)
	$(COPY) bin/est_error_thlds $(BIN_DIR)
	$(COPY) bin/fp_fn_rates $(BIN_DIR)
	$(COPY) bin/pp_embedding $(BIN_DIR)
	$(COPY) bin/pp_ref_sib_wr_ref_models $(BIN_DIR)
	$(COPY) bin/pp_wr_selected_models $(BIN_DIR)
	$(COPY) bin/rk_stats $(BIN_DIR)
	$(COPY) bin/vicut $(BIN_DIR)
	$(COPY) perl/*.pl $(BIN_DIR)

dist:
	cd .. && tar zcvfX speciateIT-0.1.tgz speciateIT/.speciateIT_exclude speciateIT

print:
	echo speciateIT-0.1

make_default: src-make_default FORCE
all: src-make_default FORCE
clean: src-clean FORCE
distclean: dist-clean FORCE

FORCE:
