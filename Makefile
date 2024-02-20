#############################################################################
#
#               Makefile for building all executables
#
#############################################################################

# Set this to your installation directory.
BIN_DIR = /usr/local/bin
UNZIP_DIR = vSpeciateDB_models

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

BUILDDIR= src/.build

all: src-make_default install unzip

src-make_default: src/Makefile_buildMC src/Makefile_buildModelTree src/Makefile_classify  src/Makefile_pp_embedding src/Makefile_vicut
	cd src && $(MAKE) -f Makefile_buildMC
	cd src && $(MAKE) -f Makefile_buildModelTree
	cd src && $(MAKE) -f Makefile_classify
	cd src && $(MAKE) -f Makefile_pp_embedding
	cd src && $(MAKE) -f Makefile_vicut

src-clean:
	$(DEL_FILE) $(BUILDDIR)/*.o

unzip:
	# Unzip files into the specified directory
	cd $(UNZIP_DIR) && unzip -Xo vSpeciateIT_V1V3.zip
	cd $(UNZIP_DIR) && unzip -Xo vSpeciateIT_V3V4.zip
	cd $(UNZIP_DIR) && unzip -Xo vSpeciateIT_V4V4.zip

install: bin/buildMC bin/buildModelTree bin/classify bin/pp_embedding bin/vicut bin/count_table.py
	$(shell $(CHK_DIR_EXISTS) $(BIN_DIR) || $(MKDIR) $(BIN_DIR))
	$(COPY) bin/buildMC $(BIN_DIR)
	$(COPY) bin/buildModelTree $(BIN_DIR)
	$(COPY) bin/classify $(BIN_DIR)
	$(COPY) bin/pp_embedding $(BIN_DIR)
	$(COPY) bin/vicut $(BIN_DIR)
	$(COPY) bin/count_table.py $(BIN_DIR)

print:
	@echo speciateIT-0.1

make_default: src-make_default FORCE
all: src-make_default FORCE
clean: src-clean FORCE
unzip: unzip FORCE
.PHONY: unzip

FORCE:
