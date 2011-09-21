#!/bin/sh

# define environmental variables and set program paths

# galaxy installation directory, all data related to the specific galaxy installation will
# be stored  here e.g.
# galaxy_dist/  python_egg_cache/  start.sh@  stop.sh@  tool_conf.xml@  tools/  tools-data/  uct_cbio/  universe_wsgi.ini@
export GALAXY_INSTALL_DIR=/usr/local/share/galaxy

# perl and python libraries
export PERL5LIB=$PERL5LIB:$GALAXY_INSTALL_DIR/uct_cbio/src/tools
export PERL5LIB=$PERL5LIB:$GALAXY_INSTALL_DIR/uct_cbio/src/tools/predict_peptide_seq/lib
export PYTHONPATH=$PYTHONPATH:$GALAXY_INSTALL_DIR/uct_cbio/src/tools/
export PYTHONPATH=$PYTHONPATH:$GALAXY_INSTALL_DIR/uct_cbio/src/tools/converters/lib
export PYTHON_EGG_CACHE=$GALAXY_INSTALL_DIR/python_egg_cache

# reverse complement
export PATH=$PATH:$GALAXY_INSTALL_DIR/uct_cbio/src/tools/dependencies/revcomp

# blast
export BLASTDB=/usr/local/share/paeano/data/blast_db
export BLASTMAT=/usr/local/share/paeano/software/blast/data
export PATH=$PATH:/usr/local/share/paeano/software/blast/bin

# basecall
export PATH=$PATH:/usr/local/share/paeano/software/phred
export PHRED_PARAMETER_FILE=/usr/local/share/paeano/software/phred/phredpar.dat

# cluster
export PATH=$PATH:/usr/local/share/paeano/software/wcd/src
export PATH=$PATH:$GALAXY_INSTALL_DIR/uct_cbio/src/tools/dependencies/clobb2

# assembly
export PATH=$PATH:/usr/local/share/paeano/software/cap3
export PATH=$PATH:/usr/local/share/paeano/software/phrap
export PATH=$PATH:$GALAXY_INSTALL_DIR/uct_cbio/src/tools/dependencies/contigimage

# protein prediction
export PATH=$PATH:/usr/local/share/paeano/software/estscan
export PATH=$PATH:/usr/local/share/paeano/software/decoder/src
export ESTSCANDIR=/usr/local/share/paeano/software/ESTScan2-2.1
export ESTSCAN_MATRICES_DIR=/usr/local/share/paeano/software/estscan/matrices
export PATH=$PATH:/usr/local/share/paeano/software/EMBOSS-6.0.1/emboss
export RRNA_DB=$GALAXY_INSTALL_DIR/uct_cbio/src/tools/predict_peptide_seq/rRNA/rRNAeukAll_nuc.fsa
export MITO_DB=$GALAXY_INSTALL_DIR/uct_cbio/src/tools/predict_peptide_seq/mito/mitoMetazoan90_pro.fsa
export GEN_CODE_FILE=$GALAXY_INSTALL_DIR/uct_cbio/src/tools/predict_peptide_seq/genetic_code/gc.prt

# blast2go
export PATH=$PATH:$GALAXY_INSTALL_DIR/uct_cbio/src/tools/dependencies/interproscan
export BLAST2GO_JAR_FILE=/usr/local/share/paeano/software/blast2go/b2g4pipe/blast2go.jar
export BLAST2GO_PROPERTIES_FILE=/usr/local/share/paeano/software/blast2go/b2g4pipe/blast2go.properties

# metagenemark
export PATH=$PATH:/usr/local/share/paeano/software/metagenemark
export METAGENEMARK_MODEL_FILE=/usr/local/share/paeano/software/metagenemark/MetaGeneMark_v1.mod

# interproscan
export PATH=$PATH:/usr/local/share/paeano/software/interproscan/iprscan/bin
export IPRSCAN_OUTPUT_DIR=/usr/local/share/paeano/software/interproscan/iprscan/tmp
export IPRSCAN_IMAGES_DIR=/usr/local/share/paeano/software/interproscan/iprscan/images

# ngs software
export PATH=$PATH:/usr/local/share/bwa
export PATH=$PATH:/usr/local/share/samtools
export PATH=$PATH:/usr/local/share/bowtie

# change to galaxy dist directory (this is where the run.sh script is)
export GALAXY_DIST_DIR=$GALAXY_INSTALL_DIR/galaxy_dist
export TEMP=$GALAXY_DIST_DIR/database/tmp
cd $GALAXY_DIST_DIR

# start daemon
sudo -u galaxy TEMP=$TEMP PERL5LIB=$PERL5LIB  PYTHONPATH=$PYTHONPATH PHRED_PARAMETER_FILE=$PHRED_PARAMETER_FILE BLASTDB=$BLASTDB BLASTMAT=$BLASTMAT PYTHON_EGG_CACHE=$PYTHON_EGG_CACHE PATH=$PATH RRNA_DB=$RRNA_DB MITO_DB=$MITO_DB GEN_CODE_FILE=$GEN_CODE_FILE ESTSCAN_MATRICES_DIR=$ESTSCAN_MATRICES_DIR BLAST2GO_JAR_FILE=$BLAST2GO_JAR_FILE BLAST2GO_PROPERTIES_FILE=$BLAST2GO_PROPERTIES_FILE METAGENEMARK_MODEL_FILE=$METAGENEMARK_MODEL_FILE IPRSCAN_OUTPUT_DIR=$IPRSCAN_OUTPUT_DIR IPRSCAN_IMAGES_DIR=$IPRSCAN_IMAGES_DIR sh multiprocess.sh --daemon
