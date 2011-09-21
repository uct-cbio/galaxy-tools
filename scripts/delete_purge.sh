#!/bin/sh

# http://bitbucket.org/galaxy/galaxy-central/wiki/Config/PurgeHistoriesAndDatasets

# at the moment the delete_datasets do not work (step 5). if a dataset needs to be removed the folder needs to be removed first.

# define environmental variables and set program paths

# galaxy installation directory, all data related to the specific galaxy installation will
# be stored  here e.g.
# galaxy_dist/  python_egg_cache/  start.sh@  stop.sh@  tool_conf.xml@  tools/  tools-data/  uct_cbio/  universe_wsgi.ini@
export GALAXY_INSTALL_DIR=/usr/local/share/galaxy

# perl and python libraries
export PYTHON_EGG_CACHE=$GALAXY_INSTALL_DIR/python_egg_cache

# change to galaxy dist directory (this is where the run.sh script is)
export GALAXY_DIST_DIR=$GALAXY_INSTALL_DIR/galaxy_dist
cd $GALAXY_DIST_DIR

# number of days to use as cutoff
export DAYS=0

# flags
# -i just give info without changes or removing
# -r remove
# -f force reload
# e.g FLAGS="-r -f"
export FLAGS="-r"

# delete and purge
# 1 - delete_userless_histories.sh 
sudo -u galaxy PYTHON_EGG_CACHE=$PYTHON_EGG_CACHE python ./scripts/cleanup_datasets/cleanup_datasets.py universe_wsgi.ini -d $DAYS -1

# 2 - purge_histories.sh
sudo -u galaxy PYTHON_EGG_CACHE=$PYTHON_EGG_CACHE python ./scripts/cleanup_datasets/cleanup_datasets.py universe_wsgi.ini -d $DAYS -2 $FLAGS

# 3 - purge_libraries.sh
sudo -u galaxy PYTHON_EGG_CACHE=$PYTHON_EGG_CACHE python ./scripts/cleanup_datasets/cleanup_datasets.py universe_wsgi.ini -d $DAYS -4 $FLAGS

# 4 - purge_folders.sh
sudo -u galaxy PYTHON_EGG_CACHE=$PYTHON_EGG_CACHE python ./scripts/cleanup_datasets/cleanup_datasets.py universe_wsgi.ini -d $DAYS -5 $FLAGS

# 5 - delete_datasets.sh If it is desired that datasets be removed before their outer container (history, library/library folder) has been deleted, 
# the delete_datasets.sh script can be used before the purge_datasets.sh script. This script may take some time to complete. 
# If this step is not desirable comment the next line
sudo -u galaxy PYTHON_EGG_CACHE=$PYTHON_EGG_CACHE python ./scripts/cleanup_datasets/cleanup_datasets.py universe_wsgi.ini -d $DAYS -6 $FLAGS

# 6 - purge_datasets.sh 
sudo -u galaxy PYTHON_EGG_CACHE=$PYTHON_EGG_CACHE python ./scripts/cleanup_datasets/cleanup_datasets.py universe_wsgi.ini -d $DAYS -3 $FLAGS




