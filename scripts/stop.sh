#!/bin/sh

# define environmental variables and set program paths

# galaxy installation directory, all data related to the specific galaxy installation will
# be stored  here e.g. galaxy_dist, uct_cbio_tools, uct_cbio_tools_data
export GALAXY_INSTALL_DIR=/usr/local/share/galaxy

# change to galaxy dist directory (this is where the run.sh script is)
export GALAXY_DIST_DIR=$GALAXY_INSTALL_DIR/galaxy_dist
cd $GALAXY_DIST_DIR

# stop daemon
sudo -u galaxy sh run.sh --stop-daemon
