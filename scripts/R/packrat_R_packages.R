# Step 1: Install and load packrat
#install.packages("packrat")
setwd("/research/labs/neurology/fryer/m239830/LBD_CWOW/scripts/R/package_manager")
library(packrat)

# Step 2: Initialize packrat
# I created a new directory to store libraries (ran into issues otherwise)
packrat::init(project = "/research/labs/neurology/fryer/m239830/LBD_CWOW/scripts/R/package_manager")

# Step 3: Turn packrat on/off, see how library path changes
packrat::on()
#.libPaths()
#packrat::off()
#.libPaths()

# Step 4: Take a snapshot to save changes
packrat::snapshot()
# initiat
packrat::init()
#packrat::bundle()

#--------- Notes 
# Basic Commands
# Remove packages
#remove.packages("plyr")
# Install local source packages
#packrat::set_opts(local.repos = "<path_to_repo>")
# Save the current state of your library
#packrat::snapshot()
# Restore the library state saved in the most recent snapshot
#packrat::restore()
# Remove unused packages from your library
#packrat::clean()
# Bundle a packrat project, for easy sharing
#packrat::bundle()
# Unbundle a packrat project, generating a project directory with libraries restored from the most recent snapshot
#packrat::unbundle()
