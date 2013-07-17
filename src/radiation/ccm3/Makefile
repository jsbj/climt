# Makes Python module 

# Rodrigo Caballero Augi (rca@geosci.uchicago.edu), 
# Department of the Geophysical Sciences, University of Chicago
# 2003

name := ccm3_radiation

OS := $(shell uname -s)
include $(CLIMT_SRC_DIR)/CompilerFlags.mk

# local flags
INCLUDE   := -I./include -I./src
FCFLAGS   := $(FCFLAGS) $(INCLUDE) 
F90FLAGS  := $(F90FLAGS) $(INCLUDE) 
CPPFLAGS  := -DSUN -DPLON=1 -DPLEV=$(KM) -DPLEVR=$(KM)

include $(CLIMT_SRC_DIR)/Rules.mk
