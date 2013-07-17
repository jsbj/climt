# Makes  Python module 

# Rodrigo Caballero Augi (rca@geosci.uchicago.edu), 
# Department of the Geophysical Sciences, University of Chicago
# 2004

name := greygas_radiation

OS := $(shell uname -s)

include $(CLIMT_SRC_DIR)/CompilerFlags.mk
# add local flags
INCLUDE   := -I./src
FCFLAGS   := $(FCFLAGS) $(INCLUDE) 
F90FLAGS  := $(F90FLAGS) $(INCLUDE) 

include $(CLIMT_SRC_DIR)/Rules.mk
