# Makes  Python module 

# Rodrigo Caballero Augi (rca@geosci.uchicago.edu), 
# Department of the Geophysical Sciences, University of Chicago
# 2004

name := axisymmetric_dynamics

OS := $(shell uname -s)
include $(CLIMT_SRC_DIR)/CompilerFlags.mk

# local flags
INCLUDE   := -I/sw/include
FCFLAGS   := $(FCFLAGS) $(INCLUDE) 
F90FLAGS  := $(F90FLAGS) $(INCLUDE) 

include $(CLIMT_SRC_DIR)/Rules.mk
