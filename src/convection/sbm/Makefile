# Makes  Python module 

# Rodrigo Caballero Augi (rca@geosci.uchicago.edu), 
# Department of the Geophysical Sciences, University of Chicago
# 2004

name := sbm_convection

OS := $(shell uname -s)
include $(CLIMT_SRC_DIR)/CompilerFlags.mk

# local flags
INCLUDE   := -I./include -I./src
FCFLAGS   := $(FCFLAGS) $(INCLUDE) 
F90FLAGS  := $(F90FLAGS) $(INCLUDE) 
CPPFLAGS  := -DIM=$(IM) -DKM=$(KM)

include $(CLIMT_SRC_DIR)/Rules.mk
