# Makes  Python module 

# Rodrigo Caballero Augi (rca@geosci.uchicago.edu), 
# Department of the Geophysical Sciences, University of Chicago
# 2004

name := ccm3_convection

OS := $(shell uname -s)
include $(CLIMT_SRC_DIR)/CompilerFlags.mk

# local flags
INCLUDE   := -I./include -I./src
FCFLAGS   := $(FCFLAGS) $(INCLUDE) 
F90FLAGS  := $(F90FLAGS) $(INCLUDE) 
CPPFLAGS  := -DIM=$(IM) -DKM=$(KM) -DALPHA -DPLON=$(IM) -DPLEV=$(KM) \
            -DPLEVR=$(KM) -DPCNST=1 -DPLAT=1 -DPTRM=42 -DPTRN=42 -DPTRK=42 -DPNATS=0

include $(CLIMT_SRC_DIR)/Rules.mk
