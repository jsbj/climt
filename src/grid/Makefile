# Makes Python module 

# Rodrigo Caballero Augi (rca@geosci.uchicago.edu), 
# Department of the Geophysical Sciences, University of Chicago
# 2004

name := grid

OS := $(shell uname -s)

include $(CLIMT_SRC_DIR)/CompilerFlags.mk
# add local flags
CPPFLAGS  := -DIM=$(IM) -DJM=$(JM) -DKM=$(KM)

include $(CLIMT_SRC_DIR)/Rules.mk
