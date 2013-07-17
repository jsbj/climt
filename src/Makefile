
# Makes and installs CliMT

# Rodrigo Caballero Augi (rca@geosci.uchicago.edu), 
# Department of the Geophysical Sciences, University of Chicago
# 2005

export

# --- User adjustable parameters:---------------------

# Some CliMT components require dimensioning information at 
# compile time. Specific components with this requirement are
# listed below in parentheses.

# No. longitude grid points 
# (Placeholder; currently not needed by any component)
IM := 1

# No. latitude grid points
# (Placeholder; currently not needed by any component)
JM := 1

# No. of levels in vertical
# (ccm3_radiation, chou_radiation, ccm3_turbulence)
KM := 30

# Root directory where you want CliMT to be installed
# (e.g. /usr/lib/python2.3/site-packages for a site-wide installation)
CLIMT_DIR  := $(HOME)/lib/my_python

# location of Python include files
# (try 'locate Python.h' to determine this directory)
PYTHONINC  := /Library/Frameworks/Python.framework/Headers/

# location of Numeric include files
# (try 'locate arrayobject.h' to determine this directory)
NUMPYINC  := /Library/Frameworks/Python.framework/Versions/2.5/lib/python2.5/site-packages/numpy/core/include/

# Fortran compiler: must be one of IBM, Intel,  PGI, GNU
FCOMPILER  :=GNU

# Full path to the python you'll be using
PYTHON    := /usr/local/bin/python

# --- No changes necessary below this line ------------------------

F2PY          := /usr/local/bin/f2py
INSTALL_DIR   := $(CLIMT_DIR)/climt
sharedlibs    := `find . -name "*.so"`
CLIMT_SRC_DIR := $(PWD)
NUMPYINC      := $(subst /Numeric,,$(NUMPYINC))

# list of all dirs that will be Made
SUBDIRS := grid thermodyn timestep \
           radiation/chou \
           radiation/ccm3 \
           radiation/greygas \
           radiation/ozone \
           radiation/insolation \
           convection/emanuel \
           convection/hard \
           convection/sbm \
           ocean/slab_ocean \
           dynamics/axisymmetric \
           turbulence/ccm3 \
           turbulence/simple \

.PHONY: all test $(SUBDIRS)

all: $(SUBDIRS)

install: all
	mkdir -p $(INSTALL_DIR)
	cp -f $(sharedlibs) ../lib/climt/*.py $(INSTALL_DIR)

$(SUBDIRS):
	$(MAKE) -C $@

test: all
	for dir in $(SUBDIRS); do  \
	  $(PYTHON) $$dir/test.py; \
	done
clean:
	find . \( -name "*.o" -o -name "*.mod" -o -name "*.so" \) | xargs rm -f
	find . \( -name "*.lst" -o -name "*.ipo"  \) | xargs rm -f
	find . \( -name "*_lib.*" -o -name "*~" -o -name ".*~" \) | xargs rm -f
	find .. \( -name "*.pyc" -o -name "*.so" -o -name ".*~" \) | xargs rm -f

