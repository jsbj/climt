# $Header: /Users/rca/cvsroot/CliMT/src/radiation/ccm3/src/crm-2.1.2-ccm-3.6/bld/Makefile,v 1.2 2005/02/17 23:33:44 rca Exp $

# Purpose: Makefile for CRM
# Allow for precedence of user-modified files, different search paths
# This Makefile requires GNU make, AT&T Berkeley-style make will not parse it

# Known problems: 
# 19990907: GNU make 3.77 has a bug in wildcard function on Solaris
# Solution: Use a different version of GNU make on that platform

# 19990907: Solaris f90 does not support -r8 -i4 so netCDF output is garbage
# Workaround: Use f77 on Solaris and #undef CRM_SRB in misc.h to remove Fortran90 code

MY_BIN_DIR := ../bin
MY_BLD_DIR := ../bld
MY_DAT_DIR := ../data
MY_DOC_DIR := ../doc
MY_DPN_DIR := ../obj
MY_INC_DIR := 
MY_LIB_DIR := 
MY_OBJ_DIR := $(MY_DPN_DIR)
MY_SRC_DIR := ../src
ifndef OPTS
 OPTS := O
endif
ifndef PVM_ARCH
 PVM_ARCH := $(shell ../bld/pvmgetarch)
endif
ifneq (,$(findstring SUN,$(PVM_ARCH)))
 CCM_TRG_FLG := SUN
endif
ifneq (,$(findstring CRAY,$(PVM_ARCH)))
 CCM_TRG_FLG := CRAY
endif
ifneq (,$(findstring SGI,$(PVM_ARCH)))
 CCM_TRG_FLG := SGI
endif
ifneq (,$(findstring LINUX,$(PVM_ARCH)))
 CCM_TRG_FLG := SGI
endif
ifndef NETCDF_INC
 NETCDF_INC := /usr/local/include # Directory containing netcdf.h
endif
ifndef NETCDF_LIB
 NETCDF_LIB := /usr/local/lib # Directory containing libnetcdf.a
endif

# NB: Do NOT add comment lines, e.g., # This is a comma, to character definitions
null :=
space := $(null) $(null)
comma := ,
newline := \n
# We replace `/' by `cszzsc', call perl, then replace `cszzxc' by `/'
# Unique character(s) to substitute for $(slash) before passing to perl regex
slash_rx := cszzsc
# Unique character(s) to replace by $(slash_rx) before passing to perl regex
slash := /
MY_OBJ_DIR_RX := $(subst $(slash),$(slash_rx),$(MY_OBJ_DIR))
MY_DPN_DIR_RX := $(subst $(slash),$(slash_rx),$(MY_DPN_DIR))

# Directories to search for source files
# misc.h params.h preproc.h are in /fs/cgd/data0/zender/ccm_dst/src/dst
MDL_PTH := \
	$(MY_SRC_DIR)/crm \
	$(MY_SRC_DIR)/ccmlsm_share \
	$(MY_SRC_DIR)/control \
	$(MY_SRC_DIR)/csm_share \
	$(MY_SRC_DIR)/dom \
	$(MY_SRC_DIR)/eul \
	$(MY_SRC_DIR)/physics \
	$(MY_SRC_DIR)/srchutil
# Find all C++, C, and Fortran files in a given directory
FIND_FNC = $(wildcard $(dir)/*.cc $(dir)/*.c $(dir)/*.F)
# 19990907: make-3.77 wildcard function does not work correctly on Solaris
FOO_LST = $(wildcard ~/.* )
# Assemble source files from all directories
SRC_LST = $(foreach dir, $(MDL_PTH),$(FIND_FNC))
# Source file names with directories removed
MDL_SRC := $(notdir $(SRC_LST))
# Dependency list for executable
MDL_OBJ := $(addprefix $(MY_OBJ_DIR)/,$(addsuffix .o, $(basename $(MDL_SRC)))) 
# Dependency (make) file for each object file
MDL_DPN := $(addprefix $(MY_DPN_DIR)/,$(addsuffix .d, $(basename $(MDL_SRC)))) 
# VPATH helps make find dependencies (which are not pathname qualified) in *.d file
VPATH := $(subst $(space),:,$(MDL_PTH))
# Prepend -I to use for compiler argument
CPP_PTH := $(foreach dir,$(MDL_PTH),-I$(dir))

# Variables having to do with binary executables created by module
MDL_BIN_TRG := crm # Binary targets
MDL_BIN_SYM_LNK := # Symbolic links
MDL_BIN_STB := $(MDL_BIN_TRG) $(MDL_BIN_SYM_LNK) # All C++ files in MY_BIN_DIR
MDL_BIN := $(addprefix $(MY_BIN_DIR)/,$(MDL_BIN_STB)) # dst_cln removes these files

# Variables having to do with header files created by module
MDL_INC_TRG := 
MDL_INC_SYM_LNK := # Symbolic links
MDL_INC_STB := $(MDL_INC_TRG) $(MDL_INC_SYM_LNK) # All header files in MY_INC_DIR
MDL_INC := $(addprefix $(MY_INC_DIR)/,$(MDL_INC_STB)) # dst_cln removes these files

# Redefine default C and C++ pattern rules
$(MY_OBJ_DIR)/%.o : %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $(MY_OBJ_DIR)/$(notdir $@)
$(MY_OBJ_DIR)/%.o : %.cc
	$(C++) $(CPPFLAGS) $(C++FLAGS) -c $< -o $(MY_OBJ_DIR)/$(notdir $@)

# Default Fortran pattern rules: CRAY and RS6K must override these rules
$(MY_OBJ_DIR)/%.o : %.F
	$(FC) $(CPPFLAGS) -c $(FFLAGS) -o $(MY_OBJ_DIR)/$(notdir $@) $<
$(MY_OBJ_DIR)/%.o : %.f
	$(FC) -c $(FFLAGS) -o $(MY_OBJ_DIR)/$(notdir $@) $<

# Automatically generate a dependency file for each source file
# $* is the stem, e.g., f
# $@ is the filename of the target, e.g., f.d
# Linux gcc may return an extra `.F' on Fortran names, e.g., `hello.F.o: hello.F'
# (.F)? gets rid of this extra `.F'
$(MY_DPN_DIR)/%.d : %.F
# Following command makes, e.g., f.d begin "f.o f.d : f.F ..."
# Since f.o is not preceded by $(MY_OBJ_DIR), objects are not recompiled when sources are touched
#	$(CPP) -M $(CPPFLAGS) $< | perl -p -e 's/$*\.F\.o/$*.o $@/g;' > $@
# Following command makes, e.g., f.d begin "/home/zender/obj/LINUX/f.o f.d : f.F ..."
# This works fairly well, but is a hack
# First pattern substitutes MY_OBJ_DIR_RX, which has placeholders for slashes
# Second pattern substitutes slashes for the placeholders
	$(CPP) -M $(CPPFLAGS) $< | perl -p -e 's/$*(\.F)?\.o/$(MY_OBJ_DIR_RX)\/$*.o $(MY_DPN_DIR_RX)\/$(notdir $@)/g;s/$(slash_rx)/\$(slash)/g' > $@
# Following command makes, e.g., f.d begin "$(MY_OBJ_DIR)/f.o f.d : f.F ..."
# This would be the ideal command but I can't get the dollar sign to show up
#	$(CPP) -M $(CPPFLAGS) $< | perl -p -e 's/$*\.F\.o/\$(dollar)MY_OBJ_DIR\/$*.o $@/g;' > $@

$(MY_DPN_DIR)/%.d : %.c
#	$(CPP) -M $(CPPFLAGS) $< | perl -p -e 's/$*\.o/$*.o $@/g;' > $@
	$(CPP) -M $(CPPFLAGS) $< | perl -p -e 's/$*\.o/$(MY_OBJ_DIR_RX)\/$*.o $(MY_DPN_DIR_RX)\/$(notdir $@)/g;s/$(slash_rx)/\$(slash)/g' > $@

$(MY_DPN_DIR)/%.d : %.cc
#	$(CPP) -M $(CPPFLAGS) $< | perl -p -e 's/$*\.o/$*.o $@/g;' > $@
	$(CPP) -M $(CPPFLAGS) $< | perl -p -e 's/$*\.o/$(MY_OBJ_DIR_RX)\/$*.o $(MY_DPN_DIR_RX)\/$(notdir $@)/g;s/$(slash_rx)/\$(slash)/g' > $@

ifeq ($(PVM_ARCH),CRAY)
C++ := g++
CC := cc
CPP := cpp
CPPFLAGS := -D$(CCM_TRG_FLG) $(CPP_PTH) -I$(NETCDF_INC)
FC := f90
LD := ld
LDFLAGS := $(MY_LIB_DIR) -L$(NETCDF_LIB) -lnetcdf -lm
LEX := lex
LINT := lint
YACC := yacc
ifeq ($(OPTS),O)
 CFLAGS := -O
 FFLAGS = -N 132
endif
ifeq ($(OPTS),D)
 CFLAGS := -g 
 FFLAGS = -g -N 132
endif
ifeq ($(OPTS),X)
 CFLAGS := -g -N 132
 FFLAGS = -g -N 132 -e i
endif
# 1997/10/21 Added -P to suppress #line # directives on Fortran files
# 1998/09/07 Added -N to prevent CPP from producing .f that complain about single quotes ('s)
# 1998/09/08 Added -I$(NETCDF_INC) to $(FC) for CCM code (e.g., initext.F) that uses Fortran includes
# NB: -N and -P are not part of CPPFLAGS because they should only used for Fortran files
$(MY_OBJ_DIR)/%.o : %.F
	$(CPP) -N -P $(CPPFLAGS) $< > $(MY_OBJ_DIR)/$(patsubst %.F,%.f,$(notdir $<))
	$(FC) -I$(NETCDF_INC) $(FFLAGS) -b $@ $(MY_OBJ_DIR)/$(patsubst %.F,%.f,$(notdir $<))
$(MY_OBJ_DIR)/%.o : %.f
	$(FC) -c $(FFLAGS) $<
	mv -f $(notdir $@) $(MY_OBJ_DIR)
endif
# endif CRAY

ifeq ($(PVM_ARCH),LINUX)
C++ := g++
CC := gcc -ansi
CPP := gcc
CPPFLAGS := -D$(CCM_TRG_FLG) $(CPP_PTH) -I$(NETCDF_INC)
FC := pgf90
LD := ld
LDFLAGS := $(MY_LIB_DIR) -L$(NETCDF_LIB) -lnetcdf -lm
LEX := lex
LINT := lint
YACC := yacc
ifeq ($(CC),gcc)
ifeq ($(OPTS),O)
 CFLAGS := -O -Wall
endif
ifeq ($(OPTS),D)
 CFLAGS := -g -Wall
endif
ifeq ($(OPTS),R)
 CFLAGS := -Wall
endif
ifeq ($(OPTS),X)
 CFLAGS := -g -O -Wall
endif
 C++FLAGS := $(CFLAGS)
endif
ifeq ($(FC),pgf90)
# NB: -Mrecursive ensures that memory is allocated on the stack
ifeq ($(OPTS),O)
 FFLAGS := -fast -mp -Mr8 -Mi4 -byteswapio -Mextend -Mnosecond_underscore -Mrecursive -Mdalign
endif
ifeq ($(OPTS),D)
 FFLAGS := -g -mp -Mr8 -Mi4 -byteswapio -Mextend -Mnosecond_underscore -Mrecursive -Mdalign
endif
ifeq ($(OPTS),R)
 FFLAGS := -mp -Mr8 -Mi4 -byteswapio -Mextend -Mnosecond_underscore -Mrecursive -Mdalign
endif
ifeq ($(OPTS),X)
 FFLAGS := -g -Mbounds -mp -Mr8 -Mi4 -byteswapio -Mextend -Mnosecond_underscore -Mrecursive -Mdalign
endif
endif
ifeq ($(FC),g77)
ifeq ($(OPTS),O)
 FFLAGS := -O -ffixed-line-length-132 -fno-second-underscore
endif
ifeq ($(OPTS),D)
 FFLAGS := -g -ffixed-line-length-132 -fno-second-underscore -fdebug-kludge
endif
ifeq ($(OPTS),R)
 FFLAGS := -ffixed-line-length-132 -fno-second-underscore -fdebug-kludge
endif
ifeq ($(OPTS),X)
 FFLAGS := -g -O -ffixed-line-length-132 -fno-second-underscore -fdebug-kludge
endif
endif
endif
# endif LINUX

ifeq ($(PVM_ARCH),RS6K)
C++ := g++
CC := gcc -ansi
CPP := /lib/cpp -P
CPPFLAGS := -D$(CCM_TRG_FLG) $(CPP_PTH) -I$(NETCDF_INC)
FC := IBM
LD := ld
LDFLAGS := $(MY_LIB_DIR) -L$(NETCDF_LIB) -lncaru -lnetcdf -lm
LEX := lex
LINT := lint
YACC := yacc
ifeq ($(OPTS),O)
 CFLAGS := -O2
 PREPROCESS.F := $(CPP) $(CPPFLAGS)
 FFLAGS := -O -NS2000 -qfixed=132
endif
ifeq ($(OPTS),D)
 CFLAGS := -g
 PREPROCESS.F := $(CPP) $(CPPFLAGS)
 FFLAGS := -g -NS2000 -qfixed=132
endif
 FFLAGS += -qREALSIZE=8 -qINTSIZE=4
$(MY_OBJ_DIR)/%.o : %.F
	$(CPP) $(CPPFLAGS) $< $(MY_OBJ_DIR)/$(basename $<).f 
	$(FC) -c $(FFLAGS) -o $(MY_OBJ_DIR)/$(notdir $@) $(MY_OBJ_DIR)/$(basename $<).f
$(MY_OBJ_DIR)/%.o : %.f
	$(FC) -c $(FFLAGS) -o $(MY_OBJ_DIR)/$(notdir $@) $<
endif
# endif RS6K

ifeq ($(PVM_ARCH),SGI5)
C++ := g++
CC := gcc -ansi
CPP := cpp
CPPFLAGS := -D$(CCM_TRG_FLG) $(CPP_PTH) -I$(NETCDF_INC)
FC := f77
LD := ld
LDFLAGS := $(MY_LIB_DIR) -L$(NETCDF_LIB) -lnetcdf -lm
LEX := lex
LINT := lint
YACC := yacc
ifeq ($(OPTS),O)
 CFLAGS := -O2
 FFLAGS := -O2 -e -extend_source
endif
ifeq ($(OPTS),D)
 CFLAGS := -g
 FFLAGS := -g -e -extend_source
endif
endif
# endif SGI5

ifeq ($(PVM_ARCH),SGI64)
C++ := c++
CC := cc
CPP := cpp
CPPFLAGS := -D$(CCM_TRG_FLG) $(CPP_PTH) -I$(NETCDF_INC)
FC := f90 -cpp
LD := ld
LDFLAGS := -64 $(MY_LIB_DIR) -L$(NETCDF_LIB) -lnetcdf
LEX := flex
LINT := lint
YACC := bison
ifeq ($(OPTS),O)
 CFLAGS := -64 -O2
 FFLAGS := -64 -O2 -extend_source
endif
ifeq ($(OPTS),D)
 CFLAGS := -64 -g
 FFLAGS := -64 -g -extend_source
endif
ifeq ($(OPTS),X)
 CFLAGS := -64 -g -trapuv
 FFLAGS := -64 -g -extend_source -check_bounds -trapuv
endif
 FFLAGS += -r8 -i4
endif
# endif SGI64

ifeq ($(PVM_ARCH),SGIMP64)
C++ := CC
CC := cc
CPP := cpp
CPPFLAGS := -D$(CCM_TRG_FLG) $(CPP_PTH) -I$(NETCDF_INC)
FC := f90 -cpp
LD := ld
LDFLAGS := -64 -mips4 $(MY_LIB_DIR) -L$(NETCDF_LIB) -lnetcdf
LEX := flex
LINT := lint
YACC := bison
ifeq ($(OPTS),O)
 CFLAGS := -64 -mips4 -O2
 FFLAGS := -64 -mips4 -O2 -extend_source
endif
ifeq ($(OPTS),R)
 CFLAGS := -64 -mips4
 FFLAGS := -64 -mips4 -extend_source
endif
ifeq ($(OPTS),D)
 CFLAGS := -64 -mips4 -g
 FFLAGS := -64 -mips4 -g -extend_source
endif
ifeq ($(OPTS),X)
 CFLAGS := -64 -mips4 -g -trapuv
 FFLAGS := -64 -mips4 -g -extend_source -check_bounds -trapuv
endif
 FFLAGS += -r8 -i4
endif
# endif SGIMP64

ifeq ($(PVM_ARCH),SUN4)
C++ := g++
CC := acc
CPP := cpp
CPPFLAGS := -D$(CCM_TRG_FLG) $(CPP_PTH) -I$(NETCDF_INC)
FC := f77
LD := ld
LDFLAGS := $(MY_LIB_DIR) -L$(NETCDF_LIB) -lnetcdf
LEX := lex
LINT := lint
YACC := yacc
ifeq ($(OPTS),O)
 CFLAGS := -O2
 FFLAGS := -fast -e -Nl99
endif
ifeq ($(OPTS),D)
 CFLAGS := -g
 FFLAGS := -g -e -Nl99
endif
 FFLAGS += -r8 -i4
endif
# endif SUN4

ifeq ($(PVM_ARCH),SUN4SOL2)
C++ := g++
CC := gcc -ansi
CPP := cpp
CPPFLAGS := -D$(CCM_TRG_FLG) $(CPP_PTH) -I$(NETCDF_INC)
FC := f77
LD := ld
LDFLAGS := $(MY_LIB_DIR) -L$(NETCDF_LIB) -lnetcdf
LEX := lex
LINT := lint
YACC := yacc
ifeq ($(OPTS),O)
 CFLAGS := -O2
 FFLAGS := -O -e
endif
ifeq ($(OPTS),D)
 CFLAGS := -g
 FFLAGS := -g -e
endif
ifeq ($(OPTS),X)
 CFLAGS := -g 
 FFLAGS := -g -e
# NB: 1998/06/01 -C (range-checking) is not supported by Sun f90
ifeq ($(FC),f77)
 FFLAGS += -C -r8 -i4
endif
endif
# NB: 1999/09/07 -r8 -i4 are not supported by Sun f90, Sun f77 now uses -dbl
ifeq ($(FC),f77)
 FFLAGS += -r8 -i4
endif
endif
# endif SUN4SOL2

ifeq ($(PVM_ARCH),SUNMP)
C++ := g++
CC := gcc -ansi
CPP := cpp
CPPFLAGS := -D$(CCM_TRG_FLG) $(CPP_PTH) -I$(NETCDF_INC)
FC := f90
LD := ld
LDFLAGS := $(MY_LIB_DIR) -L$(NETCDF_LIB) -lnetcdf
LEX := lex
LINT := lint
YACC := yacc
ifeq ($(OPTS),O)
 CFLAGS := -O2
 FFLAGS := -O -e
endif
ifeq ($(OPTS),D)
 CFLAGS := -g
 FFLAGS := -g -e
endif
ifeq ($(OPTS),X)
 CFLAGS := -g 
 FFLAGS := -g -e
# 1998/06/01 -C (range-checking) is not supported by Sun f90
ifeq ($(FC),f77)
 FFLAGS += $(FFLAGS) -C
endif
endif
# 1999/09/07 -r8 -i4 are not supported by Sun f90
ifeq ($(FC),f77)
 FFLAGS += -r8 -i4
endif
endif
# endif SUNMP

all : $(MDL_BIN_TRG)
# .PHONY tells make to remake the following non-file targets 
.PHONY: all cln dst_cln debug 
# Delete default suffixes---this should increase speed
.SUFFIXES: 
# Define suffixes which matter
.SUFFIXES: .cc .c .o .F .d

# Targets in bin
CRM_OBJ := $(MDL_OBJ)
crm:	$(CRM_OBJ)
	$(FC) -o $(MY_BIN_DIR)/$@ $(MDL_OBJ) $(LDFLAGS)
	chmod 755 $(MY_BIN_DIR)/$@
bin : $(MDL_BIN_TRG)
bin_cln: 
	rm -f $(MDL_BIN)
# Targets in bld
bld_cln :
	cd $(MY_BLD_DIR); rm -f TAGS
clean : cln
cln : dat_cln dpn_cln obj_cln
debug: dbg
dbg:
	printf "PVM_ARCH = $(PVM_ARCH)\n"
	printf "LDFLAGS = $(LDFLAGS)\n"
	printf "MDL_PTH = $(MDL_PTH)\n"
	printf "VPATH = $(VPATH)\n"
	printf "FOO_LST = $(FOO_LST)\n"
	printf "SRC_LST = $(SRC_LST)\n"
	printf "MDL_SRC = $(MDL_SRC)\n"
	printf "MDL_OBJ = $(MDL_OBJ)\n"
	printf "MDL_DPN = $(MDL_DPN)\n"
distclean : dst_cln
dst_cln : cln bld_cln
install : inc 
tags : 
	etags $(SRC_LST) Makefile
uninstall : dst_cln bin_cln inc_cln lib_cln
# Targets in dat
dat_cln :
	-cd $(MY_DAT_DIR); rm -f foo.nc
# Targets in dpn
dpn : $(MDL_DPN)
dpn_cln :
	rm -f $(MDL_DPN)
# Targets in inc
inc : $(MDL_INC)
inc_cln :
	rm -f $(MDL_INC)
# Targets in lib
lib : 
	-cd $(MY_LIB_DIR)
lib_cln :
	-cd $(MY_LIB_DIR)
# Targets in obj
obj : $(MDL_OBJ)
obj_cln :
	rm -f $(MDL_OBJ)

# Create dependency files only if they will not be immediately deleted
INCLUDE_DPN := TRUE
GOALS_WHICH_DELETE_DEPENDENCY_FILES := cln clean dir distclean dst_cln dpn_cln tags uninstall
ifeq (,$(findstring $(MAKECMDGOALS),$(GOALS_WHICH_DELETE_DEPENDENCY_FILES)))
 INCLUDE_DPN := TRUE
else
 INCLUDE_DPN := FALSE
endif
ifeq ($(INCLUDE_DPN),TRUE)
-include $(MDL_DPN)
endif

