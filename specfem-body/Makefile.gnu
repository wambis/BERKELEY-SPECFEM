#=====================================================================
#
#          S p e c f e m 3 D  G l o b e  V e r s i o n  6 . 0
#          --------------------------------------------------
#
#     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
#                        Princeton University, USA
#                and CNRS / University of Marseille, France
#                 (there are currently many more authors!)
# (c) Princeton University and CNRS / University of Marseille, April 2014
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
#=====================================================================

# Makefile.  Generated from Makefile.in by configure.

#######################################

FC = gfortran
FCFLAGS = #-g -O2
MPIFC = mpifort
MPILIBS = 

FLAGS_CHECK = -std=gnu -frange-check -O2 -fmax-errors=10 -pedantic -pedantic-errors -Waliasing -Wampersand -Wcharacter-truncation -Wline-truncation -Wsurprising -Wno-tabs -Wunderflow -ffree-line-length-512
#FLAGS_CHECK = -std=gnu -frange-check -O2 -fmax-errors=10 -pedantic -pedantic-errors -Waliasing -Wampersand -Wsurprising -Wno-tabs -Wunderflow
FCFLAGS_f90 = -J./obj -I./obj -I.  -I. -I${SETUP}

FC_MODEXT = mod
FC_MODDIR = ./obj

FCCOMPILE_CHECK = ${FC} ${FCFLAGS} $(FLAGS_CHECK)

MPIFCCOMPILE_CHECK = ${MPIFC} ${FCFLAGS} $(FLAGS_CHECK)

CC = gcc
CFLAGS = -g -O2
CPPFLAGS = -I${SETUP} 

FCLINK = $(MPIFCCOMPILE_CHECK)


#######################################
####
#### CUDA
#### with configure: ./configure --with-cuda=cuda5 CUDA_FLAGS=.. CUDA_LIB=.. CUDA_INC=.. MPI_INC=.. ..
####
#######################################

#CUDA = yes
CUDA = no

#CUDA5 = yes
CUDA5 = no

MPI_INCLUDES = 
CUDA_FLAGS = 
CUDA_INC =  -I${SETUP}
CUDA_LINK =  

#NVCC = nvcc
NVCC = gcc

# GPU architecture

# CUDA architecture / code version
# Fermi: -gencode=arch=compute_10,code=sm_10 not supported
# Tesla (default): -gencode=arch=compute_20,code=sm_20
# Geforce GT 650m: -gencode=arch=compute_30,code=sm_30
# Kepler (cuda5) : -gencode=arch=compute_35,code=sm_35
GENCODE_20 = -gencode=arch=compute_20,code=\"sm_20,compute_20\"
GENCODE_30 = -gencode=arch=compute_30,code=\"sm_30,compute_30\"
GENCODE_35 = -gencode=arch=compute_35,code=\"sm_35,compute_35\"

# CUDA version 5.x
##GENCODE = $(GENCODE_35)
# CUDA version 4.x
#GENCODE = $(GENCODE_20)

# CUDA flags and linking
#NVCC_FLAGS_BASE = $(CUDA_FLAGS) $(CUDA_INC) $(MPI_INCLUDES)
##NVCC_FLAGS = $(NVCC_FLAGS_BASE) -dc $(GENCODE)
#NVCC_FLAGS = $(NVCC_FLAGS_BASE) -DUSE_OLDER_CUDA4_GPU $(GENCODE)

##NVCCLINK_BASE = $(NVCC) $(CUDA_INC) $(MPI_INCLUDES)
##NVCCLINK = $(NVCCLINK_BASE) -dlink $(GENCODE)
#NVCCLINK = $(NVCCLINK_BASE) -DUSE_OLDER_CUDA4_GPU $(GENCODE)

NVCC_FLAGS =
NVCCLINK = $(NVCC) $(NVCC_FLAGS)

#######################################
####
#### VTK
#### with configure: ./configure --enable-vtk ..
####
#######################################

#VTK = yes
VTK = no

CPPFLAGS += 
LDFLAGS += 
MPILIBS += 


#######################################
####
#### ADIOS
#### with configure: ./configure --with-adios ADIOS_CONFIG=..
####
#######################################

#ADIOS = yes
ADIOS = no

FCFLAGS += 
MPILIBS += 


#######################################
####
#### FORCE_VECTORIZATION
#### with configure: ./configure --with-vec ..
####
#######################################
#FORCE_VECTORIZATION = yes
FORCE_VECTORIZATION = no


#######################################
####
#### directories
####
#######################################

## compilation directories
# B : build directory
B = .
# E : executables directory
E = $B/bin
# O : objects directory
O = $B/obj
# S_TOP : source file root directory
S_TOP = .
# setup file directory
SETUP = $B/setup
# output file directory
OUTPUT = $B/OUTPUT_FILES


#######################################
####
#### targets
####
#######################################

# code subdirectories
SUBDIRS = \
	shared \
	create_header_file \
	compute_optimized_dumping_undo_att \
	cuda \
	meshfem3D \
	specfem3D \
	auxiliaries

# default targets
DEFAULT = \
	xcreate_header_file \
	xcompute_optimized_dumping_undo_att \
	xmeshfem3D \
	xspecfem3D \
	xcombine_vol_data \
	xcombine_vol_data_vtk \
	xcombine_surf_data \
	xcombine_AVS_DX \
	xconvolve_source_timefunction \
	xcreate_movie_AVS_DX \
	xcreate_movie_GMT_global \
	$(EMPTY_MACRO)

all: default

default: $(DEFAULT)

backup:
	cp -rp src setup DATA/Par_file* Makefile go_mesher* go_solver* mymachines bak

ifdef CLEAN
clean:
	@echo "cleaning by CLEAN defined"
	-rm -f $(foreach dir, $(CLEAN), $($(dir)_OBJECTS) $($(dir)_MODULES) $($(dir)_SHARED_OBJECTS) $($(dir)_TARGETS)) $O/*
else
clean:
	@echo "cleaning by CLEAN not defined"
	-rm -f $(foreach dir, $(SUBDIRS), $($(dir)_OBJECTS) $($(dir)_MODULES) $($(dir)_TARGETS)) $O/*
endif

help:
	@echo "usage: make [executable]"
	@echo ""
	@echo "supported executables:"
	@echo "    xmeshfem3D"
	@echo "    xspecfem3D"
	@echo "    xcreate_header_file"
	@echo "    xcompute_optimized_dumping_undo_att"
	@echo "    xcombine_vol_data"
	@echo "    xcombine_vol_data_vtk"
	@echo "    xcombine_vol_data_adios"
	@echo "    xcombine_surf_data"
	@echo "    xcombine_AVS_DX"
	@echo "    xconvolve_source_timefunction"
	@echo "    xcreate_movie_AVS_DX"
	@echo "    xcreate_movie_GMT_global"
	@echo ""

.PHONY: all default backup clean help

#######################################


# Get dependencies and rules for building stuff
include $(patsubst %, ${S_TOP}/src/%/rules.mk, $(SUBDIRS))


#######################################

##
## Shortcuts
##

# Shortcut for: <prog>/<xprog> -> bin/<xprog>
define target_shortcut
$(patsubst $E/%, %, $(1)): $(1)
.PHONY: $(patsubst $E/%, %, $(1))
$(patsubst $E/x%, %, $(1)): $(1)
.PHONY: $(patsubst $E/x%, %, $(1))
endef

# Shortcut for: dir -> src/dir/<targets in here>
define shortcut
$(1): $($(1)_TARGETS)
.PHONY: $(1)
$$(foreach target, $$(filter $E/%,$$($(1)_TARGETS)), $$(eval $$(call target_shortcut,$$(target))))
endef

$(foreach dir, $(SUBDIRS), $(eval $(call shortcut,$(dir))))

# Other old shortcuts
bak: backup
mesh: $E/xmeshfem3D
spec: $E/xspecfem3D
.PHONY: bak mesh spec

