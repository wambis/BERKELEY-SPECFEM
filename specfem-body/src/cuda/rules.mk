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

#######################################

cuda_TARGETS = \
	$(cuda_OBJECTS) \
	$(EMPTY_MACRO)

cuda_OBJECTS = \
	$O/assemble_MPI_scalar_cuda.cuda.o \
	$O/assemble_MPI_vector_cuda.cuda.o \
	$O/check_fields_cuda.cuda.o \
	$O/compute_add_sources_elastic_cuda.cuda.o \
	$O/compute_coupling_cuda.cuda.o \
	$O/compute_forces_crust_mantle_cuda.cuda.o \
	$O/compute_forces_inner_core_cuda.cuda.o \
	$O/compute_forces_outer_core_cuda.cuda.o \
	$O/compute_kernels_cuda.cuda.o \
	$O/compute_stacey_acoustic_cuda.cuda.o \
	$O/compute_stacey_elastic_cuda.cuda.o \
	$O/initialize_cuda.cuda.o \
	$O/noise_tomography_cuda.cuda.o \
	$O/prepare_mesh_constants_cuda.cuda.o \
	$O/transfer_fields_cuda.cuda.o \
	$O/update_displacement_cuda.cuda.o \
	$O/write_seismograms_cuda.cuda.o \
	$O/save_and_compare_cpu_vs_gpu.cudacc.o \
	$(EMPTY_MACRO)

cuda_STUBS = \
	$O/specfem3D_gpu_cuda_method_stubs.cudacc.o \
	$(EMPTY_MACRO)

cuda_DEVICE_OBJ = \
	$O/cuda_device_obj.o \
	$(EMPTY_MACRO)


#######################################

## compilation directories
S := ${S_TOP}/src/cuda
$(cuda_OBJECTS): S = ${S_TOP}/src/cuda

####
#### rule for each .o file below
####

###
### CUDA compilation
###

$O/%.cuda.o: $S/%.cu ${SETUP}/config.h $S/mesh_constants_cuda.h $S/prepare_constants_cuda.h
	${NVCC} -c $< -o $@ $(NVCC_FLAGS)

$O/%.cudacc.o: $S/%.c ${SETUP}/config.h
	${CC} -c $(CPPFLAGS) $(CFLAGS) $(MPI_INCLUDES) -o $@ $<

