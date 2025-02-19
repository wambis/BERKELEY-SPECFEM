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

specfem3D_TARGETS = \
	$E/xspecfem3D \
	$(EMPTY_MACRO)

specfem3D_OBJECTS = \
	$O/assemble_MPI_scalar.solver.o \
	$O/assemble_MPI_vector.solver.o \
	$O/comp_source_spectrum.solver.o \
	$O/comp_source_time_function.solver.o \
	$O/compute_adj_source_frechet.solver.o \
	$O/convert_time.solver.o \
	$O/define_derivation_matrices.solver.o \
	$O/file_io_threads.cc.o \
	$O/get_backazimuth.solver.o \
	$O/get_cmt.solver.o \
	$O/get_event_info.solver.o \
	$O/netlib_specfun_erf.solver.o \
	$O/recompute_jacobian.solver.o \
	$O/ucb_dfour.solver.o \
	$(EMPTY_MACRO)

# solver objects with statically allocated arrays; dependent upon
# values_from_mesher.h
specfem3D_OBJECTS += \
	$O/asdf_data.solverstatic_module.o \
	$O/specfem3D_par.solverstatic_module.o \
	$O/mirror.solverstatic_module.o \
	$O/ucb_heaviside.solverstatic_module.o \
	$O/write_seismograms.solverstatic.o \
	$O/check_stability.solverstatic.o \
	$O/compute_add_sources.solverstatic.o \
	$O/compute_arrays_source.solverstatic.o \
	$O/compute_boundary_kernel.solverstatic.o \
	$O/compute_coupling.solverstatic.o \
	$O/compute_element.solverstatic.o \
	$O/compute_element_att_memory.solverstatic.o \
	$O/compute_element_strain.solverstatic.o \
	$O/compute_forces_acoustic_calling_routine.solverstatic.o \
	$O/compute_forces_viscoelastic_calling_routine.solverstatic.o \
	$O/compute_forces_crust_mantle_noDev.solverstatic.o \
	$O/compute_forces_crust_mantle_Dev.solverstatic.o \
	$O/compute_forces_inner_core_noDev.solverstatic.o \
	$O/compute_forces_inner_core_Dev.solverstatic.o \
	$O/compute_forces_outer_core_noDev.solverstatic.o \
	$O/compute_forces_outer_core_Dev.solverstatic.o \
	$O/compute_kernels.solverstatic.o \
	$O/compute_seismograms.solverstatic.o \
	$O/compute_stacey_crust_mantle.solverstatic.o \
	$O/compute_stacey_outer_core.solverstatic.o \
	$O/finalize_simulation.solverstatic.o \
	$O/get_attenuation.solverstatic.o \
	$O/initialize_simulation.solverstatic.o \
	$O/iterate_time.solverstatic.o \
	$O/iterate_time_undoatt.solverstatic.o \
	$O/locate_receivers.solverstatic.o \
	$O/locate_regular_points.solverstatic.o \
	$O/locate_sources.solverstatic.o \
	$O/multiply_arrays_source.solverstatic.o \
	$O/noise_tomography.solverstatic.o \
	$O/prepare_timerun.solverstatic.o \
	$O/read_adjoint_sources.solverstatic.o \
	$O/read_arrays_solver.solverstatic.o \
	$O/read_forward_arrays.solverstatic.o \
	$O/read_mesh_databases.solverstatic.o \
	$O/read_topography_bathymetry.solverstatic.o \
	$O/save_forward_arrays.solverstatic.o \
	$O/save_kernels.solverstatic.o \
	$O/save_regular_kernels.solverstatic.o \
	$O/setup_GLL_points.solverstatic.o \
	$O/setup_sources_receivers.solverstatic.o \
	$O/specfem3D.solverstatic.o \
	$O/update_displacement_LDDRK.solverstatic.o \
	$O/update_displacement_Newmark.solverstatic.o \
	$O/write_movie_output.solverstatic.o \
	$O/write_movie_volume.solverstatic.o \
	$O/write_movie_surface.solverstatic.o \
	$O/write_output_ASCII.solverstatic.o \
	$O/write_output_SAC.solverstatic.o \
	$(EMPTY_MACRO)

specfem3D_MODULES = \
	$(FC_MODDIR)/asdf_data.$(FC_MODEXT) \
	$(FC_MODDIR)/mirror.$(FC_MODEXT) \
	$(FC_MODDIR)/ucb_heaviside.$(FC_MODEXT) \
	$(FC_MODDIR)/constants_solver.$(FC_MODEXT) \
	$(FC_MODDIR)/specfem_par.$(FC_MODEXT) \
	$(FC_MODDIR)/specfem_par_crustmantle.$(FC_MODEXT) \
	$(FC_MODDIR)/specfem_par_innercore.$(FC_MODEXT) \
	$(FC_MODDIR)/specfem_par_outercore.$(FC_MODEXT) \
	$(FC_MODDIR)/specfem_par_movie.$(FC_MODEXT) \
	$(FC_MODDIR)/write_seismograms_mod.$(FC_MODEXT) \
	$(EMPTY_MACRO)

# These files come from the shared directory
specfem3D_SHARED_OBJECTS = \
	$O/shared_par.shared_module.o \
	$O/auto_ner.shared.o \
	$O/broadcast_computed_parameters.shared.o \
	$O/calendar.shared.o \
	$O/count_elements.shared.o \
	$O/count_number_of_sources.shared.o \
	$O/count_points.shared.o \
	$O/create_name_database.shared.o \
	$O/define_all_layers.shared.o \
	$O/euler_angles.shared.o \
	$O/exit_mpi.shared.o \
	$O/force_ftz.cc.o \
	$O/get_model_parameters.shared.o \
	$O/get_timestep_and_layers.shared.o \
	$O/get_value_parameters.shared.o \
	$O/gll_library.shared.o \
	$O/hex_nodes.shared.o \
	$O/intgrl.shared.o \
	$O/lagrange_poly.shared.o \
	$O/make_ellipticity.shared.o \
	$O/make_gravity.shared.o \
	$O/model_prem.shared.o \
	$O/model_topo_bathy.shared.o \
	$O/parallel.sharedmpi.o \
	$O/param_reader.cc.o \
	$O/read_compute_parameters.shared.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$O/reduce.shared.o \
	$O/rthetaphi_xyz.shared.o \
	$O/spline_routines.shared.o \
	$O/write_c_binary.cc.o \
	$O/write_VTK_file.shared.o \
	$(EMPTY_MACRO)

###
### CUDA
###

cuda_specfem3D_OBJECTS = \
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

cuda_specfem3D_STUBS = \
	$O/specfem3D_gpu_cuda_method_stubs.cudacc.o \
	$(EMPTY_MACRO)

cuda_specfem3D_DEVICE_OBJ = \
	$O/cuda_device_obj.o \
	$(EMPTY_MACRO)

ifeq ($(CUDA),yes)
specfem3D_OBJECTS += $(cuda_specfem3D_OBJECTS)
ifeq ($(CUDA5),yes)
specfem3D_OBJECTS += $(cuda_specfem3D_DEVICE_OBJ)
endif
else
specfem3D_OBJECTS += $(cuda_specfem3D_STUBS)
endif

###
### ADIOS
###

adios_specfem3D_OBJECTS = \
	$O/read_arrays_solver_adios.solverstatic_adios.o \
	$O/read_attenuation_adios.solverstatic_adios.o \
	$O/read_forward_arrays_adios.solverstatic_adios.o \
	$O/read_mesh_databases_adios.solverstatic_adios.o \
	$O/save_forward_arrays_adios.solverstatic_adios.o \
	$O/save_kernels_adios.solverstatic_adios.o \
	$O/write_specfem_adios_header.solverstatic_adios.o \
	$O/write_output_ASDF.solverstatic_adios.o \
	$(EMPTY_MACRO)

adios_specfem3D_SHARED_OBJECTS = \
	$O/adios_helpers_definitions.shared_adios_module.o \
	$O/adios_helpers_writers.shared_adios_module.o \
	$O/adios_helpers.shared_adios.o \
	$O/adios_manager.shared_adios.o \
	$O/asdf_helpers_definitions.shared_adios_module.o \
	$O/asdf_helpers_writers.shared_adios_module.o \
	$O/asdf_helpers.shared_adios.o \
	$(EMPTY_MACRO)

adios_specfem3D_STUBS = \
	$(EMPTY_MACRO)

adios_specfem3D_SHARED_STUBS = \
	$O/adios_method_stubs.cc.o \
	$(EMPTY_MACRO)

# conditional adios linking
ifeq ($(ADIOS),yes)
specfem3D_OBJECTS += $(adios_specfem3D_OBJECTS)
specfem3D_SHARED_OBJECTS += $(adios_specfem3D_SHARED_OBJECTS)
else
specfem3D_OBJECTS += $(adios_specfem3D_STUBS)
specfem3D_SHARED_OBJECTS += $(adios_specfem3D_SHARED_STUBS)
endif

###
### VTK
###

vtk_specfem3D_OBJECTS = \
  $O/visual_vtk.visualcc.o \
	$(EMPTY_MACRO)

vtk_specfem3D_STUBS = \
	$O/visual_vtk_stubs.visualc.o \
	$(EMPTY_MACRO)

# conditional adios linking
ifeq ($(VTK),yes)
specfem3D_OBJECTS += $(vtk_specfem3D_OBJECTS)
else
specfem3D_OBJECTS += $(vtk_specfem3D_STUBS)
endif


#######################################

####
#### rules for executables
####

ifeq ($(CUDA),yes)
## cuda version

ifeq ($(CUDA5),yes)

## cuda 5 version
${E}/xspecfem3D: $(specfem3D_OBJECTS) $(specfem3D_SHARED_OBJECTS)
	@echo ""
	@echo "building xspecfem3D with CUDA 5 support"
	@echo ""
	${FCLINK} -o $@ $+ $(LDFLAGS) $(MPILIBS) $(CUDA_LINK) $(LIBS)
	@echo ""

else

## cuda 4 version
${E}/xspecfem3D: $(specfem3D_OBJECTS) $(specfem3D_SHARED_OBJECTS)
	@echo ""
	@echo "building xspecfem3D with CUDA 4 support"
	@echo ""
	${FCLINK} -o $@ $+ $(LDFLAGS) $(MPILIBS) $(CUDA_LINK) $(LIBS)
	@echo ""

endif

else

## non-cuda version
${E}/xspecfem3D: $(specfem3D_OBJECTS) $(specfem3D_SHARED_OBJECTS)
	@echo ""
	@echo "building xspecfem3D without CUDA support"
	@echo ""
## use MPI here
## DK DK add OpenMP compiler flag here if needed
#	${MPIFCCOMPILE_CHECK} -qsmp=omp -o $@ $+ $(LDFLAGS) $(MPILIBS) $(LIBS)
	${MPIFCCOMPILE_CHECK} -o $@ $+ $(LDFLAGS) $(MPILIBS) $(LIBS)
	@echo ""

endif

#######################################

## compilation directories
S := ${S_TOP}/src/specfem3D
$(specfem3D_OBJECTS): S = ${S_TOP}/src/specfem3D

####
#### rule for each .o file below
####

###
### additional dependencies
###

$O/write_output_ASDF.solverstatic_adios.o: $O/asdf_helpers.shared_adios.o $O/asdf_helpers_writers.shared_adios.o $O/asdf_data.solverstatic_module.o

$O/write_seismograms.solverstatic.o: $O/asdf_data.solverstatic_module.o

$O/compute_arrays_source.solverstatic.o: $O/write_seismograms.solverstatic.o
$O/iterate_time.solverstatic.o: $O/write_seismograms.solverstatic.o
$O/iterate_time_undoatt.solverstatic.o: $O/write_seismograms.solverstatic.o
$O/locate_receivers.solverstatic.o: $O/write_seismograms.solverstatic.o
$O/read_adjoint_sources.solverstatic.o: $O/write_seismograms.solverstatic.o
$O/comp_source_time_function.solver.o: $O/ucb_heaviside.solverstatic_module.o

###
### specfem3D - optimized flags and dependence on values from mesher here
###
$O/%.solverstatic_module.o: $S/%.F90 ${OUTPUT}/values_from_mesher.h $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.solverstatic_module.o: $S/%.f90 ${OUTPUT}/values_from_mesher.h $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.solverstatic.o: $S/%.f90 ${OUTPUT}/values_from_mesher.h $O/shared_par.shared_module.o $O/specfem3D_par.solverstatic_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.solverstatic.o: $S/%.F90 ${OUTPUT}/values_from_mesher.h $O/shared_par.shared_module.o $O/specfem3D_par.solverstatic_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.solverstatic_openmp.o: $S/%.f90 ${OUTPUT}/values_from_mesher.h $O/shared_par.shared_module.o $O/specfem3D_par.solverstatic_module.o
## DK DK add OpenMP compiler flag here if needed
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -qsmp=omp -o $@ $<


$O/%.solverstatic_adios.o: $S/%.f90 ${OUTPUT}/values_from_mesher.h $O/shared_par.shared_module.o $O/specfem3D_par.solverstatic_module.o $O/adios_helpers.shared_adios.o
	${MPIFCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.solverstatic_adios.o: $S/%.F90 ${OUTPUT}/values_from_mesher.h $O/shared_par.shared_module.o $O/specfem3D_par.solverstatic_module.o $O/adios_helpers.shared_adios.o
	${MPIFCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

###
### no dependence on values from mesher here
###

$O/%.solver.o: $S/%.f90 $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.solver.o: $S/%.F90 $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.cc.o: $S/%.c ${SETUP}/config.h
	${CC} -c $(CPPFLAGS) $(CFLAGS) -o $@ $< 

###
### CUDA 5 only
###

$(cuda_specfem3D_DEVICE_OBJ): $(cuda_OBJECTS)
	${NVCCLINK} -o $(cuda_specfem3D_DEVICE_OBJ) $(cuda_OBJECTS)

###
### VTK compilation
###

$O/%.visualcc.o: $S/%.cpp ${SETUP}/config.h
	${CC} -c $(CPPFLAGS) $(MPI_INCLUDES) -o $@ $<

$O/%.visualc.o: $S/%.c ${SETUP}/config.h
	${CC} -c $(CPPFLAGS) $(MPI_INCLUDES) -o $@ $<

