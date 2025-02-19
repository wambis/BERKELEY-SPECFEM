!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  6 . 0
!          --------------------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

! setup/constants.h.  Generated from constants.h.in by configure.

!
!--- user can modify parameters below
!

!
! solver in single or double precision depending on the machine (4 or 8 bytes)
!
!  ALSO CHANGE FILE precision.h ACCORDINGLY
!
  integer, parameter :: SIZE_REAL = 4, SIZE_DOUBLE = 8

! usually the size of integer and logical variables is the same as regular single-precision real variable
  integer, parameter :: SIZE_INTEGER = SIZE_REAL
  integer, parameter :: SIZE_LOGICAL = SIZE_REAL

! set to SIZE_REAL to run in single precision
! set to SIZE_DOUBLE to run in double precision (increases memory size by 2)
  integer, parameter :: CUSTOM_REAL = SIZE_REAL

! ********************************************************************************************************   
logical, parameter :: SAVE_MIRRORS = .false. 
    
!*********************************************************************************************************


! if files on a local path on each node are also seen as global with same path
! set to .true. typically on a machine with a common (shared) file system, e.g. LUSTRE, GPFS or NFS-mounted /home
! set to .false. typically on a cluster of nodes with local disks used to store the mesh (not very common these days)
! if running on a cluster of nodes with local disks, also customize global path
! to local files in create_serial_name_database.f90 ("20 format ...")
! Flag is used only when one checks the mesh with the serial codes
! ("xcheck_buffers_1D" etc.), ignore it if you do not plan to use them
  logical, parameter :: LOCAL_PATH_IS_ALSO_GLOBAL = .true.

! input, output and main MPI I/O files
! note: careful with these unit numbers, we mostly use units in the 40-50 range.
!       cray fortran e.g. reserves 0,5,6 (standard error,input,output units) and 100-102
  integer, parameter :: ISTANDARD_OUTPUT = 6
  integer, parameter :: IIN = 40,IOUT = 41
! uncomment this to write messages to a text file
  integer, parameter :: IMAIN = 42
! uncomment this to write messages to the screen (slows down the code)
! integer, parameter :: IMAIN = ISTANDARD_OUTPUT
! I/O unit for noise data files
  integer, parameter :: IIN_NOISE = 43,IOUT_NOISE = 44
! I/O unit for adjoint source files
  integer, parameter :: IIN_ADJ = 45
! local file unit for output of buffers
  integer, parameter :: IOUT_BUFFERS = 46
! I/O unit for source and receiver vtk file
  integer, parameter :: IOUT_VTK = 47
! I/O unit for sac files
  integer, parameter :: IOUT_SAC = 48

! I/O Bluegene serial access
  logical, parameter :: IO_BLUEGENE_SERIAL = .false.

! R_EARTH is the radius of the bottom of the oceans (radius of Earth in m)
  double precision, parameter :: R_EARTH = 6371000.d0
! uncomment line below for PREM with oceans
! double precision, parameter :: R_EARTH = 6368000.d0

! radius of the Earth in km
  double precision, parameter :: R_EARTH_KM = R_EARTH / 1000.d0

! average density in the full Earth to normalize equation
  double precision, parameter :: RHOAV = 5514.3d0


!!-----------------------------------------------------------
!!
!! for topography/bathymetry model
!!
!!-----------------------------------------------------------
!! (uncomment desired resolution)

!!--- ETOPO5 5-minute model, smoothed Harvard version
!! size of topography and bathymetry file
!  integer, parameter :: NX_BATHY = 4320,NY_BATHY = 2160
!! resolution of topography file in minutes
!  integer, parameter :: RESOLUTION_TOPO_FILE = 5
!! pathname of the topography file
!  character (len=*), parameter :: PATHNAME_TOPO_FILE = 'DATA/topo_bathy/topo_bathy_etopo5_smoothed_Harvard.dat'

!! DK DK for Roland_Sylvain this below could be changed if needed (i.e. using a non-smoothed version for instance)
!---  ETOPO4 4-minute model created by subsampling and smoothing etopo-2
! size of topography and bathymetry file
!  integer, parameter :: NX_BATHY = 5400,NY_BATHY = 2700
! resolution of topography file in minutes
!  integer, parameter :: RESOLUTION_TOPO_FILE = 4
! pathname of the topography file
!  character (len=*), parameter :: PATHNAME_TOPO_FILE = 'DATA/topo_bathy/topo_bathy_etopo4_smoothed_window_7.dat'

!!--- ETOPO2 2-minute model, not implemented yet
!! size of topography and bathymetry file
!  integer, parameter :: NX_BATHY = 10800,NY_BATHY = 5400
!! resolution of topography file in minutes
!  integer, parameter :: RESOLUTION_TOPO_FILE = 2
!! pathname of the topography file
!  character (len=*), parameter :: PATHNAME_TOPO_FILE = 'DATA/topo_bathy/topo_bathy_etopo2_smoothed_window7.dat'

!!--- ETOPO1 1-minute model, implemented now, but data file must be created first
!! size of topography and bathymetry file
!  integer, parameter :: NX_BATHY = 21600,NY_BATHY = 10800
!! resolution of topography file in minutes
!  integer, parameter :: RESOLUTION_TOPO_FILE = 1
!! pathname of the topography file (un-smoothed)
!  character (len=*), parameter :: PATHNAME_TOPO_FILE = 'DATA/topo_bathy/ETOPO1.xyz'


!!--- ETOPO1 5 - Berkeley filtered
!! size of topography and bathymetry file
  integer, parameter :: NX_BATHY = 360,NY_BATHY = 180
!! resolution of topography file in minutes
  integer, parameter :: RESOLUTION_TOPO_FILE = 2
!! pathname of the topography file (un-smoothed)
  character (len=*), parameter :: PATHNAME_TOPO_FILE = 'DATA/berkeley_model/ETOPO5_1x1_filtre.dat'

! For the reference ellipsoid to convert geographic latitudes to geocentric:
!
! From Dahlen and Tromp (1998): "Spherically-symmetric Earth models all have the same hydrostatic
! surface ellipticity 1/299.8. This is 0.5 percent smaller than observed flattening of best-fitting ellipsoid 1/298.3.
! The discrepancy is referred to as the "excess equatorial bulge of the Earth",
! an early discovery of artificial satellite geodesy."
!
! From Paul Melchior, IUGG General Assembly, Vienna, Austria, August 1991 Union lecture,
! available at http://www.agu.org/books/sp/v035/SP035p0047/SP035p0047.pdf :
! "It turns out that the spheroidal models constructed on the basis of the spherically-symmetric models (PREM, 1066A)
! by using the Clairaut differential equation to calculate the flattening in function of the radius vector imply hydrostaticity.
! These have surface ellipticity 1/299.8 and a corresponding dynamical flattening of .0033 (PREM).
! The actual ellipticty of the Earth for a best-fitting ellipsoid is 1/298.3 with a corresponding dynamical flattening of .0034."
!
! Thus, flattening f = 1/299.8 is what is used in SPECFEM3D_GLOBE, as it should.
! And thus eccentricity squared e^2 = 1 - (1-f)^2 = 1 - (1 - 1/299.8)^2 = 0.00665998813529,
! and the correction factor used in the code to convert geographic latitudes to geocentric
! is 1 - e^2 = (1-f)^2 = (1 - 1/299.8)^2 = 0.9933400118647.
!
! As a comparison, the classical World Geodetic System reference ellipsoid WGS 84
! (see e.g. http://en.wikipedia.org/wiki/World_Geodetic_System) has f = 1/298.2572236.
  double precision, parameter :: FLATTENING_F = 1.d0 / 299.8d0
  double precision, parameter :: ONE_MINUS_F_SQUARED = (1.d0 - FLATTENING_F)**2

! Use GLL points to capture TOPOGRAPHY and ELLIPTICITY (experimental feature, currently does not work, DO NOT USE)
  logical,parameter :: USE_GLL = .false.

! maximum depth of the oceans in trenches and height of topo in mountains
! to avoid taking into account spurious oscillations in global model ETOPO
  logical, parameter :: USE_MAXIMUM_HEIGHT_TOPO = .false.
  integer, parameter :: MAXIMUM_HEIGHT_TOPO = +20000
  logical, parameter :: USE_MAXIMUM_DEPTH_OCEANS = .false.
  integer, parameter :: MAXIMUM_DEPTH_OCEANS = -20000

! minimum thickness in meters to include the effect of the oceans and topo
  double precision, parameter :: MINIMUM_THICKNESS_3D_OCEANS = 100.d0


!!-----------------------------------------------------------
!!
!! for crustal model
!!
!!-----------------------------------------------------------
!! (uncomment desired model)

! crustal model constants
  integer, parameter :: ICRUST_CRUST1 = 1
  integer, parameter :: ICRUST_CRUST2 = 2
  integer, parameter :: ICRUST_CRUSTMAPS = 3
  integer, parameter :: ICRUST_EPCRUST = 4
  integer, parameter :: ICRUST_BERKELEY = 5

! increase smoothing for critical regions  (increases mesh stability)
  logical, parameter :: SMOOTH_CRUST = .true.

! use sedimentary layers in crustal model
  logical, parameter :: INCLUDE_SEDIMENTS_CRUST = .true.
  double precision, parameter :: MINIMUM_SEDIMENT_THICKNESS = 2.d0 ! minimim thickness in km

! default crustal model
! (used as default when CRUSTAL flag is set for simulation)
!!-- uncomment for using Crust1.0
!  integer, parameter :: ITYPE_CRUSTAL_MODEL = ICRUST_CRUST1
!!-- uncomment for using Crust2.0
!  integer, parameter :: ITYPE_CRUSTAL_MODEL = ICRUST_CRUST2
!!-- uncomment for using General Crustmaps instead
!  integer, parameter :: ITYPE_CRUSTAL_MODEL = ICRUST_CRUSTMAPS
!!-- uncomment for using EPcrust instead (European regional model)
!  integer, parameter :: ITYPE_CRUSTAL_MODEL = ICRUST_EPCRUST
!!-- uncomment for using BERKELEY homogenized crust instead (for use with SEMUCB model)
  integer, parameter :: ITYPE_CRUSTAL_MODEL = ICRUST_BERKELEY


!!-----------------------------------------------------------
!!
!! Gauss-Lobatto-Legendre resolution
!!
!!-----------------------------------------------------------

! number of GLL points in each direction of an element (degree plus one)
  integer, parameter :: NGLLX = 5
  integer, parameter :: NGLLY = NGLLX
  integer, parameter :: NGLLZ = NGLLX


!!-----------------------------------------------------------
!!
!! source/receiver setup
!!
!!-----------------------------------------------------------

! flag to exclude elements that are too far from target in source detection
  logical, parameter :: USE_DISTANCE_CRITERION = .true.

! flag to display detailed information about location of stations
  logical, parameter :: DISPLAY_DETAILS_STATIONS = .false.

! maximum length of station and network name for receivers
  integer, parameter :: MAX_LENGTH_STATION_NAME = 32
  integer, parameter :: MAX_LENGTH_NETWORK_NAME = 8

! we mimic a triangle of half duration equal to half_duration_triangle
! using a Gaussian having a very close shape, as explained in Figure 4.2
! of the manual. This source decay rate to mimic an equivalent triangle
! was found by trial and error
  double precision, parameter :: SOURCE_DECAY_MIMIC_TRIANGLE = 1.628d0

! maximum number of sources and receivers to locate simultaneously
  integer, parameter :: NSOURCES_SUBSET_MAX = 100
  integer, parameter :: NREC_SUBSET_MAX = 200

! use a force source located exactly at a grid point instead of a CMTSOLUTION source
! this can be useful e.g. for asteroid impact simulations
! in which the source is a vertical force, normal force, impact etc.
  logical, parameter :: USE_FORCE_POINT_SOURCE = .false.
  double precision, parameter :: FACTOR_FORCE_SOURCE = 1.d15
  integer, parameter :: COMPONENT_FORCE_SOURCE = 3  ! takes direction in comp E/N/Z = 1/2/3

! use this t0 as earliest starting time rather than the automatically calculated one
! (must be positive and bigger than the automatically one to be effective;
!  simulation will start at t = - t0)
  double precision, parameter :: USER_T0 = 0.0d0
! source-time function is a UCB-style filtered heaviside
  logical, parameter :: STF_IS_UCB_HEAVISIDE = .true. 
! distance threshold (in km) above which we consider that a receiver
! is located outside the mesh and therefore excluded from the station list
  double precision, parameter :: THRESHOLD_EXCLUDE_STATION = 50.d0


!!-----------------------------------------------------------
!!
!! mesh optimization
!!
!!-----------------------------------------------------------

! the first doubling is implemented right below the Moho
! it seems optimal to implement the three other doublings at these depths
! in the mantle
  double precision, parameter :: DEPTH_SECOND_DOUBLING_OPTIMAL = 1650000.d0
! in the outer core
  double precision, parameter :: DEPTH_THIRD_DOUBLING_OPTIMAL  = 3860000.d0
! in the outer core
  double precision, parameter :: DEPTH_FOURTH_DOUBLING_OPTIMAL = 5000000.d0

! Boundary Mesh -- save Moho, 400, 670 km discontinuity topology files (in
! the mesher) and use them for the computation of boundary kernel (in the solver)
  logical, parameter :: SAVE_BOUNDARY_MESH = .false.

!!-----------------------------------------------------------
!!
!! GPU optimization
!!
!!-----------------------------------------------------------
! added these parameters for the GPU version of the solver

! asynchronuous memcopy between CPU and GPU
  logical, parameter :: GPU_ASYNC_COPY = .true.

! mesh coloring
! add mesh coloring for the GPU + MPI implementation
  logical, parameter :: USE_MESH_COLORING_GPU = .false.
  integer, parameter :: MAX_NUMBER_OF_COLORS = 1000
! enhanced coloring:
! using Droux algorithm
! try several times with one more color before giving up
  logical, parameter :: USE_DROUX_OPTIMIZATION = .false.
  integer, parameter :: MAX_NB_TRIES_OF_DROUX_1993 = 15
! using balancing algorithm
! postprocess the colors to balance them if Droux (1993) algorithm is not used
  logical, parameter :: BALANCE_COLORS_SIMPLE_ALGO = .true.


!!-----------------------------------------------------------
!!
!! ADIOS Related values
!!
!!-----------------------------------------------------------

  integer, parameter :: ADIOS_BUFFER_SIZE_IN_MB = 200
  character(len=256), parameter :: ADIOS_TRANSPORT_METHOD = 'MPI'


!!-----------------------------------------------------------
!!
!! adjoint kernel outputs
!!
!!-----------------------------------------------------------

! asynchronuous reading of adjoint sources
  logical, parameter :: IO_ASYNC_COPY = .true.

! regular kernel parameters
  character(len=*), parameter :: PATHNAME_KL_REG = 'DATA/kl_reg_grid.txt'
  integer, parameter :: NM_KL_REG_LAYER = 100
  integer, parameter :: NM_KL_REG_PTS = 300000
  real, parameter :: KL_REG_MIN_LON = 0.0
  real, parameter :: KL_REG_MAX_LON = 360.0
  real, parameter :: KL_REG_MIN_LAT = -90.0
  real, parameter :: KL_REG_MAX_LAT = +90.0


!!-----------------------------------------------------------
!!
!! new version compatibility
!!
!!-----------------------------------------------------------
! old version 5.1.5 uses full 3d attenuation arrays (set to .true.), custom_real for attenuation arrays (Qmu_store, tau_e_store)
! new version uses optional full 3d attenuation
  logical, parameter :: USE_OLD_VERSION_5_1_5_FORMAT = .false.


!!-----------------------------------------------------------
!!
!! time stamp information
!!
!!-----------------------------------------------------------

! print date and time estimate of end of run in another country,
! in addition to local time.
! For instance: the code runs at Caltech in California but the person
! running the code is connected remotely from France, which has 9 hours more.
! The time difference with that remote location can be positive or negative
  logical, parameter :: ADD_TIME_ESTIMATE_ELSEWHERE = .false.
  integer, parameter :: HOURS_TIME_DIFFERENCE = +9
  integer, parameter :: MINUTES_TIME_DIFFERENCE = +0


!!-----------------------------------------------------------
!!
!! movie outputs
!!
!!-----------------------------------------------------------

! runs external bash script after storing movie files
! (useful for postprocessing during simulation time)
  logical, parameter :: RUN_EXTERNAL_MOVIE_SCRIPT = .false.
  character(len=256) :: MOVIE_SCRIPT_NAME = "./tar_movie_files.sh"


!!-----------------------------------------------------------
!!
!! debugging flags
!!
!!-----------------------------------------------------------

! flags to actually assemble with MPI or not
! and to actually match fluid and solid regions of the Earth or not
! should always be set to true except when debugging code
  logical, parameter :: ACTUALLY_ASSEMBLE_MPI_SLICES = .true.
  logical, parameter :: ACTUALLY_ASSEMBLE_MPI_CHUNKS = .true.
  logical, parameter :: ACTUALLY_COUPLE_FLUID_CMB = .true.
  logical, parameter :: ACTUALLY_COUPLE_FLUID_ICB = .true.

! flag to turn off the conversion of geographic to geocentric coordinates for
! the seismic source and the stations; i.e. assume a perfect sphere, which
! can be useful for benchmarks of a spherical Earth with fictitious sources and stations
  logical, parameter :: ASSUME_PERFECT_SPHERE = .false.

! flags to do benchmark runs to measure scaling of the code
! for a limited number of time steps only, setting the initial field to 1
! to make sure gradual underflow trapping does not slow down the code
  logical, parameter :: DO_BENCHMARK_RUN_ONLY = .false.
  integer, parameter :: NSTEP_FOR_BENCHMARK = 300
  logical, parameter :: SET_INITIAL_FIELD_TO_1_IN_BENCH = .true.

!!-----------------------------------------------------------
!!
!! for Roland_Sylvain integrals
!!
!!-----------------------------------------------------------

!! DK DK for Roland_Sylvain
  logical, parameter :: ROLAND_SYLVAIN = .false.

! altitude of the observation point
  double precision, parameter :: altitude_of_observation_point = 255.d0 ! in km

! standard gravity at the surface of the Earth
  double precision, parameter :: STANDARD_GRAVITY_EARTH = 9.80665d0 ! in m.s-2

  double precision, parameter :: SQUARE_ROOT_OF_THREE = 1.73205080756887729352d0

! position of the observation point in non-dimensionalized value

!! DK DK along X only
! double precision, parameter :: x_observation = (R_EARTH_KM + altitude_of_observation_point) / R_EARTH_KM
! double precision, parameter :: y_observation = 0.d0
! double precision, parameter :: z_observation = 0.d0

!! DK DK along diagonal
! double precision, parameter :: x_observation = (R_EARTH_KM + altitude_of_observation_point) / R_EARTH_KM / SQUARE_ROOT_OF_THREE
! double precision, parameter :: y_observation = x_observation
! double precision, parameter :: z_observation = x_observation

!! DK DK point randomly chosen in space, not at the right elevation we want
  double precision, parameter :: x_observation = (R_EARTH_KM + altitude_of_observation_point) / R_EARTH_KM / SQUARE_ROOT_OF_THREE
  double precision, parameter :: y_observation = x_observation * 1.12d0
  double precision, parameter :: z_observation = x_observation * 1.38d0

!------------------------------------------------------
!----------- do not modify anything below -------------
!------------------------------------------------------

! on some processors (e.g. some Intel chips) it is necessary to suppress underflows
! by using a small initial field instead of zero
  logical, parameter :: FIX_UNDERFLOW_PROBLEM = .true.

! some useful constants
  double precision, parameter :: PI = 3.141592653589793d0
  double precision, parameter :: TWO_PI = 2.d0 * PI
  double precision, parameter :: PI_OVER_FOUR = PI / 4.d0
  double precision, parameter :: PI_OVER_TWO = PI / 2.0d0

! to convert angles from degrees to radians
  double precision, parameter :: DEGREES_TO_RADIANS = PI / 180.d0
! to convert angles from radians to degrees
  double precision, parameter :: RADIANS_TO_DEGREES = 180.d0 / PI

! 3-D simulation
  integer, parameter :: NDIM = 3

! dimension of the boundaries of the slices
  integer, parameter :: NDIM2D = 2

! number of nodes for 2D and 3D shape functions for hexahedra with 27 nodes
  integer, parameter :: NGNOD = 27, NGNOD2D = 9

! Deville routines optimized for NGLLX = NGLLY = NGLLZ = 5
  integer, parameter :: m1 = NGLLX, m2 = NGLLX * NGLLY
  integer, parameter :: NGLLCUBE = NGLLX * NGLLY * NGLLZ

! mid-points inside a GLL element
  integer, parameter :: MIDX = (NGLLX+1)/2
  integer, parameter :: MIDY = (NGLLY+1)/2
  integer, parameter :: MIDZ = (NGLLZ+1)/2

! gravitational constant in S.I. units i.e. in m3 kg-1 s-2, or equivalently in N.(m/kg)^2
!! DK DK April 2014: switched to the 2010 Committee on Data for Science and Technology (CODATA) recommended value
!! DK DK see e.g. http://www.physics.nist.gov/cgi-bin/cuu/Value?bg
!! DK DK and http://en.wikipedia.org/wiki/Gravitational_constant
! double precision, parameter :: GRAV = 6.6723d-11
  double precision, parameter :: GRAV = 6.67384d-11

! a few useful constants
  double precision, parameter :: ZERO = 0.d0,ONE = 1.d0,TWO = 2.d0,HALF = 0.5d0
  double precision, parameter :: ONE_HALF = HALF, ONE_FOURTH = 0.25d0, ONE_EIGHTH = 0.125d0, ONE_SIXTEENTH = 0.0625d0

  real(kind=CUSTOM_REAL), parameter :: &
    ONE_THIRD   = 1._CUSTOM_REAL/3._CUSTOM_REAL, &
    TWO_THIRDS  = 2._CUSTOM_REAL/3._CUSTOM_REAL, &
    FOUR_THIRDS = 4._CUSTOM_REAL/3._CUSTOM_REAL

! very large and very small values
  double precision, parameter :: HUGEVAL = 1.d+30,TINYVAL = 1.d-9

! very large real value declared independently of the machine
  real(kind=CUSTOM_REAL), parameter :: HUGEVAL_SNGL = 1.e+30_CUSTOM_REAL

! very large integer value
  integer, parameter :: HUGEINT = 100000000

! normalized radius of free surface
  double precision, parameter :: R_UNIT_SPHERE = ONE

! fixed thickness of 3 km for PREM oceans
  double precision, parameter :: THICKNESS_OCEANS_PREM = 3000.d0 / R_EARTH

! shortest radius at which crust is implemented (80 km depth)
! to be constistent with the D80 discontinuity, we impose the crust only above it
  double precision, parameter :: R_DEEPEST_CRUST = (R_EARTH - 80000.d0) / R_EARTH

! definition of an Eotvos compared to S.I. units.
! The unit of gravity gradient is the Eotvos, which is equivalent to 1e-9 s-2 (or 1e-4 mGal/m).
! A person walking at a distance of 2 meters provides a gravity gradient signal of approximately one Eotvos.
! Mountains can create signals of several hundred Eotvos.
  double precision, parameter :: SI_UNITS_TO_EOTVOS = 1.d+9

! maximum number of chunks (full sphere)
  integer, parameter :: NCHUNKS_MAX = 6

! define block type based upon chunk number (between 1 and 6)
! do not change this numbering, chunk AB must be number 1 for central cube
  integer, parameter :: CHUNK_AB = 1
  integer, parameter :: CHUNK_AC = 2
  integer, parameter :: CHUNK_BC = 3
  integer, parameter :: CHUNK_AC_ANTIPODE = 4
  integer, parameter :: CHUNK_BC_ANTIPODE = 5
  integer, parameter :: CHUNK_AB_ANTIPODE = 6

! maximum number of regions in the mesh
  integer, parameter :: MAX_NUM_REGIONS = 3

! define flag for regions of the global Earth mesh
  integer, parameter :: IREGION_CRUST_MANTLE = 1
  integer, parameter :: IREGION_OUTER_CORE = 2
  integer, parameter :: IREGION_INNER_CORE = 3

! define flag for elements
  integer, parameter :: IFLAG_CRUST = 1

  integer, parameter :: IFLAG_80_MOHO = 2
  integer, parameter :: IFLAG_220_80 = 3
  integer, parameter :: IFLAG_670_220 = 4
  integer, parameter :: IFLAG_MANTLE_NORMAL = 5

  integer, parameter :: IFLAG_OUTER_CORE_NORMAL = 6

  integer, parameter :: IFLAG_INNER_CORE_NORMAL = 7
  integer, parameter :: IFLAG_MIDDLE_CENTRAL_CUBE = 8
  integer, parameter :: IFLAG_BOTTOM_CENTRAL_CUBE = 9
  integer, parameter :: IFLAG_TOP_CENTRAL_CUBE = 10
  integer, parameter :: IFLAG_IN_FICTITIOUS_CUBE = 11

  integer, parameter :: NSPEC2D_XI_SUPERBRICK = 8
  integer, parameter :: NSPEC2D_ETA_SUPERBRICK = 8
  integer, parameter :: NSPEC2D_XI_SUPERBRICK_1L = 6
  integer, parameter :: NSPEC2D_ETA_SUPERBRICK_1L = 6

! dummy flag used for mesh display purposes only
  integer, parameter :: IFLAG_DUMMY = 100

! max number of layers that are used in the radial direction to build the full mesh
  integer, parameter :: MAX_NUMBER_OF_MESH_LAYERS = 15

! define number of spectral elements and points in basic symmetric mesh doubling superbrick
  integer, parameter :: NSPEC_DOUBLING_SUPERBRICK = 32
  integer, parameter :: NGLOB_DOUBLING_SUPERBRICK = 67
  integer, parameter :: NSPEC_SUPERBRICK_1L = 28
  integer, parameter :: NGLOB_SUPERBRICK_1L = 58
  integer, parameter :: NGNOD_EIGHT_CORNERS = 8

! define flag for reference 1D Earth model
  integer, parameter :: REFERENCE_MODEL_PREM   = 1
  integer, parameter :: REFERENCE_MODEL_IASP91 = 2
  integer, parameter :: REFERENCE_MODEL_1066A  = 3
  integer, parameter :: REFERENCE_MODEL_AK135F_NO_MUD = 4
  integer, parameter :: REFERENCE_MODEL_1DREF = 5
  integer, parameter :: REFERENCE_MODEL_JP1D  = 6
  integer, parameter :: REFERENCE_MODEL_SEA1D = 7
  integer, parameter :: REFERENCE_MODEL_BERKELEY = 8

! define flag for 3D Earth model
  integer, parameter :: THREE_D_MODEL_S20RTS   = 1
  integer, parameter :: THREE_D_MODEL_S362ANI   = 2
  integer, parameter :: THREE_D_MODEL_S362WMANI = 3
  integer, parameter :: THREE_D_MODEL_S362ANI_PREM  = 4
  integer, parameter :: THREE_D_MODEL_S29EA  = 5
  integer, parameter :: THREE_D_MODEL_SEA99_JP3D  = 6
  integer, parameter :: THREE_D_MODEL_SEA99  = 7
  integer, parameter :: THREE_D_MODEL_JP3D  = 8
  integer, parameter :: THREE_D_MODEL_PPM  = 9     ! format for point profile models
  integer, parameter :: THREE_D_MODEL_GLL  = 10    ! format for iterations with GLL mesh
  integer, parameter :: THREE_D_MODEL_S40RTS = 11
  integer, parameter :: THREE_D_MODEL_GAPP2  = 12
  integer, parameter :: THREE_D_MODEL_BERKELEY  = 13

! number of standard linear solids for attenuation
  integer, parameter :: N_SLS = 3

! computation of standard linear solids in meshfem3D
! ATTENUATION_COMP_RESOLUTION: Number of Digits after decimal
! ATTENUATION_COMP_MAXIMUM:    Maximum Q Value
  integer, parameter :: ATTENUATION_COMP_RESOLUTION = 1
  integer, parameter :: ATTENUATION_COMP_MAXIMUM    = 5000

! for determination of the attenuation period range
! if this is set to .true. then the hardcoded values will be used
! otherwise they are computed automatically from the Number of elements
! This *may* be a useful parameter for Benchmarking against older versions
  logical, parameter :: ATTENUATION_RANGE_PREDEFINED = .false.

! flag for the four edges of each slice and for the bottom edge
  integer, parameter :: XI_MIN  = 1
  integer, parameter :: XI_MAX  = 2
  integer, parameter :: ETA_MIN = 3
  integer, parameter :: ETA_MAX = 4
  integer, parameter :: BOTTOM  = 5

! flags to select the right corner in each slice
  integer, parameter :: ILOWERLOWER = 1
  integer, parameter :: ILOWERUPPER = 2
  integer, parameter :: IUPPERLOWER = 3
  integer, parameter :: IUPPERUPPER = 4

! number of points in each AVS or OpenDX quadrangular cell for movies
  integer, parameter :: NGNOD2D_AVS_DX = 4

! number of faces a given slice can share with other slices
! this is at most 2, except when there is only once slice per chunk
! in which case it is 4
  integer, parameter :: NUMFACES_SHARED = 4

! number of corners a given slice can share with other slices
! this is at most 1, except when there is only once slice per chunk
! in which case it is 4
  integer, parameter :: NUMCORNERS_SHARED = 4

! number of slaves per corner
  integer, parameter :: NUMSLAVES = 2

! number of layers in PREM
  integer, parameter :: NR = 640

! smallest real number on many machines =  1.1754944E-38
! largest real number on many machines =  3.4028235E+38
! small negligible initial value to avoid very slow underflow trapping
! but not too small to avoid trapping on velocity and acceleration in Newmark
  real(kind=CUSTOM_REAL), parameter :: VERYSMALLVAL = 1.E-24_CUSTOM_REAL

! displacement threshold above which we consider that the code became unstable
  real(kind=CUSTOM_REAL), parameter :: STABILITY_THRESHOLD = 1.E+25_CUSTOM_REAL

! geometrical tolerance for boundary detection
  double precision, parameter :: SMALLVAL = 0.00001d0

! small tolerance for conversion from x y z to r theta phi
  double precision, parameter :: SMALL_VAL_ANGLE = 1.d-10

! geometry tolerance parameter to calculate number of independent grid points
! sensitive to actual size of model, assumes reference sphere of radius 1
! this is an absolute value for normalized coordinates in the Earth
  double precision, parameter :: SMALLVALTOL = 1.d-10

! do not use tags for MPI messages, use dummy tag instead
  integer, parameter :: itag = 0,itag2 = 0

! for the Gauss-Lobatto-Legendre points and weights
  double precision, parameter :: GAUSSALPHA = 0.d0,GAUSSBETA = 0.d0

! number of lines per source in CMTSOLUTION file
  integer, parameter :: NLINES_PER_CMTSOLUTION_SOURCE = 13

! number of iterations to solve the non linear system for xi and eta
  integer, parameter :: NUM_ITER = 4

! number of hours per day for rotation rate of the Earth
  double precision, parameter :: HOURS_PER_DAY = 24.d0
  double precision, parameter :: SECONDS_PER_HOUR = 3600.d0

! for lookup table for gravity every 100 m in radial direction of Earth model
  integer, parameter :: NRAD_GRAVITY = 70000

!!-----------------------------------------------------------
!!
!! model constants
!!
!!-----------------------------------------------------------

! The meaningful range of Zhao et al. (1994) model is as follows:
!        latitude : 32 - 45 N
!        longitude: 130-145 E
!        depth    : 0  - 500 km
! The deepest Moho beneath Japan is 40 km
  double precision,parameter :: LAT_MAX = 45.d0
  double precision,parameter :: LAT_MIN = 32.d0
  double precision,parameter :: LON_MAX = 145.d0
  double precision,parameter :: LON_MIN = 130.d0
  double precision,parameter :: DEP_MAX = 500.d0

! parameters for GLL model (used for iterative model inversions)
  character(len=*), parameter :: PATHNAME_GLL_modeldir = 'DATA/GLL/'

  integer, parameter :: GLL_REFERENCE_1D_MODEL = REFERENCE_MODEL_1DREF
  integer, parameter :: GLL_REFERENCE_MODEL = THREE_D_MODEL_S29EA

!!-- uncomment for using S362ANI as reference model  (Hejun Zhu)
!  integer, parameter :: GLL_REFERENCE_1D_MODEL = REFERENCE_MODEL_1DREF
!  integer, parameter :: GLL_REFERENCE_MODEL = THREE_D_MODEL_S362ANI

!!-----------------------------------------------------------
!!
!! crustal stretching
!!
!!-----------------------------------------------------------

! for the stretching of crustal elements in the case of 3D models
! (values are chosen for 3D models to have RMOHO_FICTICIOUS at 35 km
!  and RMIDDLE_CRUST to become 15 km with stretching function stretch_tab)
  double precision, parameter :: MAX_RATIO_CRUST_STRETCHING = 0.75d0
  double precision, parameter :: RMOHO_STRETCH_ADJUSTEMENT = 5000.d0 ! moho up to 35km
  double precision, parameter :: R80_STRETCH_ADJUSTEMENT = -40000.d0 ! r80 down to 120km

! adapted regional moho stretching
! 1 chunk simulations, 3-layer crust
  logical, parameter :: REGIONAL_MOHO_MESH = .false.
  logical, parameter :: REGIONAL_MOHO_MESH_EUROPE = .false. ! used only for fixing time step
  logical, parameter :: REGIONAL_MOHO_MESH_ASIA = .false.   ! used only for fixing time step
  logical, parameter :: HONOR_DEEP_MOHO = .false.
!!-- uncomment for e.g. Europe case, where deep moho is rare
!  double precision, parameter :: RMOHO_STRETCH_ADJUSTEMENT = -15000.d0  ! moho mesh boundary down to 55km
!!-- uncomment for deep moho cases, e.g. Asia case (Himalayan moho)
!  double precision, parameter :: RMOHO_STRETCH_ADJUSTEMENT = -20000.d0  ! moho mesh boundary down to 60km

!!-----------------------------------------------------------
!!
!! mesh tweaking
!!
!!-----------------------------------------------------------

! to suppress the crustal layers
! (replaced by an extension of the mantle: R_EARTH is not modified, but no more crustal doubling)
  logical, parameter :: SUPPRESS_CRUSTAL_MESH = .false.

! to inflate the central cube (set to 0.d0 for a non-inflated cube)
  double precision, parameter :: CENTRAL_CUBE_INFLATE_FACTOR = 0.41d0

! to add a fourth doubling at the bottom of the outer core
  logical, parameter :: ADD_4TH_DOUBLING = .false.

! parameters to cut the doubling brick

! this to cut the superbrick: 3 possibilities, 4 cases max / possibility
! three possibilities: (cut in xi and eta) or (cut in xi) or (cut in eta)
! case 1: (ximin and etamin) or ximin or etamin
! case 2: (ximin and etamax) or ximax or etamax
! case 3: ximax and etamin
! case 4: ximax and etamax
  integer, parameter :: NB_CUT_CASE = 4

! corner 1: ximin and etamin
! corner 2: ximax and etamin
! corner 3: ximax and etamax
! corner 4: ximin and etamax
  integer, parameter :: NB_SQUARE_CORNERS = 4

! two possibilities: xi or eta
! face 1: ximin or etamin
! face 2: ximax or etamax
  integer, parameter :: NB_SQUARE_EDGES_ONEDIR = 2

! this for the geometry of the basic doubling brick
  integer, parameter :: NSPEC_DOUBLING_BASICBRICK = 8
  integer, parameter :: NGLOB_DOUBLING_BASICBRICK = 27

!!-----------------------------------------------------------
!!
!! for LDDRK high-order time scheme
!!
!!-----------------------------------------------------------
  integer, parameter :: NSTAGE = 6

  real(kind=CUSTOM_REAL), parameter, dimension(NSTAGE) :: ALPHA_LDDRK = &
    (/0.0_CUSTOM_REAL,-0.737101392796_CUSTOM_REAL, -1.634740794341_CUSTOM_REAL,&
      -0.744739003780_CUSTOM_REAL,-1.469897351522_CUSTOM_REAL,-2.813971388035_CUSTOM_REAL/)

  real(kind=CUSTOM_REAL), parameter, dimension(NSTAGE) :: BETA_LDDRK = &
    (/0.032918605146_CUSTOM_REAL,0.823256998200_CUSTOM_REAL,0.381530948900_CUSTOM_REAL,&
      0.200092213184_CUSTOM_REAL,1.718581042715_CUSTOM_REAL,0.27_CUSTOM_REAL/)

  real(kind=CUSTOM_REAL), parameter, dimension(NSTAGE) :: C_LDDRK = &
    (/0.0_CUSTOM_REAL,0.032918605146_CUSTOM_REAL,0.249351723343_CUSTOM_REAL,&
      0.466911705055_CUSTOM_REAL,0.582030414044_CUSTOM_REAL,0.847252983783_CUSTOM_REAL/)

