 
 !
 ! this is the parameter file for static compilation of the solver
 !
 ! mesh statistics:
 ! ---------------
 !
 !
 ! number of chunks =            6
 !
 ! these statistics include the central cube
 !
 ! number of processors =          384
 !
 ! maximum number of points per region =       682293
 !
 ! on NEC SX, make sure "loopcnt=" parameter
 ! in Makefile is greater than max vector length =      2046879
 !
 ! total elements per slice =        11760
 ! total points per slice =       789875
 !
 ! total for full 6-chunk mesh:
 ! ---------------------------
 !
 ! exact total number of spectral elements in entire mesh = 
 !    4352000.00000000     
 ! approximate total number of points in entire mesh = 
 !    292578555.000000     
 ! approximate total number of degrees of freedom in entire mesh = 
 !    826407921.000000     
 !
 ! resolution of the mesh at the surface:
 ! -------------------------------------
 !
 ! spectral elements along a great circle =         1024
 ! GLL points along a great circle =         4096
 ! average distance between points in degrees =   8.7890625E-02
 ! average distance between points in km =    9.772991    
 ! average size of a spectral element in km =    39.09196    
 !
 ! number of time steps =        22200
 !
 ! number of seismic sources =            1
 !
 
 ! approximate static memory needed by the solver:
 ! ----------------------------------------------
 !
 ! (lower bound, usually the real amount used is 5% to 10% higher)
 !
 ! (you can get a more precise estimate of the size used per MPI process
 !  by typing "size -d bin/xspecfem3D"
 !  after compiling the code with the DATA/Par_file you plan to use)
 !
 ! size of static arrays per slice =    134.058684000000       MB
 !                                 =    127.848323822021       MiB
 !                                 =   0.134058684000000       GB
 !                                 =   0.124851878732443       GiB
 !
 ! (should be below to 80% or 90% of the memory installed per core)
 ! (if significantly more, the job will not run by lack of memory )
 ! (note that if significantly less, you waste a significant amount
 !  of memory per processor core)
 ! (but that can be perfectly acceptable if you can afford it and
 !  want faster results by using more cores)
 !
 ! size of static arrays for all slices =    51.4785346560000       GB
 !                                      =    47.9431214332581       GiB
 !                                      =   5.147853465600000E-002  TB
 !                                      =   4.681945452466607E-002  TiB
 !
 
 integer, parameter :: NEX_XI_VAL =          256
 integer, parameter :: NEX_ETA_VAL =          256
 
 integer, parameter :: NSPEC_CRUST_MANTLE =        10240
 integer, parameter :: NSPEC_OUTER_CORE =          960
 integer, parameter :: NSPEC_INNER_CORE =          560
 
 integer, parameter :: NGLOB_CRUST_MANTLE =       682293
 integer, parameter :: NGLOB_OUTER_CORE =        66833
 integer, parameter :: NGLOB_INNER_CORE =        40749
 
 integer, parameter :: NSPECMAX_ANISO_IC =            1
 
 integer, parameter :: NSPECMAX_ISO_MANTLE =        10240
 integer, parameter :: NSPECMAX_TISO_MANTLE =        10240
 integer, parameter :: NSPECMAX_ANISO_MANTLE =            1
 
 integer, parameter :: NSPEC_CRUST_MANTLE_ATTENUATION =            1
 integer, parameter :: NSPEC_INNER_CORE_ATTENUATION =            1
 
 integer, parameter :: NSPEC_CRUST_MANTLE_STR_OR_ATT =            1
 integer, parameter :: NSPEC_INNER_CORE_STR_OR_ATT =            1
 
 integer, parameter :: NSPEC_CRUST_MANTLE_STR_AND_ATT =            1
 integer, parameter :: NSPEC_INNER_CORE_STR_AND_ATT =            1
 
 integer, parameter :: NSPEC_CRUST_MANTLE_STRAIN_ONLY =            1
 integer, parameter :: NSPEC_INNER_CORE_STRAIN_ONLY =            1
 
 integer, parameter :: NSPEC_CRUST_MANTLE_ADJOINT =            1
 integer, parameter :: NSPEC_OUTER_CORE_ADJOINT =            1
 integer, parameter :: NSPEC_INNER_CORE_ADJOINT =            1
 integer, parameter :: NGLOB_CRUST_MANTLE_ADJOINT =            1
 integer, parameter :: NGLOB_OUTER_CORE_ADJOINT =            1
 integer, parameter :: NGLOB_INNER_CORE_ADJOINT =            1
 integer, parameter :: NSPEC_OUTER_CORE_ROT_ADJOINT =            1
 
 integer, parameter :: NSPEC_CRUST_MANTLE_STACEY =            1
 integer, parameter :: NSPEC_OUTER_CORE_STACEY =            1
 
 integer, parameter :: NGLOB_CRUST_MANTLE_OCEANS =            1
 
 logical, parameter :: TRANSVERSE_ISOTROPY_VAL = .true.
 
 logical, parameter :: ANISOTROPIC_3D_MANTLE_VAL = .false.
 
 logical, parameter :: ANISOTROPIC_INNER_CORE_VAL = .false.
 
 logical, parameter :: ATTENUATION_VAL = .false.
 
 logical, parameter :: ATTENUATION_3D_VAL = .false.
 
 logical, parameter :: ELLIPTICITY_VAL = .false.
 
 logical, parameter :: GRAVITY_VAL = .true.
 
 logical, parameter :: OCEANS_VAL = .false.
 
 integer, parameter :: NX_BATHY_VAL = 1
 integer, parameter :: NY_BATHY_VAL = 1
 
 logical, parameter :: ROTATION_VAL = .false.
 integer, parameter :: NSPEC_OUTER_CORE_ROTATION =            1
 
 logical, parameter :: PARTIAL_PHYS_DISPERSION_ONLY_VAL = .false.
 
 integer, parameter :: NPROC_XI_VAL =            8
 integer, parameter :: NPROC_ETA_VAL =            8
 integer, parameter :: NCHUNKS_VAL =            6
 integer, parameter :: NPROCTOT_VAL =          384
 
 integer, parameter :: ATT1_VAL =            1
 integer, parameter :: ATT2_VAL =            1
 integer, parameter :: ATT3_VAL =            1
 integer, parameter :: ATT4_VAL =            1
 integer, parameter :: ATT5_VAL =            1
 
 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_CM =          560
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_CM =          560
 integer, parameter :: NSPEC2D_BOTTOM_CM =           64
 integer, parameter :: NSPEC2D_TOP_CM =         1024
 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_IC =          140
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_IC =          140
 integer, parameter :: NSPEC2D_BOTTOM_IC =           16
 integer, parameter :: NSPEC2D_TOP_IC =           16
 integer, parameter :: NSPEC2DMAX_XMIN_XMAX_OC =          144
 integer, parameter :: NSPEC2DMAX_YMIN_YMAX_OC =          144
 integer, parameter :: NSPEC2D_BOTTOM_OC =           16
 integer, parameter :: NSPEC2D_TOP_OC =           64
 integer, parameter :: NSPEC2D_MOHO =            1
 integer, parameter :: NSPEC2D_400 =            1
 integer, parameter :: NSPEC2D_670 =            1
 integer, parameter :: NSPEC2D_CMB =            1
 integer, parameter :: NSPEC2D_ICB =            1
 
 logical, parameter :: USE_DEVILLE_PRODUCTS_VAL = .true.
 integer, parameter :: NSPEC_CRUST_MANTLE_3DMOVIE = 1
 integer, parameter :: NGLOB_CRUST_MANTLE_3DMOVIE = 1
 
 integer, parameter :: NSPEC_OUTER_CORE_3DMOVIE = 1
 integer, parameter :: NM_KL_REG_PTS_VAL = 1
 
 integer, parameter :: NGLOB_XY_CM =            1
 integer, parameter :: NGLOB_XY_IC =            1
 
 logical, parameter :: ATTENUATION_1D_WITH_3D_STORAGE_VAL = .false.
 
 logical, parameter :: FORCE_VECTORIZATION_VAL = .false.
