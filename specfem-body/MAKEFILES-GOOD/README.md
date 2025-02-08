This is a SPECFEM version capable of reading and interpreting SEMUCB-type models. 

I/ Compilation
Makefile.intel is made for intel compilers and intel version of mpi (impi). GCC can be required as well as there are a few C routines.
Makefile.gnu is for GNU compilers and openmpi
Once the proper modules are loaded, rename the chosen makefile to Makefile and just run `make` and it *should* compile smoothly.

II/ Preparation of the simulation
Everything you will need to set up is in the DATA/ directory. You should be able to see the Par_file that governs the characteristics of your simulation. Refer to the onlone manual to decide how to set up NEX_XI and NEX_ETA, as those will determine the minimal wavelength of your simulation and the number of processors you can use for the run.

In this DATA/ directory, you can set up the source (CMTSOLUTION file) and the stations (STATIONS file) and the background Earth model in the berkeley_model subdirectory.

Finally, you will find the source code in the src/ directory.

!!! Make sure you compile the code after each time you change the Par_file !!!

III/ Run
To run specfem, you must run the mesher first and then the solver. You can see an exemple of a slurm file in this current directory (sub_specfem.job).

For any further question or enquiry, you can contact me at sevan.adourian@berkeley.edu (en francais, ca marche aussi!)
