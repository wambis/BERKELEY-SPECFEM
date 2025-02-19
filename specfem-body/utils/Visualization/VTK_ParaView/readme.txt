--------------------------------
readme
--------------------------------

ParaView, ParaView is an open-source, multi-platform application designed to visualize data sets
http://www.paraview.org/



- contains mesh2vtu/:

  mesh2vtu/ -- it requires the installation of vtk package. 
               once vtk is installed, you can set the vars 
               in the Makefile and compile all the programs

  How to collect and visualize 3D/2D data generated by the parallel simulation?

  0. run global_slice_number(2).pl to figure out the slices you need to collect
  
  1. copy database using the corresponding scripts in collect_database/
  
  2. combine all the databases using the program combine_vol_data.f90 or combine_surf_data.f90
  
  3. convert the mesh file to vtu file using mesh2vtu/
  
  4. visualizing in Paraview. After loading in the .vtu file, you can generate the source-
     receiver cross-section using the normals given by global_slice_number.pl

  Yes, it is a pretty painful process. You can certainly streamline them using a master bash
  script. To make life easier, you can use VTK scripts instead of Paraview in step 4.

  Good luck!

  Qinya, May 2007, Caltech



- example procedure to make a movie plot:

  1. run specfem with the movie options (see Par_file):

     MOVIE_SURFACE = .true.
     MOVIE_COARSE  = .false.
     
     and adjust the time steps NSTEP_BETWEEN_FRAMES   
     
     this creates binary files in directory OUTPUT_FILES/ like: moviedata000100,...

     
  2. convert binary files to AVS-files (readable by Paraview):

     in SPECFEM3D_GLOBE:  > make xcreate_movie_AVS_DX
                            
                            choose option 2 for individual files
                            
     outputs files like: OUTPUT_FILES/AVS_movie_000100.inp, ...

     
  3. make a default state-file with paraview:
   
     run > paraview &
     
     load files: under -> File / Open... select e.g. AVS_movie_009000.inp & AVS_continent_boundary.inp
     align and orient image, save state file under -> File / Save State as...
           
     save under: paraview_movie.pvsm

                    
  4. create for each AVS-inp file a corresponding state-file:

     > ls -1 AVS_movie_0*.inp | gawk '{print "sed -e \"s/AVS_movie_009000.inp/"$1"/g\"< paraview_movie.pvsm > "$1".pvsm"}' | sh
     
     and create a new file tmp.in:
     
     >  ls -1 AVS_movie_0*.inp.pvsm > tmp.in

     
  5. run the script:

     > ./paraview-run.py   
     
     which renders jpeg files: paraview_movie.000001.jpg,...
       

  6. make mpeg-movie (and resize images) using command from ImageMagick (http://www.imagemagick.org)

     > convert -delay 8 -resize 400x320 paraview_movie.*.jpg movie.mpg

