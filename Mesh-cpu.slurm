#!/bin/bash -l
#SBATCH -J mesh_UCB
#SBATCH --nodes=12
#SBATCH --ntasks-per-node=32
#SBATCH --ntasks=384
#SBATCH --output log-mesh-%J.out
#SBATCH --error  log-mesh-%J.err
#SBATCH --dependency=afterok:6369253
##SBATCH --time=00:30:00
#SBATCH --time=01:30:00


#cd $SLURM_SUBMIT_DIR

module purge
module load intel/19.0/64/19.0.5.281 intel-mpi/intel/2019.7/64
module list
# change directory to build
#cd ../specfem
#cd ../specfem-body

#Run Mesher
ulimit -s unlimited
srun  ./bin/xmeshfem3D
#srun  --mpi=pmix  ./bin/xmeshfem3D
#srun /tigress/mw1685/Share-FILES/TEST-Berkeley-Model/specfem/bin/xmeshfem3D
#srun  --mpi=pmix -n 64 ../specfem/bin/xmeshfem3D
