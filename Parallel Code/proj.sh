#!/bin/sh
#SBATCH --partition=general-compute
#SBATCH --nodes=12
#SBATCH --time=00:30:00
#SBATCH --ntasks-per-node=12 
#SBATCH --exclusive 
#SBATCH --job-name="HPC_Proj"
#SBATCH --output=HPC_Proj.out
#SBATCH --mail-user=sbhatt2@buffalo.edu
#SBATCH --mail-type=END
#SBATCH --requeue
SLURM_NPROCS=`srun --nodes=${SLURM_NNODES} bash -c 'hostname' |wc -l`
echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURM_NPROCS"=$SLURM_NTASKs
echo "SLURMTMPDIR="$SLURMTMPDIR

cd $SLURM_SUBMIT_DIR
echo "working directory = "$SLURM_SUBMIT_DIR


module load intel/16.0
module load intel-mpi/5.1.1
module load mkl/11.3
module list

ulimit -s unlimited


echo -e "HPC Project: \n"
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so


make

for((i=12;i<=$SLURM_NPROCS;i+=12)); do
	echo -e "\n For $i Processes: \n"
	srun -n $i ./lap
done
 

echo "Completed!"	
