#!/bin/bash

#SBATCH --job-name=GRChombo_Proca


## redirect output
#SBATCH --output="/home/hd/hd_hd/hd_pb293/JobScripts/GRChombo/SlurmOut/%x_%j.out"
#SBATCH --error="/home/hd/hd_hd/hd_pb293/JobScripts/GRChombo/SlurmOut/%x_%j.err"

## send email on updates
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE
#SBATCH --mail-user=shaundbfell@gmail.com

#############################################
## COMMENTED: set when calling sbatch
## An example command:  sbatch -J GRChombo_Proca -p multiple_il --mem=100gb -N 20 -t 48:00:00 --ntasks-per-node=8 --cpus-per-task=8 ~/JobScripts/GRChombo/ProcaSuperradiance/generic.slurm params_Dynamic_HighRes_AH.txt        
## -p is the queue
## -J is the job name
## -N is the number of compute nodes
## -t is the walltime limit
#############################################

###SBATCH --partition=multiple_il
###SBATCH --time=2-0
#### 12 nodes * 80 cores = ~1000 cores
###SBATCH --nodes=12
#### number of MPI ranks for node
#### I want 2 OMP threads per task
###SBATCH --ntasks-per-node=40
#### number of cpus per rank. Equals OMP threads per MPI rank
###SBATCH --cpus-per-task=2
#### memory per cpu
###SBATCH --mem-per-cpu=1950mb


export environment="testing"  ## values should be either 'prod' for production or 'testing' for testing

# CD to directory with exectable and params folder
cd "/home/hd/hd_hd/hd_pb293/Documents/Github/ProblemsWithProca/proca-on-kerr/Examples/GeneralizedProca"

module load compiler/gnu/10.2 mpi/openmpi lib/hdf5/1.12.2-gnu-10.2-openmpi-4.1 numlib/gsl/2.6-gnu-10.2 numlib/petsc/3.17.2-gnu-10.2-openmpi-4.1
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PETSC_DIR}/lib

# NO OVERLOAD allowed
export MPIRUN_OPTIONS="--bind-to core --map-by socket:PE=${OMP_NUM_THREADS}"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export NUM_CORES=${SLURM_NTASKS}*${OMP_NUM_THREADS}
## as suggested by Baden-Wurttemburg Cluster support ticket
export OMPI_MCA_btl_openib_if_exclude=mlx5_2


##### GRCHOMBO stuff now

# some magic here: replace the output path with a dynamically generated one with the job id
OUTPATH="/pfs/work7/workspace/scratch/hd_pb293-WS_GRChombo/${environment}/${SLURM_JOB_NAME}_${SLURM_JOB_ID}"
mkdir $OUTPATH

# PARAMETER FILE as 1st argument to the script, JUST NAME
PARAMFILE=$1

# replace the output_path line in params.txt
sed -i 's|.*output_path.*|output_path = '"\"$OUTPATH\""'|' params/$PARAMFILE



##prepare and copy visualization session files to newly created output directory

#path to generic files
VISITSESSIONFILE="/home/hd/hd_hd/hd_pb293/JobScripts/GRChombo/ProcaSuperradiance/visit.session"
VISITSESSIONGUIFILE="/home/hd/hd_hd/hd_pb293/JobScripts/GRChombo/ProcaSuperradiance/visit.session.gui"

#copy generic files to output path
cp $VISITSESSIONFILE $OUTPATH
cp $VISITSESSIONGUIFILE $OUTPATH

#replace generic parts of files with slurm job id
sed -i "s/XXXXXXXX/${SLURM_JOB_ID}/g" ${OUTPATH}/visit.session
sed -i "s/XXXXXXXX/${SLURM_JOB_ID}/g" ${OUTPATH}/visit.session.gui

## Create GRChombo directories already. This is necessary since we dynamically create the output directory, but restarting from a checkpoint file doesnt do this.
#HDF5SUBPATH=cat params/${PARAMFILE} | grep -E "^hdf5_subpath" | sed 's/[^"]*"\([^"]*\).*/\1/'
#POUTSUBPATH=cat params/${PARAMFILE} | grep -E "^pout_subpath" | sed 's/[^"]*"\([^"]*\).*/\1/'
#DATASUBPATH=cat params/${PARAMFILE} | grep -E "^data_subpath" | sed 's/[^"]*"\([^"]*\).*/\1/'
#mkdir "${OUTPATH}/${HDF5SUBPATH}"
#mkdir "${OUTPATH}/${POUTSUBPATH}"
#mkdir "${OUTPATH}/${DATASUBPATH}"
#echo "hdf5 subpath: ${HDF5SUBPATH}"
#echo "pout subpath: ${POUTSUBPATH}"
#echo "data subpath: ${DATASUBPATH}"



# cp parameters file to the output dir
cp params/$PARAMFILE $OUTPATH/params.txt

#export EXECUTABLE="./Main_GeneralizedProca3d.Linux.64.mpicxx.gfortran.DEBUG.MPI.OPENMPCC.ex ${PARAMFILE}"
export EXECUTABLE="./Main_GeneralizedProca3d.Linux.64.mpicxx.gfortran.MPI.OPENMPCC.ex params/${PARAMFILE}"

echo ""
echo ""
echo ""
echo "[ job ID = ${SLURM_JOB_ID} ]"
echo "[ running on $SLURM_PARTITION on nodes: $SLURM_NODELIST ]"
echo "[ with ${NUM_CORES} cores for ${SLURM_NTASKS} MPI tasks and ${OMP_NUM_THREADS} OMP threads each ]"
echo "[ now: $(date) ]"
echo ""
echo ""
echo "Parameter file:"
for line in "$(cat params/${PARAMFILE})"
do
	echo "$line";
done
echo ""
echo ""
echo ""

startexe="mpirun -n ${SLURM_NTASKS} ${MPIRUN_OPTIONS} ${EXECUTABLE}"

echo $startexe
exec $startexe
echo "simulation finished"
exit 0

