#!/bin/bash

#SBATCH --job-name=GRChombo_Proca_Fixed

#############################################
## COMMENTED: set when calling sbatch
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

## redirect output
#SBATCH --output="/home/hd/hd_hd/hd_pb293/JobScripts/GRChombo/SlurmOut/%x_%j.out"
#SBATCH --error="/home/hd/hd_hd/hd_pb293/JobScripts/GRChombo/SlurmOut/%x_%j.err"

## send email on updates
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE
#SBATCH --mail-user=shaundbfell@gmail.com

# CD to directory with exectable and params.txt file
cd /home/hd/hd_hd/hd_pb293/Documents/Github/ProblemsWithProca/proca-on-kerr/Examples/GeneralizedProca_FixedBG

module load compiler/gnu/10.2 mpi/openmpi lib/hdf5/1.12.2-gnu-10.2-openmpi-4.1 numlib/gsl/2.6-gnu-10.2 numlib/petsc/3.17.2-gnu-10.2-openmpi-4.1
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PETSC_DIR}/lib

# NO OVERLOAD allowed
export MPIRUN_OPTIONS="--bind-to core --map-by socket:PE=${OMP_NUM_THREADS} --report-bindings"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export NUM_CORES=${SLURM_NTASKS}*${OMP_NUM_THREADS}
## as suggested by BWCluster support ticket
export OMPI_MCA_btl_openib_if_exclude=mlx5_2

# GRCHOMBO stuff now

# some magic here: replace the output path with a dynamically generated one with the job id
OUTPATH="/pfs/work7/workspace/scratch/hd_pb293-WS_GRChombo/testing/${SLURM_JOB_NAME}_${SLURM_JOB_ID}"
mkdir $OUTPATH

# PARAMETER FILE as 1st argument to the script, JUST NAME
PARAMFILE=$1

# replace the output_path line in params.txt
sed -i 's|.*output_path.*|output_path = '"\"$OUTPATH\""'|' $PARAMFILE



##prepare and copy visualization session files to newly created output directory

#path to generic files
VISITSESSIONFILE="/home/hd/hd_hd/hd_pb293/JobScripts/GRChombo/ProcaSuperradiance/visit.session"
VISITSESSIONGUIFILE="/home/hd/hd_hd/hd_pb293/JobScripts/GRChombo/ProcaSuperradiance/visit.session.gui"

#copy generic files to output path
cp $VISITSESSIONFILE $OUTPATH
cp $VISITSESSIONGUIFILE $OUTPATH

#replace generic parts of files with slurm job id
sed -i "s/GRChombo_Proca_XXXXXXXX/Fixed_GRChombo_Proca_${SLURM_JOB_ID}/g" ${OUTPATH}/visit.session
sed -i "s/GRChombo_Proca_XXXXXXXX/Fixed_GRChombo_Proca_${SLURM_JOB_ID}/g" ${OUTPATH}/visit.session.gui
sed -i "s/GeneralizedProca/GeneralizedProcaFixed/g" ${OUTPATH}/visit.session
sed -i "s/GeneralizedProca/GeneralizedProcaFixed/g" ${OUTPATH}/visit.session.gui


# cp parameters file to the output dir
cp $PARAMFILE $OUTPATH/params.txt

#export EXECUTABLE="./Main_GeneralizedProca3d.Linux.64.mpicxx.gfortran.DEBUG.MPI.OPENMPCC.ex ${PARAMFILE}"
export EXECUTABLE="./Main_GeneralizedProca3d.Linux.64.mpicxx.gfortran.MPI.OPENMPCC.ex ${PARAMFILE}"

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
for line in "$(cat ${PARAMFILE})"
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

