#!/bin/bash
#SBATCH --job-name=05500         # create a short name for your job
#SBATCH --nodes=1                   # node count
#SBATCH --ntasks-per-node=96        # number of tasks per node
#SBATCH --cpus-per-task=1           # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=4G            # memory per cpu-core (4G is default)
#SBATCH --time=23:59:00             # total run time limit (HH:MM:SS)
#
module load intel-mpi/gcc/2021.3.1
#
ulimit -s unlimited
#
RHO_RATIO=816.32;
MU_RATIO=25.0;
L0D0_RATIO=8.0;
D0=200.e-6;
sigma=0.1;
amp_force=5500;
FORCED=1.;
MAX_LEVEL=8;
MIN_LEVEL=5;
end_sim=1000.0;
t_phy=2.5e-4;
single_ph=1;
#
# Extract SBATCH time
#
SBATCH_TIME=$(grep -oP '(?<=#SBATCH --time=)[0-9:]*' $0)
#
srun hit_sp -m $SBATCH_TIME $RHO_RATIO $MU_RATIO $L0D0_RATIO 
                            $D0 $sigma $amp_force $FORCED
			    $MAX_LEVEL $MIN_LEVEL $end_sim $t_phy $single_ph > out.log 2>&1
#
