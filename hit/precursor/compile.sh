#
module load intel-mpi/gcc/2021.3.1
rm -rf hit_sp && mkdir -p hit_sp
#
# compile
#
#CC99='mpicc -std=c99' qcc -grid=octree -events -autolink -O2 -Wall -pipe -D_FORTIFY_SOURCE=2 -D_MPI=1 -o hit_sp/hit_sp hit.c -L$BASILISK/gl -lfb_tiny -lm
CC99='mpicc -std=c99' qcc -grid=octree -events -autolink -O2 -Wall -pipe -D_FORTIFY_SOURCE=2 -D_MPI=1 -o hit_sp/hit_sp hit.c -L$BASILISK/gl -lfb_tiny -lm
