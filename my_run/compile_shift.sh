#
mkdir -p fixedflux_shift
#
# compile
#
CC99='mpicc -std=c99' qcc -grid=quadtree -events -catch -autolink -progress -O2 -g -Wall -pipe -Ddimension=2 -D_FORTIFY_SOURCE=2 -D_MPI=1 -o fixedflux_shift/fixedflux_shift fixedflux_shift.c -L$BASILISK/gl -lfb_tiny -lm

