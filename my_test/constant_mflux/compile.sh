#
mkdir -p fixedflux_periodic_dr1000
#
# compile
#
#CC99='mpicc -std=c99' qcc -grid=quadtree -events -catch -autolink -progress -O2 -g -Wall -pipe -Ddimension=2 -D_FORTIFY_SOURCE=2 -D_MPI=1 -o fixedflux_periodic_dr1000/fixedflux_periodic_dr1000 fixedflux_periodic_dr1000.c -L$BASILISK/gl -lfb_tiny -lm
CC99='mpicc -std=c99' qcc -grid=quadtree -events -autolink -progress -O2 -g -Wall -pipe -Ddimension=2 -D_FORTIFY_SOURCE=2 -D_MPI=1 -o fixedflux_periodic_dr1000/fixedflux_periodic_dr1000 fixedflux_periodic_dr1000.c -L$BASILISK/gl -lfb_tiny -lm

