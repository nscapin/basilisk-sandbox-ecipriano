#
mkdir -p forcedbi
#
# compile
#
CC99='mpicc -std=c99' qcc -grid=quadtree -events -autolink -progress -O2 -Wall -pipe -D_FORTIFY_SOURCE=2 -D_MPI=1 -o forcedbi/forcedbi forcedbi.c -L$BASILISK/gl -lfb_tiny -lm
#qcc -grid=quadtree -events -autolink -progress -O2 -g -Wall -pipe -D_FORTIFY_SOURCE=2 -o water/water water.c -L$BASILISK/gl -lfb_tiny -lm
#qcc -grid=quadtree -events -autolink -progress -O2 -Wall -pipe -Ddimension=2 -D_FORTIFY_SOURCE=2 -D_MPI=0 -o water/water water.c -L$BASILISK/gl -lfb_tiny -lm
#qcc -grid=quadtree -events -autolink -Wall -g -DTRASH=1 -o water/water water.c -L$BASILISK/gl -lfb_tiny -lm
#qcc -grid=quadtree -events -autolink -Wall -g -catch -o water/water water.c -L$BASILISK/gl -lfb_tiny -lm
#qcc -grid=quadtree -events -autolink -progress -Wall -O2 -o water/water water.c -L$BASILISK/gl -lfb_tiny -lm


