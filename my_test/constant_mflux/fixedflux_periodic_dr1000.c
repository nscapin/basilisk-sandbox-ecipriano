/**
# Fixed Flux Droplet Evaporation

Evaporation of a liquid droplet with a fixed vaporization
flowrate. In these conditions, the droplet consumption
is obtained from a mass balance on the liquid droplet and
the analytic solution can be easily derived:
$$
  \dfrac{dR}{dt} = \dfrac{\dot{m}}{\rho_l}
$$

We want to test and compare the droplet consumption with
the exact solution, considering the presence of the Stefan
flow. The simulation setup used here was adapted from
[Malan et al, 2021](#malan2021geometric).
The animation shows the consuption of the liquid droplet,
that maintains a perfectly spherical shape throughtout the
entire lifetime, as well as the radial velocity field
caused by the phase change.
![Evolution of the interface position and the velocity field](fixedflux/movie.mp4)(width="400" height="400")
*/

/**
## Phase Change Setup

We divide the mass balance by the density of the gas phase,
according with the non-dimensional treatment reported in
the [paper](#malan2021geometric). There is no need to shift
the volume expansion term because we do not transport the
temperature or mass fractions in this test case.
*/

#define NOSHIFTING
#define CLOSED_DOMAIN

/**
## Simulation Setup

We use the centered Navier-Stokes equations solver with
the evaporation source term in the projection step. The
double pressure velocity coupling approach is adopted to
obtain a divergence-free velocity which can be used for
the VOF advection equation. The evaporation model is
combined with the fixed flux mechanism, that imposes
a constant vaporization flowrate.
*/

#include "../../src/navier-stokes/centered-evaporation.h"
#include "../../src/navier-stokes/centered-doubled.h"
#include "two-phase.h"
//#include "../../src/navier-stokes/conserving-evaporation.h"
#include "tension.h"
#include "../../src/evaporation.h"
#include "../../src/fixedflux.h"
#include "view.h"

/**
### Boundary Conditions

Here we have two options:

 1. Outflow boundary conditions are imposed on every
    domain boundary since the droplet is initially placed
    at the center of the domain. Because the
    centered-doubled method is used, the boundary conditions
    must be imposed also for the *extended* fields.
 2. Periodic or closed domanins: we set periodicity for all fields 
    and in all physical boundaries. Note a compensation term is 
    added to the rhs of the continuity equation  

*/

#ifdef CLOSED_DOMAIN
#else
  u.n[top] = neumann (0.);
  u.t[top] = neumann (0.);
  p[top] = dirichlet (0.);
  uext.n[top] = neumann (0.);
  uext.t[top] = neumann (0.);
  pext[top] = dirichlet (0.);
  
  u.n[bottom] = neumann (0.);
  u.t[bottom] = neumann (0.);
  p[bottom] = dirichlet (0.);
  uext.n[bottom] = neumann (0.);
  uext.t[bottom] = neumann (0.);
  pext[bottom] = dirichlet (0.);
  
  u.n[left] = neumann (0.);
  u.t[left] = neumann (0.);
  p[left] = dirichlet (0.);
  uext.n[left] = neumann (0.);
  uext.t[left] = neumann (0.);
  pext[left] = dirichlet (0.);
  
  u.n[right] = neumann (0.);
  u.t[right] = neumann (0.);
  p[right] = dirichlet (0.);
  uext.n[right] = neumann (0.);
  uext.t[right] = neumann (0.);
  pext[right] = dirichlet (0.);
#endif

/**
### Parameters Data
*/

int maxlevel;
double R0 = 0.125;
double mEvapVal = -0.10;
double RE = 25.0;
double WE = 0.1;
double MU_RATIO  = 50.0;
double RHO_RATIO = 1000.0;
double f_r = 500.0;

int main(void) {

  /**
  Set domain dimension. */

  L0 = 8.0*R0;
  
  /**
  We define some reference quantities. */

  rho2 = 1.0;
  double u_ref = fabs(mEvapVal/rho2);
  
  /**
  We set the liquid density and the viscosity values. */

  rho1 = rho2*RHO_RATIO;
  mu2  = rho2*u_ref*(2.0*R0)/RE;
  mu1  = mu2*MU_RATIO;

  /**
  The surface tension is set to zero in order to assess
  the sphericity of the droplet without being influenced
  by the surface tension. */

  f.sigma = 0.;

  /**
  We make lists with the levels of refinement and the
  corresponding maximum time steps. We need to impose the
  time steps because the surface tension is set to zero. */

  double mllist[] = {5, 6, 7, 8};
  double dtlist[] = {2.5e-3, 1.25e-3, 6.25e-4, 3.125e-4};
  //double mllist[] = {6};
  //double dtlist[] = {1.25e-3};

  /**
  We modify the Poisson equation solver tolerance and
  we create a folder where the facet outputs will be
  written. */

  system ("mkdir -p facets");
  TOLERANCE = 1.e-4;

#ifdef CLOSED_DOMAIN
  foreach_dimension() {
    periodic (right);
  }
#endif

  /**
  We perform the simulation at different levels of
  refinement. */

  for (int sim = 0; sim < 3; sim++) {
    maxlevel = mllist[sim];
    DT = dtlist[sim];
    TOLERANCE = 1.e-4;
#if TREE
    init_grid (1 << 4);
#else
    init_grid (1 << maxlevel);
#endif
    run();
  }

}

#define circle(x, y, R) (sq(R) - sq(x - 0.5*L0) - sq(y - 0.5*L0))

/**
We initialize the volume fraction field.
*/

event init (i = 0) {
#if TREE
  fprintf(stderr, "Simulation run at level on a tree grid %d\n", maxlevel);
  refine (sq(x-L0*0.5) + sq(y-L0*0.5) < sq(1.5*R0) && level < maxlevel);
#else
  fprintf(stderr, "Simulation run at level on a fixed grid %d\n", maxlevel);
#endif
  fraction (f, circle (x, y, R0));
}

/**
## Post-Processing

The following lines of code are for post-processing purposes.
*/

/**
### Exact Solution

We write a function that computes the analytic
solution for the droplet volume.
*/

double exact (double time) {
  return pi*sq (R0 + mEvapVal*time/rho1);
}

/**
### Output Files

We write the droplet volume from the numerical simulation
as well as the exact solution.
*/

event log_simulation (i++) {

  double tot_vol = 0, tot_div = 0, tot_div_ex = 0;
  foreach(reduction(+:tot_vol),reduction(+:tot_div),reduction(+:tot_div_ex)) {
    tot_vol += dv();
    foreach_dimension() {
      tot_div += dv()*(uf.x[1]-uf.x[0])/Delta;
    }
    tot_div_ex += stefanflow[]*dv();
  }
  tot_div    /= tot_vol;
  tot_div_ex /= tot_vol;

  char name[80];
  sprintf (name, "log_simulation.out");
  static FILE * fp = fopen (name, "w");
  fprintf (fp, "%8E %8E %8E %8E %8E %8E\n", t, 1.0*i, dt, tot_vol, tot_div, tot_div_ex);
  fflush (fp);

}

event output (i++) {
  char name[80];
  sprintf (name, "OutputData-%d.out", maxlevel);

  double dropvol = statsf(f).sum;
  double relerr = fabs (exact(t) - dropvol) / exact(t);

  static FILE * fp = fopen (name, "w");
  fprintf (fp, "%8E %8E %8E %8E\n", t, dropvol, exact(t), relerr);
  fflush (fp);
}

/**
### Facets

We write the vof facets in such a way that it works also
in parallel: every processor writes its own facets in
a different file and the files are then gathered in a
single facets output file, ready to be plotted.
*/

event output_facets (t += 1*f_r) {
  char names[80];
  sprintf (names, "interfaces%d", pid());
  FILE * fpf = fopen (names,"w");
  output_facets (f, fpf);
  fclose (fpf); fflush (fpf);
  char command[80];
  sprintf(command, "LC_ALL=C cat interfa* > facets/facets-%d-%.1f", maxlevel, t);
  system(command);
}

/**
### Movie and picture

We write the animation with the volume fraction field
and the velocity vectors.
*/

event movie (t += 0.125*f_r) {
  clear();
  view (tx = -0.5, ty = -0.5, width = 1000, height = 1000);
  draw_vof ("f", lw = 1.5);
  squares ("u.x", linear = false, n = {0,1}, alpha = -L0/2.0, map = cool_warm, spread = -1);
  box();
  char movie[99];
  sprintf (movie, "movie_%02d.ppm", maxlevel);
  {
    static FILE * fp = fopen (movie, "a");
    save (fp = fp);
  }
}

event output_pic (t += 1*f_r) {
  clear();
  view (tx = -0.5, ty = -0.5, width = 1000, height = 1000);
  draw_vof ("f", lw = 1.6);
  squares ("u.x", linear = false, n = {0,1}, alpha = -L0/2.0, map = cool_warm, spread = -1);
  box();
  cells (n = {0,1}, alpha = -L0/2.0);
  char init_0[99];
  sprintf (init_0, "grid_vel_%09d_%02d.ppm", i, maxlevel);
  save(init_0);
}

event dump (i += 100) {
  char dname[100];
  sprintf (dname, "dump_eva_%d.bin", i);
  dump (dname);
}


#if TREE
event adapt (i++) {
  double uemax = 1e-2;
  adapt_wavelet_leave_interface ({u.x,u.y}, {f}, (double[]){uemax,uemax,1e-3}, maxlevel, 5, 1);
}
#endif

event end_simulation (t = 2*f_r) {
  return 1;
}

/**
## Results

The droplet is consumed by the evaporation is a smooth way,
and the comparison between the numerical and the analytic
droplet show a very good agreement and the convergence to
the exact solution.

~~~gnuplot Droplet Volume
reset
set xlabel "t [s]"
set ylabel "Droplet Volume [m^3]"
set size square
set grid

plot "OutputData-6" every 50 u 1:3 w p t "Analytic", \
     "OutputData-4" u 1:2 w l t "LEVEL 4", \
     "OutputData-5" u 1:2 w l t "LEVEL 5", \
     "OutputData-6" u 1:2 w l t "LEVEL 6"
~~~

~~~gnuplot Relative Errors
reset

stats "OutputData-4" using 4 nooutput name "LEVEL4"
stats "OutputData-5" using 4 nooutput name "LEVEL5"
stats "OutputData-6" using 4 nooutput name "LEVEL6"

set print "Errors.csv"

print sprintf ("%d %.12f", 2**4, LEVEL4_mean)
print sprintf ("%d %.12f", 2**5, LEVEL5_mean)
print sprintf ("%d %.12f", 2**6, LEVEL6_mean)

unset print

reset
set xlabel "N"
set ylabel "Relative Error"

set logscale x 2
set logscale y

set xr[8:128]
set size square
set grid

plot "Errors.csv" w p pt 8 title "Results", \
  10*x**(-1) title "1^{st} order", \
  30*x**(-2) title "2^{nd} order"
~~~

~~~gnuplot Droplet Sphericity
reset
set size square
set xrange[0.5:0.75]
set yrange[0.5:0.75]
set grid

array r[3]
r[1] = 0.18
r[2] = 0.13
r[3] = 0.08

set print "Circles.dat"
do for [i=1:3] {
  print sprintf ("%g %g %g", 0.5, 0.5, r[i])
}
unset print

set style fill transparent solid 0.2 noborder

p \
  "facets/facets-4-1.0" w l lw 2 lc 1 t "LEVEL 4", \
  "facets/facets-4-2.0" w l lw 2 lc 1 notitle, \
  "facets/facets-4-3.0" w l lw 2 lc 1 notitle, \
  "facets/facets-5-1.0" w l lw 2 lc 2 t "LEVEL 5", \
  "facets/facets-5-2.0" w l lw 2 lc 2 notitle, \
  "facets/facets-5-3.0" w l lw 2 lc 2 notitle, \
  "facets/facets-6-1.0" w l lw 2 lc 3 t "LEVEL 6", \
  "facets/facets-6-2.0" w l lw 2 lc 3 notitle, \
  "facets/facets-6-3.0" w l lw 2 lc 3 notitle, \
  "Circles.dat" u 1:2:3:1 w circles t "Analytic"
~~~

## References

~~~bib
@article{malan2021geometric,
  title={A geometric VOF method for interface resolved phase change and conservative thermal energy advection},
  author={Malan, LC and Malan, Arnaud G and Zaleski, St{\'e}phane and Rousseau, PG},
  journal={Journal of Computational Physics},
  volume={426},
  pages={109920},
  year={2021},
  publisher={Elsevier}
}
~~~
*/
