/**
# Non-Isothermal Evaporation of a water droplet

In this test case we want to simulate the evaporation of a
n-heptane droplet in nitrogen. The droplet is pure in a non
isothermal environment. Therefore, we want to solve both
the chemical species mass fractions and the temperature
field. At the beginning of the simulation, the droplet
is initialized at 363K, while the ambient temperature is 565K.
The heat conduction from the environment heats up the liquid
droplet increasing the vapor pressure value and, therefore,
increasing the evaporation rate of the droplet. The
interface temperature tends to a plateau, given by the
interplay between the heat conduction from the environment
and the evaporation process that cools down the interface.

![Evolution of the temperature field (left) and the n-heptane mass fraction (right)](c7pathak/movie.mp4)(height=400 width=900)
*/

/**
## Phase Change Setup

We define the number of gas and liquid species in the domain,
we filter the density field to reduce problems related with
the strong density ratio. The interface temperature must be
obtained from the root-finding procedure with the *fsolve()*
function, using the GSL interface. The Antoine equation is
used according to the thermodynamic equilibrium implemented
in [Pathak et al., 2018](#pathak2018steady), which proposed
this test case. */

#define CLOSED_DOMAIN

#define NGS 2
#define NLS 1

#define FILTERED
#define SOLVE_TEMPERATURE
#define USE_GSL 0
#define USE_ANTOINE

/**
## Simulation Setup

We use the centered solver with the evaporation source term
in the projection step. The extended velocity is obtained
from the doubled pressure-velocity coupling. We use the
evaporation model together with the multiomponent phase
change mechanism. */

#include "../../src/navier-stokes/centered-evaporation.h"
#include "../../src/navier-stokes/centered-doubled.h"
#include "two-phase.h"
//#include "../../src/navier-stokes/conserving-evaporation.h"
#include "tension.h"
#include "../../src/evaporation.h"
#include "../../src/multicomponent.h"
#include "view.h"

/**
### Data for multicomponent model

We define the data required by the multicomponent phase
change mechanism, including the solution of the temperature
field. The equilibrium constant *inKeq* is ignored when
*USE_ANTOINE* or *USE_CLAPEYRON* is used. */

char* gas_species[NGS] = {"H2O", "N2"};
char* liq_species[NLS] = {"H2O"};
char* inert_species[1] = {"N2"};
double gas_start[NGS] = {0., 1.};
double liq_start[NLS] = {1.};
double inDmix1[NLS] = {0.};
double inDmix2[NGS] = {2.49e-5, 2.49e-5};
double inKeq[NLS] = {0.};
double Tboil[NLS] = {373.15};

double lambda1 = 0.60700;
double lambda2 = 0.02624;
double dhev = 2.44e6;
double cp1 = 4137.00;
double cp2 = 1040.00;
double TL0 = 278.15;
double TG0 = 301.15;

/**
### Boundary conditions

Outflow boundary conditions are set at the top and right
sides of the domain. */

#ifdef CLOSED_DOMAIN
#else
u.n[top] = neumann (0.);
u.t[top] = neumann (0.);
p[top] = dirichlet (0.);
uext.n[top] = neumann (0.);
uext.t[top] = neumann (0.);
pext[top] = dirichlet (0.);

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
p[right] = dirichlet (0.);
uext.n[right] = neumann (0.);
uext.t[right] = neumann (0.);
pext[right] = dirichlet (0.);

u.n[left] = neumann (0.);
u.t[left] = neumann (0.);
p[left] = dirichlet (0.);
uext.n[left] = neumann (0.);
uext.t[left] = neumann (0.);
pext[left] = dirichlet (0.);

u.n[bottom] = neumann (0.);
u.t[bottom] = neumann (0.);
p[bottom] = dirichlet (0.);
uext.n[bottom] = neumann (0.);
uext.t[bottom] = neumann (0.);
pext[bottom] = dirichlet (0.);
#endif

/**
### Simulation Data

We declare the maximum and minimum levels of refinement,
the initial radius and diameter, and the radius from the
numerical simulation. */

int maxlevel, minlevel = 2;
double D0 = 5.e-6, effective_radius0;

int main (void) {
  
  /**
  We set the material properties of the fluids: air and water. */

  rho1 = 1000.0; rho2 = 1.125;
  mu1  = 4.45e-4; mu2 = 1.78e-5;
  Pref = 101325.0;

  /**
  We change the dimension of the domain as a function
  of the initial diameter of the droplet. */

  L0 = 2.*D0;
  origin (0.0,0.0);

  /**
  We change the surface tension coefficient. and we
  decrease the tolerance of the Poisson solver. */

  f.sigma = 0.01;

#ifdef CLOSED_DOMAIN
  foreach_dimension() {
    periodic (right);
  }
#endif

  /**
  We run the simulation at different maximum
  levels of refinement. */

  init_grid (1 << 4);
  for (maxlevel = 5; maxlevel <= 7; maxlevel++) {
    init_grid (1 << maxlevel);
    run();
  }
}

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

/**
We initialize the volume fraction field and we compute the
initial radius of the droplet. We don't use the value D0
because for small errors of initialization the squared
diameter decay would not start from 1. */

event init (i = 0) {

  refine (sq(x-0.5*L0) + sq(y-0.5*L0) < sq(0.80*D0) && 
	  sq(x-0.5*L0) + sq(y-0.5*L0) > sq(0.20*D0) &&
	  level < maxlevel);
  fraction (f, circle (x-0.5*L0, y-0.5*L0, 0.5*D0));
  effective_radius0 = pow(3.*statsf(f).sum, 1./3.);

  /**
  We set the molecular weights of the chemial species
  involved in the simulation (by default inMW=1). */

  inMW[0] = 18.0; inMW[1] = 29.0;

  /**
  The proper Antoine equation function must be set
  to the attribute *antoine* of the liquid phase mass
  fraction fields. */

  scalar YL = YLList[0];
  YL.antoine = antoine_heptane;

}

/**
We use the same boundary conditions used by
[Pathak at al., 2018](#pathak2018steady). */

event bcs (i = 0) {

#ifdef CLOSED_DOMAIN
#else
  scalar H2O = YGList[0];
  scalar N2  = YGList[1];

  H2O[top] = dirichlet (0.);
  H2O[right] = dirichlet (0.);
  H2O[left] = dirichlet (0.);
  H2O[bottom] = dirichlet (0.);

  N2[top] = dirichlet (1.);
  N2[right] = dirichlet (1.);
  N2[left] = dirichlet (1.);
  N2[bottom] = dirichlet (1.);

  TG[top] = dirichlet (TG0);
  TG[right] = dirichlet (TG0);
  TG[left] = dirichlet (TG0);
  TG[bottom] = dirichlet (TG0);
  /*
  TL[top] = dirichlet (TG0);
  TL[right] = dirichlet (TG0);
  TL[left] = dirichlet (TG0);
  TL[bottom] = dirichlet (TG0);
  */
#endif

}

/**
We adapt the grid according to the mass fractions of the
mass fraction of n-heptane, the temperature, and the
velocity field. */

event adapt (i++) {
  scalar H2O = YList[0];
  adapt_wavelet_leave_interface ({H2O,T,u.x,u.y}, {f},
      (double[]){1.e-3,1.e-1,1.e-1,1.e-1,1.e-3}, maxlevel, minlevel, 1);
}

/**
## Post-Processing

The following lines of code are for post-processing purposes. */

/**
### Output Files

We write on a file the squared diameter decay and the
dimensionless time. */

event output_data (i++) {
  char name[80];
  sprintf (name, "OutputData-%d", maxlevel);
  static FILE * fp = fopen (name, "w");

  scalar YInt_H2O = YGIntList[0];
  scalar Y_H2O = YList[0];
  double effective_radius = pow(3.*statsf(f).sum, 1./3.);
  double d_over_d02 = sq (effective_radius / effective_radius0);

  double TIntavg = avg_interface (TInt, f);
  double YIntavg = avg_interface (YInt_H2O, f);
  double Tavg = avg_interface (T, f);
  double Yavg = avg_interface (Y_H2O, f);

  fprintf (fp, "%g %g %g %g %g %g %g\n",
      t, effective_radius, d_over_d02, TIntavg, YIntavg, Tavg, Yavg);
  fflush (fp);
}

/**
### Temperature and Mass Fraction Profiles

We write on a file the temperature and mass fraction
profiles at different time instants. */

event profiles (t = {3.29e-6, 3.e-5, 1.05e-4, 1.5e-4}) {
  char name[80];
  sprintf (name, "Profiles-%d", maxlevel);

  /**
  We create an array with the temperature and mass
  fraction profiles for each processor. */

  scalar H20 = YList[0];

  Array * arrtemps = array_new();
  Array * arrmassf = array_new();
  for (double x = 0.; x < L0; x += 0.5*L0/(1 << maxlevel)) {
    double valt = interpolate (T, x, L0/2.0);
    double valm = interpolate (H20, x, L0/2.0);
    valt = (valt == nodata) ? 0. : valt;
    valm = (valm == nodata) ? 0. : valm;
    array_append (arrtemps, &valt, sizeof(double));
    array_append (arrmassf, &valm, sizeof(double));
  }
  double * temps = (double *)arrtemps->p;
  double * massf = (double *)arrmassf->p;

  /**
  We sum each element of the arrays in every processor. */

  @if _MPI
  int size = arrtemps->len/sizeof(double);
  MPI_Allreduce (MPI_IN_PLACE, temps, size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, massf, size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  @endif

  /**
  The master node writes the profiles on a file. */

  if (pid() == 0) {
    char name[80];
    sprintf (name, "Profiles-%d", maxlevel);
    FILE * fpp = fopen (name, "a");
    int count = 0;
    for (double x = 0.; x < L0; x += 0.5*L0/(1 << maxlevel)) {
      fprintf (fpp, "%g %g %g\n", x, temps[count], massf[count]);
      count++;
    }
    fprintf (fpp, "\n\n");
    fclose (fpp);
  }
  array_free (arrtemps);
  array_free (arrmassf);

}

/**
### Movie

We write the animation with the evolution of the
n-heptane mass fraction, the interface position
and the temperature field. */

#define POPEN(name, mode) fopen (name ".ppm", mode)

//event movie (t += 2.e-6; t <= 1.6e-4) {
//event movie (t += 2.e-6; t <= 3.2e-4) {
event movie (i += 10) {
  clear();
  box();
  view (tx = -0.5, ty = -0.5);
  draw_vof ("f");
  squares ("H2O", min = 0., max = 1., linear = true);
  {
    static FILE * fp = POPEN ("movie_sca", "a");
    save (fp = fp);
  }
  clear();
  box();
  view (tx = -0.5, ty = -0.5);
  draw_vof ("f");
  squares ("T", min = statsf(T).min, max = TG0, linear = true);
  {
    static FILE * fp = POPEN ("movie_tmp", "a");
    save (fp = fp);
  }
}

/**
## Results

The numerical results are compared with the results obtained
by [Pathak et al., 2018](#pathak2018steady) using a radially
symmetric model, integrated using an ODE solver.

~~~gnuplot Evolution of the squared diameter decay
reset
set xlabel "t [s]"
set ylabel "(D/D_0)^2"
set key top right
set grid

plot "../data/pathak-heptane-T563-diam.csv" w p ps 2 t "Pathank et al., 2018", \
     "OutputData-6" u 1:3 w l lw 2 t "LEVEL 6"
~~~

~~~gnuplot Evolution of the interface temperature
reset
set xlabel "t [s]"
set ylabel "Interface Temperature [K]"
set key bottom right
set grid

plot "../data/pathak-heptane-T563-temp.csv" w p ps 2 t "Pathank et al., 2018", \
     "OutputData-6" u 1:4 w l lw 2 t "LEVEL 6"
~~~

~~~gnuplot Evolution of the temperature profiles
reset
set xlabel "radius [m] x10^{6}"
set ylabel "Temperature [K]"
set key bottom right
set grid

plot "../data/pathak-heptane-T563-Tprofile-329e-6.csv" w p pt 8 lc 1 t "time = 3.29x10^{-6} s", \
     "../data/pathak-heptane-T563-Tprofile-3e-5.csv"   w p pt 8 lc 2 t "time = 3.00x10^{-5} s", \
     "../data/pathak-heptane-T563-Tprofile-105e-4.csv" w p pt 8 lc 3 t "time = 1.05x10^{-4} s", \
     "../data/pathak-heptane-T563-Tprofile-150e-4.csv" w p pt 8 lc 4 t "time = 1.50x10^{-4} s", \
     "Profiles-6" index 0 u ($1*1e+6):2 w l lw 2 lc 1 notitle, \
     "Profiles-6" index 1 u ($1*1e+6):2 w l lw 2 lc 2 notitle, \
     "Profiles-6" index 2 u ($1*1e+6):2 w l lw 2 lc 3 notitle, \
     "Profiles-6" index 3 u ($1*1e+6):2 w l lw 2 lc 4 notitle
~~~

~~~gnuplot Evolution of the n-heptane mass fraction profiles
reset
set xlabel "radius [m] x10^{6}"
set ylabel "Mass Fraction [-]"
set key top right
set grid

plot "../data/pathak-heptane-T563-Yprofile-329e-6.csv" w p pt 8 lc 1 t "time = 3.29x10^{-6} s", \
     "../data/pathak-heptane-T563-Yprofile-3e-5.csv"   w p pt 8 lc 2 t "time = 3.00x10^{-5} s", \
     "../data/pathak-heptane-T563-Yprofile-105e-4.csv" w p pt 8 lc 3 t "time = 1.05x10^{-4} s", \
     "../data/pathak-heptane-T563-Yprofile-150e-4.csv" w p pt 8 lc 4 t "time = 1.50x10^{-4} s", \
     "Profiles-6" index 0 u ($1*1e+6):3 w l lw 2 lc 1 notitle, \
     "Profiles-6" index 1 u ($1*1e+6):3 w l lw 2 lc 2 notitle, \
     "Profiles-6" index 2 u ($1*1e+6):3 w l lw 2 lc 3 notitle, \
     "Profiles-6" index 3 u ($1*1e+6):3 w l lw 2 lc 4 notitle
~~~

## References

~~~bib
@article{pathak2018steady,
  title={Steady-state and transient solutions to drop evaporation in a finite domain: Alternative benchmarks to the d2 law},
  author={Pathak, Ashish and Raessi, Mehdi},
  journal={International Journal of Heat and Mass Transfer},
  volume={127},
  pages={1147--1158},
  year={2018},
  publisher={Elsevier}
}
~~~
*/
