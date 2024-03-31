/**
# Isothermal Evaporation of a Pure Liquid Droplet

In this test case we want to simulate the evaporation of a
pure droplet in an isothermal environment. Even if the
droplet is pure, we use the multicomponent model to test
that, if the multiple species have identical properties,
they behave like a pure phase.

Since the evaporation is isothermal we don't need to use
*CLAPEYRON* or *ANTOINE*, but we just fix the value of the
thermodynamic equilibrium constant. For a pure droplet,
this is equivaluent to fixing the value of the chemical
species mass fraction on the gas-phase side of the
interface.

The animation shows the map of the mass fraction field in
gas phase, at different levels of refinement. The simulation
setup was borrowed from [Pathak et al., 2018](#pathak2018steady).

![Evolution of the mass fraction field](pureisothermal/movie.mp4)(height=500 width=500)
*/

/**
## Phase Change Setup

We define the number of gas and liquid species in the domain,
we use the default evaporation setup. */

#define NGS 2
#define NLS 1
#define CLOSED_DOMAIN

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
change mechanism, without including the solution of the
temperature field. The equilibrium constant *inKeq* is
fixed and we assume that the droplet is made of 4 chemical
species with identical properties. */

char * gas_species[NGS] = {"A", "I"};
char * liq_species[NLS] = {"A"};
char * inert_species[1] = {"I"};
double gas_start[NGS] = {0.0, 1.0};
double liq_start[NLS] = {1.0};
double inDmix1[NLS] = {0.};
double inDmix2[NGS] = {2.e-3, 2.0e-3};
double inKeq[NLS] = {0.667};

/**
### Boundary conditions

Outflow boundary conditions are set at the top and right
sides of the domain. */

#ifdef CLOSED_DOMAIN
#else
u.n[top] = neumann (0.);
u.t[top] = neumann (0.);
p[top] = dirichlet (0.);
pf[top] = dirichlet_face (0.);
uext.n[top] = neumann (0.);
uext.t[top] = neumann (0.);
pext[top] = dirichlet (0.);
pfext[top] = dirichlet_face (0.);

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
p[right] = dirichlet (0.);
pf[right] = dirichlet_face (0.);
uext.n[right] = neumann (0.);
uext.t[right] = neumann (0.);
pext[right] = dirichlet (0.);
pfext[right] = dirichlet_face (0.);

u.n[bottom] = neumann (0.);
u.t[bottom] = neumann (0.);
p[bottom] = dirichlet (0.);
pf[bottom] = dirichlet_face (0.);
uext.n[bottom] = neumann (0.);
uext.t[bottom] = neumann (0.);
pext[bottom] = dirichlet (0.);
pfext[bottom] = dirichlet_face (0.);

u.n[left] = neumann (0.);
u.t[left] = neumann (0.);
p[left] = dirichlet (0.);
pf[left] = dirichlet_face (0.);
uext.n[left] = neumann (0.);
uext.t[left] = neumann (0.);
pext[left] = dirichlet (0.);
pfext[left] = dirichlet_face (0.);
#endif

/**
### Simulation Data

We declare the maximum and minimum levels of refinement,
the initial radius and diameter, and the radius from the
numerical simulation. */

int maxlevel, minlevel = 5;
double D0 = 0.4e-3, XC, YC;
double effective_radius0;

int main (void) {
  /**
  We set the material properties of the fluids.
  The density ratio is small, according to the
  value used in [Pathak et al., 2018](#pathak2018steady)
  in order to speed up the simulation. */

  rho1 = 10.; rho2 = 1.;
  mu2 = 1.e-3; mu1 = 1.e-4;

  /**
  We change the dimension of the domain as a function
  of the initial diameter of the droplet. */

  L0 = 2.*D0;
  origin(0.0,0.0);

  /**
  We change the surface tension coefficient. and we
  decrease the tolerance of the Poisson solver, and
  the maximum allowed time step. */

#ifdef CLOSED_DOMAIN

  /**
  Set bc. */

  foreach_dimension() {
    periodic(right);
  }

#endif


  f.sigma = 0.01;
  DT = 5.e-8;

  /**
  We run the simulation at different maximum
  levels of refinement. */

  for (maxlevel = 5; maxlevel <= 7; maxlevel++) {
    init_grid (1 << maxlevel);
    run();
  }
}

#define circle(x,y,R) (sq(R) - sq(x) - sq(y))

/**
We initialize the volume fraction field and we compute the
initial radius of the droplet. We don't use the value D0
because for small errors of initialization the squared
diameter decay would not start from 1. */

event init (i = 0) {
  fraction (f, circle(x-0.5*L0,y-0.5*L0,0.5*D0));
  effective_radius0 = sqrt (4./pi*statsf(f).sum);
}

/**
We set the boundary conditions for the mass fraction
fields in liquid phase and for the inert. */

event bcs (i = 0) {

#ifdef CLOSED_DOMAIN
#else
  for (scalar YG in YGList) {
    YG[top] = dirichlet (0.);
    YG[right] = dirichlet (0.);
    YG[bottom] = dirichlet (0.);
    YG[left] = dirichlet (0.);
  }

  scalar Inert = YGList[1];

  Inert[top] = dirichlet (1.);
  Inert[right] = dirichlet (1.);
  Inert[bottom] = dirichlet (1.);
  Inert[left] = dirichlet (1.);
#endif

}

/**
We adapt the grid according to the mass fractions of the
mass fraction of the chemical species and the velocity field.
*/

#if TREE
event adapt (i++) {
  scalar A = YList[0];
  adapt_wavelet_leave_interface ({A,u}, {f},
      (double[]){1.e-3,1.e-2,1.e-2}, maxlevel, minlevel, 1);
}
#endif

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

  double effective_radius = sqrt (4./pi*statsf(f).sum);
  double d_over_d02 = sq (effective_radius / effective_radius0);
  double tad = t*inDmix2[0]/sq(2.*effective_radius0);

  fprintf (fp, "%g %g %g\n", t, tad, d_over_d02);
  fflush (fp);
}

/**
### Mass Fraction Profiles

We write on a file the temperature and mass fraction
profiles at different time instants. */

event profiles (t = {1.03e-5, 6.03e-5, 1.40e-4}) {
  char name[80];
  sprintf (name, "Profiles-%d", maxlevel);
  FILE * fp = fopen (name, "a");

  /**
  We reconstruct the total mass fraction summing the
  contribution of the different liquid species. */

  scalar Ysum[];
  foreach() {
    Ysum[] = 0.;
    foreach_elem (YList, jj) {
      if (jj != inertIndex) {
        scalar Y = YList[jj];
        Ysum[] += Y[];
      }
    }
  }

  /**
  We create an array with the mass fraction profiles
  for each processor. */

  Array * arrmassf = array_new();
  for (double x = 0.; x < L0; x += 0.5*L0/(1 << maxlevel)) {
    double valm = interpolate (Ysum, x, L0/2.0);
    valm = (valm == nodata) ? 0. : valm;
    array_append (arrmassf, &valm, sizeof(double));
  }
  double * massf = (double *)arrmassf->p;

  /**
  We sum each element of the arrays in every processor. */

  @if _MPI
  int size = arrmassf->len/sizeof(double);
  MPI_Allreduce (MPI_IN_PLACE, massf, size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  @endif

  /**
  The master node writes the profiles on a file. */

  if (pid() == 0) {
    int count = 0;
    for (double x = 0.; x < L0; x += 0.5*L0/(1 << maxlevel)) {
      fprintf (fp, "%g %g\n", x, massf[count]);
      count++;
    }
    fprintf (fp, "\n\n");
    fflush (fp);
  }
  fclose (fp);
  array_free (arrmassf);
}

/**
### Movie

We write the animation with the evolution of the
chemical species mass fractions, the interface position
and the temperature field. */

event movie (t += 0.5e-5, t <= 1.5e-4) {
  clear ();
  view (tx = -0.5, ty = -0.5);
  draw_vof ("f", lw = 1.5);
  squares ("A_G", min = 0., max = liq_start[0]*inKeq[0], linear = true);
  char name[80];
  sprintf (name,"movie_%d.ppm", maxlevel);
  {
    static FILE * fp = fopen (name, "a");
    save (fp = fp);
  }
}
