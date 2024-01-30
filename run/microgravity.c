/**
# Thermal Expansion of a Liquid Droplet

In this test case, we study the thermal expansion of a liquid droplet of
n-heptane, at different ambient temperatures and at different levels of
refinement. The aim is to evaluate the convergence of the
[multicomponent.h](multicomponent.h) solver with variable properties.

The evaporation module is used suppressing the phase change, in order to focus
on the thermal expansion only, and to avoid evaporation.

![Evolution of the temperature field](expansion/movie.mp4)
*/

/**
## Default Simulation Data

The following data can be overwritten using compilation flags
in order to study the sensitivity to these parameters running
different simulations in parallel. */

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

#ifndef TEMPERATURE
# define TEMPERATURE 773.
#endif

#ifndef PRESSURE
# define PRESSURE 10.
#endif

#ifndef DIAMETER
# define DIAMETER 1.e-3
#endif

#ifndef FUEL
# define FUEL NC7H16
#endif

#ifndef INERT
# define INERT N2
#endif

#ifndef KINFOLDER
# define KINFOLDER evaporation/n-heptane-in-nitrogen
#endif

/**
## Phase Change Setup

We define the number of gas and liquid species in the domain,
and we initialize all the properties necessary for the multicomponent phase
change model. The properties are set to null values because they are
overwritten by the variable properties formulation, which computes all the
physical properties as a function of the thermodynamic state of the mixture. */

#define NGS 2
#define NLS 1

char* gas_species[NGS] = {TOSTRING(FUEL), TOSTRING(INERT)};
char* liq_species[NLS] = {TOSTRING(FUEL)};
char* inert_species[1] = {TOSTRING(INERT)};
double gas_start[NGS] = {0., 1.};
double liq_start[NLS] = {1.};
double inDmix1[NLS] = {0.};
double inDmix2[NGS] = {0., 0.};
double inKeq[NLS] = {0.};

double lambda1 = 0.;
double lambda2 = 0.;
double dhev = 0.;
double cp1 = 0.;
double cp2 = 0.;

/**
We set the initial temperature of the liquid and of the gas phase. */

double TL0 = 300.;
double TG0 = TEMPERATURE;

/**
We solve the temperature field, with variable properties, and we reduce the
tolerance for the calculation of the variable properties. The interfacial
temperature is not computed from the jump conditon, it is just set to the gas
phase temperature value. */

#define SOLVE_TEMPERATURE
#define VARPROP
#define USE_GSL 0
#define FSOLVE_ABSTOL 1.e-3
#define USE_ANTOINE_OPENSMOKE
#define FICK_CORRECTED
#define MOLAR_DIFFUSION
#define MASS_DIFFUSION_ENTHALPY
#define RADIATION_INTERFACE 0.93

/**
## Simulation Setup

We use the centered solver with the evaporation source term
in the projection step. The extended velocity is obtained
from the doubled pressure-velocity coupling. We use the
evaporation model together with the multiomponent phase
change mechanism.

We use the centered solver with the divergence source term in the projection
step. The calculation of the extended velocity can be skipped, because no phase
change is present. OpenSMOKE++ is used for the variable properties calculation. */

#include "axi.h"
#include "navier-stokes/centered-evaporation.h"
#include "navier-stokes/centered-doubled.h"
#include "opensmoke-properties.h"
#include "two-phase-varprop.h"
#include "tension.h"
#include "evaporation-varprop.h"
#include "multicomponent.h"
#include "view.h"

/**
### Boundary conditions

Outflow boundary conditions are set at the top and right
sides of the domain. */

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

/**
### Simulation Data

We declare the maximum and minimum levels of refinement,
the initial radius and diameter, and the radius from the
numerical simulation, and additional data for post-processing. */

int maxlevel, minlevel = 2;
double D0 = DIAMETER, effective_radius0, d_over_d02 = 1., tad = 0.;
double volumecorr = 0., volume0 = 0.;

vector ur[];

int main (void) {
  /**
  We set the kinetics folder, which defines the species
  of the simulation, and it is used by OpenSMOKE++ for the
  calculation of the thermodynamic and transport properties. */

  kinfolder = TOSTRING(KINFOLDER);

  /**
  We set additional data for the simulation. */

  rho1 = 0.; rho2 = 0.;
  mu1 = 0.; mu2 = 0.;
  Pref = PRESSURE*101325.;

  /**
  We change the dimension of the domain as a function
  of the initial diameter of the droplet. */

  double RR = 7.986462e+01;
  L0 = 0.5*RR*D0;

  /**
  We change the surface tension coefficient. */

  f.sigma = 0.03;

  /**
  We run the simulation at different maximum
  levels of refinement. */

  for (maxlevel = 10; maxlevel <= 12; maxlevel++) {
    init_grid (1 << 9);
    run();
  }
}

#define circle(x,y,R)(sq(R) - sq(x) - sq(y))

/**
We initialize the volume fraction field and we compute the
initial radius of the droplet. We don't use the value D0
because for small errors of initialization the squared
diameter decay would not start from 1. */

double mLiq0 = 0.;

event init (i = 0) {
  refine (circle (x, y, 4.*D0) > 0. && level < maxlevel);
  fraction (f, circle (x, y, 0.5*D0));

  /**
  We compute initial variables useful for post-processing. */

  effective_radius0 = pow (3.*statsf(f).sum, 1./3.);
  volume0 = 4./3.*pi*pow (effective_radius0, 3.);

  foreach (reduction(+:mLiq0))
    mLiq0 += rho1v[]*f[]*dv();

  /**
  We set the molecular weights of the chemial species
  involved in the simulation (by default inMW=1). */

  foreach_elem (YGList, jj) {
    inMW[jj] = OpenSMOKE_MW (jj);
    fprintf (stdout, "inMW[%d] = %g\n", jj, inMW[jj]);
  }
}

/**
We use the same boundary conditions used by
[Pathak at al., 2018](#pathak2018steady). */

event bcs (i = 0) {
  scalar fuel  = YGList[OpenSMOKE_IndexOfSpecies (TOSTRING(FUEL))];
  scalar inert = YGList[OpenSMOKE_IndexOfSpecies (TOSTRING(INERT))];

  fuel[top] = dirichlet (0.);
  fuel[right] = dirichlet (0.);

  inert[top] = dirichlet (1.);
  inert[right] = dirichlet (1.);

  TG[top] = dirichlet (TG0);
  TG[right] = dirichlet (TG0);
}

/**
We adapt the grid according to the mass fractions of the
mass fraction of n-heptane, the temperature, and the
velocity field. */

#if TREE
event adapt (i++) {
  scalar fuel = YList[OpenSMOKE_IndexOfSpecies (TOSTRING(FUEL))];
  adapt_wavelet_leave_interface ({fuel,T,u.x,u.y}, {f},
      (double[]){1.e-2,1.e-1,1.e-1,1.e-1}, maxlevel, minlevel, 1);
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

  double effective_radius = pow (3.*statsf(f).sum, 1./3.);
  double d_over_d02_old = d_over_d02;
  double tad_old = tad;

  d_over_d02 = sq (effective_radius / effective_radius0);
  tad = t/sq(D0*1e3);

  /**
  The vaporization rate is computed according to the formula
  in Liu & Avedisian, 2011, pag. 777 bottom. */

  double kv = 0.;
  if (i > 1)
    kv = fabs ((d_over_d02 - d_over_d02_old)/(tad - tad_old));

  double mLiq = 0.;
  foreach(reduction(+:mLiq))
    mLiq += rho1v[]*f[]*dv();

  fprintf (fp, "%g %g %g %g %g %g\n", t, tad, effective_radius,
      d_over_d02, mLiq/mLiq0, kv);
}

/**
We reconstruct a one-field discontinuous velocity field for
visualization. */

event end_timestep (i++) {
  foreach()
    foreach_dimension()
      ur.x[] = f[]*uext.x[] + (1. - f[])*u.x[];
}

/**
### Movie

We write the animation with the evolution of the
n-heptane mass fraction, the interface position
and the temperature field. */

event pictures (t = {0.1, 1., 2.}) {
  char name[80];
  sprintf (name, "picturefields-%.1f", t);

  FILE * fp = fopen (name, "w");
  output_field ({f,T,u.x,u.y}, fp, linear = true);
  fclose (fp);
}

#if MOVIE
event movie (t += 0.01; t <= 10.) {
  clear();
  box();
  view (fov = 3, samples = 2);
  draw_vof ("f", lw = 1.5);
  squares ("T", min = TL0, max = TG0, linear = true);
  save ("temperature.mp4");

  clear();
  box();
  view (fov = 3, samples = 2);
  draw_vof ("f", lw = 1.5);
  squares (TOSTRING(FUEL), min = 0., max = 1., linear = true);
  save ("fuel.mp4");
}
#endif

/**
### Snapshots

Output dump files for restore or post-processing. */

#if DUMP
event snapshots (t += 0.1) {
  char name[80];
  sprintf (name, "snapshots-%g", t);
  dump (name);
}
#endif

/**
### Stopping Condition

We stop the simulation when the droplet is almost fully consumed. */

event stop (i++) {
  if (d_over_d02 <= 0.05)
    return 1;
}

/**
## Results

~~~gnuplot Effect of Temperature on the Square Diameter Decay at P = 0.1 MPa
reset
set grid
set key top right
set xlabel "t/D_0^2 [s/mm^2]"
set ylabel "(D/D_0)^2 [-]"

plot "../microgravity-T471K-P1atm/OutputData-10" u 2:4 w l lw 2 lc 1 t "", \
     "../microgravity-T555K-P1atm/OutputData-10" u 2:4 w l lw 2 lc 2 t "", \
     "../microgravity-T648K-P1atm/OutputData-10" u 2:4 w l lw 2 lc 3 t "", \
     "../microgravity-T741K-P1atm/OutputData-10" u 2:4 w l lw 2 lc 4 t "", \
     "/home/chimica2/ecipriano/ExperimentalData/Nomura/nomura-heptane-1atm-471K.csv" w p lc 1 t "471 K", \
     "/home/chimica2/ecipriano/ExperimentalData/Nomura/nomura-heptane-1atm-555K.csv" w p lc 2 t "555 K", \
     "/home/chimica2/ecipriano/ExperimentalData/Nomura/nomura-heptane-1atm-648K.csv" w p lc 3 t "648 K", \
     "/home/chimica2/ecipriano/ExperimentalData/Nomura/nomura-heptane-1atm-741K.csv" w p lc 4 t "741 K"
~~~

~~~gnuplot Effect of Temperature on the Square Diameter Decay at P = 0.5 MPa
reset
set grid
set key top right
set xlabel "t/D_0^2 [s/mm^2]"
set ylabel "(D/D_0)^2 [-]"

plot "../microgravity-T468K-P5atm/OutputData-10" u 2:4 w l lw 2 lc 1 t "", \
     "../microgravity-T556K-P5atm/OutputData-10" u 2:4 w l lw 2 lc 2 t "", \
     "../microgravity-T655K-P5atm/OutputData-10" u 2:4 w l lw 2 lc 3 t "", \
     "../microgravity-T749K-P5atm/OutputData-10" u 2:4 w l lw 2 lc 4 t "", \
     "/home/chimica2/ecipriano/ExperimentalData/Nomura/nomura-heptane-5atm-468K.csv" w p lc 1 t "468 K", \
     "/home/chimica2/ecipriano/ExperimentalData/Nomura/nomura-heptane-5atm-556K.csv" w p lc 2 t "556 K", \
     "/home/chimica2/ecipriano/ExperimentalData/Nomura/nomura-heptane-5atm-655K.csv" w p lc 3 t "655 K", \
     "/home/chimica2/ecipriano/ExperimentalData/Nomura/nomura-heptane-5atm-749K.csv" w p lc 4 t "749 K"
~~~

~~~gnuplot Effect of Temperature on the Square Diameter Decay at P = 1.0 MPa
reset
set grid
set key top right
set xlabel "t/D_0^2 [s/mm^2]"
set ylabel "(D/D_0)^2 [-]"

plot "../microgravity-T466K-P10atm/OutputData-10" u 2:4 w l lw 2 lc 1 t "", \
     "../microgravity-T508K-P10atm/OutputData-10" u 2:4 w l lw 2 lc 2 t "", \
     "../microgravity-T669K-P10atm/OutputData-10" u 2:4 w l lw 2 lc 3 t "", \
     "../microgravity-T765K-P10atm/OutputData-10" u 2:4 w l lw 2 lc 4 t "", \
     "/home/chimica2/ecipriano/ExperimentalData/Nomura/nomura-heptane-10atm-466K.csv" w p lc 1 t "466 K", \
     "/home/chimica2/ecipriano/ExperimentalData/Nomura/nomura-heptane-10atm-508K.csv" w p lc 2 t "508 K", \
     "/home/chimica2/ecipriano/ExperimentalData/Nomura/nomura-heptane-10atm-669K.csv" w p lc 3 t "669 K", \
     "/home/chimica2/ecipriano/ExperimentalData/Nomura/nomura-heptane-10atm-765K.csv" w p lc 4 t "765 K"
~~~

~~~gnuplot Effect of Temperature on the Square Diameter Decay at P = 2.0 MPa
reset
set grid
set key top right
set xlabel "t/D_0^2 [s/mm^2]"
set ylabel "(D/D_0)^2 [-]"

plot "../microgravity-T452K-P20atm/OutputData-10" u 2:4 w l lw 2 lc 1 t "", \
     "../microgravity-T511K-P20atm/OutputData-10" u 2:4 w l lw 2 lc 2 t "", \
     "../microgravity-T656K-P20atm/OutputData-10" u 2:4 w l lw 2 lc 3 t "", \
     "../microgravity-T746K-P20atm/OutputData-10" u 2:4 w l lw 2 lc 4 t "", \
     "/home/chimica2/ecipriano/ExperimentalData/Nomura/nomura-heptane-20atm-452K.csv" w p lc 1 t "452 K", \
     "/home/chimica2/ecipriano/ExperimentalData/Nomura/nomura-heptane-20atm-511K.csv" w p lc 2 t "511 K", \
     "/home/chimica2/ecipriano/ExperimentalData/Nomura/nomura-heptane-20atm-656K.csv" w p lc 3 t "656 K", \
     "/home/chimica2/ecipriano/ExperimentalData/Nomura/nomura-heptane-20atm-746K.csv" w p lc 4 t "746 K"
~~~

~~~gnuplot Effect of Pressure on the Square Diameter Decay at T = 650 K
reset
set grid
set key top right
set xlabel "t/D_0^2 [s/mm^2]"
set ylabel "(D/D_0)^2 [-]"

plot "../microgravity-T648K-P1atm/OutputData-10"  u 2:4 w l lw 2 lc 1 t "", \
     "../microgravity-T655K-P5atm/OutputData-10"  u 2:4 w l lw 2 lc 2 t "", \
     "../microgravity-T661K-P10atm/OutputData-10" u 2:4 w l lw 2 lc 3 t "", \
     "../microgravity-T656K-P20atm/OutputData-10" u 2:4 w l lw 2 lc 4 t "", \
     "/home/chimica2/ecipriano/ExperimentalData/Nomura/nomura-heptane-1atm-648K.csv"  w p lc 1 t "0.1 MPa", \
     "/home/chimica2/ecipriano/ExperimentalData/Nomura/nomura-heptane-5atm-655K.csv"  w p lc 2 t "0.5 MPa", \
     "/home/chimica2/ecipriano/ExperimentalData/Nomura/nomura-heptane-10atm-661K.csv" w p lc 3 t "1.0 MPa", \
     "/home/chimica2/ecipriano/ExperimentalData/Nomura/nomura-heptane-20atm-656K.csv" w p lc 4 t "2.0 MPa"
~~~

*/

