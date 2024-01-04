/**
# Shape oscillation of an inviscid droplet

This test case is discussed in [Popinet, 2009](/src/references.bib#popinet2009).

A two-dimensional elliptical droplet (density ratio 1/1000) is
released in a large domain. Under the effect of surface-tension the
shape of the droplet oscillates around its (circular) equilibrium
shape. The fluids inside and outside the droplet are inviscid so
ideally no damping of the oscillations should occur. As illustrated on
the figures some damping occurs in the simulation due to numerical
dissipation.

This simulation is also a stringent test case of the accuracy of the
surface tension representation as no explicit viscosity can damp
eventual parasitic currents. 

We use either the momentum-conserving or standard Navier--Stokes
solver with VOF interface tracking and surface tension. */

#if MOMENTUM

# include "momentum.h"
# include "tension.h"
# define cf

#elif COMPRESSIBLE

# include "compressible/two-phase.h"
# include "compressible/Mie-Gruneisen.h"
# include "tension.h"
# include "compressible/tension.h"

double rho1, rho2;

#else // standard centered Navier--Stokes solver

# if GFM
#   include "navier-stokes/centered-gfm.h"
# else
#   include "navier-stokes/centered.h"
# endif
# define FILTERED 1
# include "two-phase.h"
# if !GFM
# include "tension.h"
# endif

#endif // standard centered Navier--Stokes solver

/**
The diameter of the droplet is 0.2. */

#define D 0.2

/**
We will vary the level of refinement to study convergence. */

FILE * fp = NULL;
int LEVEL;

int main()
{

  /**
  The density is variable. */

  rho1 = 1, rho2 = 1e-3;

#if COMPRESSIBLE
  /**
  The parameters of the EOS for the liquid and gas. */
  
  PI1 = 300.;
  gamma1 = 7.14;
  gamma2 = 1.4;
#if 0 // does not work
  PI2 = PI1;
  gamma2 = gamma1;
#endif
#endif

#if 0
  double nu1 = 1e-4, nu2 = 1e-4;
  mu1 = nu1*rho1;
  mu2 = nu2*rho2;
#endif
  
  /**
  The surface tension is unity. Decreasing the tolerance on the
  Poisson solve improves the results. We cleanup existing files and
  vary the level of refinement. */

  f.sigma = 1.;
#if GFM
  gfm.f = f;
  gfm.sigma = 1.;
#endif
  L0 = 0.5 [0];
  TOLERANCE = 1e-4 [*];
  remove ("error");
  remove ("laplace");
  for (LEVEL = 4; LEVEL <= 7; LEVEL++) {
    N = 1 << LEVEL;
    
    /**
    We open a file indexed by the level to store the time evolution of
    the kinetic energy. */

    char name[80];
    sprintf (name, "k-%d", LEVEL);
    fp = fopen (name, "w");
    run();
    fclose (fp);
  }

  /**
  We use *grep* to filter the lines generated by gnuplot containing
  the results of the fits (see below). */

  system ("grep ^fit out >> log");
}

event init (i = 0) {

  /**
  We initialise the shape of the interface, a slightly elliptic droplet. */

  fraction (f, D/2.*(1. + 0.05*cos(2.*atan2(y,x))) - sqrt(sq(x) + sq(y)));

#if COMPRESSIBLE
  /**
  For compressible flows we initialize the densities, pressure and
  energies. */
  
  // fixme: The initial condition for pressure/energy does not respect
  // the Laplace pressure because the curvature is not constant
  double p0L = 1, p0 = p0L + f.sigma/D*2;
  foreach() {
    frho1[] = rho1*f[];
    frho2[] = rho2*(1. - f[]);
    p[]   = p0*f[] + p0L*(1. - f[]);
    fE1[] = f[]*(p0/(gamma1 - 1.) + PI1*gamma1/(gamma1 - 1.));
    fE2[] = (1. - f[])*p0L/(gamma2 - 1.);
  }
#endif // COMPRESSIBLE
}

/**
At each timestep we output the kinetic energy. */

event logfile (i++; t <= 1) {
  double ke = 0.;
  foreach (reduction(+:ke))
#if MOMENTUM || COMPRESSIBLE
    ke += dv()*(sq(q.x[]) + sq(q.y[]))/rho[];
#else
    ke += dv()*(sq(u.x[]) + sq(u.y[]))*rho(sf[]);
#endif
  fprintf (fp, "%g %g %d\n", t, ke, mgp.i);
  fflush (fp);
}

/**
At the end of the simulation, we use gnuplot to fit a function of the form
$$
k(t) = ae^{-bt}(1-\cos(ct))
$$
to the kinetic energy. This gives estimates of the oscillation
pulsation *c* and of the damping *b*.

We also compute the relative error on the pulsation, using the
theoretical value $\omega_0$ as reference. */

event fit (t = end) {
  FILE * fp = popen ("gnuplot 2>&1", "w");
  fprintf (fp, 
           "k(t)=a*exp(-b*t)*(1.-cos(c*t))\n"
           "a = 3e-4\n"
           "b = 1.5\n"
	   "\n"
           "D = %g\n"
           "n = 2.\n"
           "sigma = 1.\n"
           "rhol = 1.\n"
           "rhog = 1./1000.\n"
           "r0 = D/2.\n"
           "omega0 = sqrt((n**3-n)*sigma/((rhol+rhog)*r0**3))\n"
	   "\n"
           "c = 2.*omega0\n"
           "fit k(x) 'k-%d' via a,b,c\n"
	   "level = %d\n"
	   "res = D/%g*2.**level\n"
	   "print sprintf (\"fit %%g %%.6f %%.2f %%.0f\\n\", res, a, b, c, D)\n"
	   "\n"
	   "set table 'fit-%d'\n"
	   "plot [0:1] 2.*a*exp(-b*x)\n"
	   "unset table\n"
	   "\n"
	   "set print 'error' append\n"
	   "print res, c/2./omega0-1., D\n"
	   "\n"
	   "set print 'laplace' append\n"
	   "empirical_constant = 30.\n"
	   "print res, (1./(b**2.*D**3.))*empirical_constant**2, D\n"
	   "\n",
	   D, LEVEL, LEVEL, L0, LEVEL);
  pclose (fp);
}

#if TREE
event adapt (i++) {
#if MOMENTUM || COMPRESSIBLE
  vector u[];
  foreach()
    foreach_dimension()
      u.x[] = q.x[]/rho[];
#else
# if GFM
  scalar fa[];
  foreach() {
    if (interfacial_ghost (point, f))
      fa[] = rand();
    else
      fa[] = 0.;
  }
  adapt_wavelet ({fa,u}, (double[]){5e-3,1e-3,1e-3}, LEVEL);
# else
  adapt_wavelet ({f,u}, (double[]){5e-3,1e-3,1e-3}, LEVEL);
# endif
#endif
}
#endif

/**
## Results

~~~gnuplot Evolution of the kinetic energy as a function of time for the spatial resolutions (number of grid points per diameter) indicated in the legend. The black lines are fitted decreasing exponential functions.
set xlabel 'Time'
set ylabel 'Kinetic energy'
set logscale y
plot [0:1][8e-5:]'k-8' t "51.2" w l, 'k-7' t "25.6" w l,               \
  'k-6' t "12.8" w l, 'k-5' t "6.4" w l,			       \
  'fit-8' t "" w l lt 7, 'fit-7' t "" w l lt 7, 'fit-6' t "" w l lt 7, \
  'fit-5' t "" w l lt 7
~~~

~~~gnuplot Relative error in the oscillation frequency as a function of resolution.
set xlabel 'Diameter (grid points)'
set ylabel 'Frequency error (%)'
set logscale x 2
unset grid
set xzeroaxis
set key spacing 1.5 top right
ftitle(a,b,c) = sprintf("%.0f/x^{%4.2f} (%s)", exp(a), -b, c)
f(x)=a+b*x
fit f(x) 'error' u (log($1)):(log(abs($2)*100.)) via a,b
plot 'error' u ($1):(abs($2)*100.) t "" w p pt 5 ps 1, \
      exp(f(log(x))) t ftitle(a,b,"standard GFM")
~~~

The amount of numerical damping can be estimated by computing an
equivalent viscosity. With viscosity, kinetic energy is expected to
decrease as:
$$\exp(-C\nu/D^2t)$$
where $C$ is a constant, $\nu$ the viscosity and $D$ the droplet
diameter. Using curve fitting the damping coefficient $b=C\nu/D^2$
can be estimated (black curves on Figure \ref{kinetic}). An
equivalent Laplace number can then be computed as:
$$La=\frac{\sigma D}{\rho\nu^2}=\frac{\sigma C^2}{\rho b^2 D^3}$$
The equivalent Laplace number depends on spatial resolution as
illustrated below.

~~~gnuplot Equivalent Laplace number estimated from the numerical damping of kinetic energy.
set xlabel 'Diameter (grid points)'
set ylabel 'Equivalent Laplace number'
set grid
set key bottom right
plot 'laplace' t "standard GFM" w p pt 5 ps 1
~~~

## See also

* [Same test with Gerris](http://gerris.dalembert.upmc.fr/gerris/tests/tests/oscillation.html)
*/
