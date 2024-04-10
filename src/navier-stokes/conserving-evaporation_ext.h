/**
# Momentum-conserving advection of velocity

This file implements momentum-conserving VOF advection of the velocity
components for the [two-phase](/src/two-phase.h) Navier--Stokes
solver.

On trees, we first define refinement and restriction functions which
guarantee conservation of each component of the total momentum. Note
that these functions do not guarantee conservation of momentum for
each phase. */

vector q1ext[], q2ext[];

#if TREE
static void momentum_refine_ext (Point point, scalar u) {
  refine_bilinear (point, u);
  double rhou = 0.;
  foreach_child()
    //rhou += cm[]*rhov[]*u[];
    rhou += cm[]*rho(f[])*u[];
  double du = u[] - rhou/((1 << dimension)*(cm[] + SEPS)*rho(f[]));
  foreach_child()
    u[] += du;
}

static void momentum_restriction_ext (Point point, scalar u)
{
  double rhou = 0.;
  foreach_child()
    rhou += cm[]*rho(f[])*u[];
  u[] = rhou/((1 << dimension)*(cm[] + SEPS)*rho(f[]));
}
#endif // TREE

/**
We switch-off the default advection scheme of the [centered
solver](centered.h). */

event defaults (i = 0)
{
  stokes = true;

  foreach_dimension() {
    q1ext.x.nodump = true;
    q2ext.x.nodump = true;
  }

#if TREE

  /**
  On trees, the refinement and restriction functions above rely on the
  volume fraction field *f* being refined/restricted before the
  components of velocity. To ensure this, we move *f* to the front of
  the field list (*all*). */

  int i = 0;
  while (all[i].i != f.i) i++;
  while (i > 0 && all[i].i)
    all[i] = all[i-1], i--;
  all[i] = f;
    
  /**
  We then set the refinement and restriction functions for the
  components of the velocity field. The boundary conditions on
  $\mathbf{u}$ now depend on those on $f$. */
  
  foreach_dimension() {
    uext.x.refine = uext.x.prolongation = momentum_refine_ext;    // don't work well
    uext.x.restriction = momentum_restriction_ext;             // same
    uext.x.depends = list_add (uext.x.depends, f);
  }
#endif
}

/**
We need to overload the stability event so that the CFL is taken into
account (because we set stokes to true). */

event stability (i++)
  dtmax = timestep (uf, dtmax);

/**
We will transport the two components of the momentum, $q_1=f \rho_1
\mathbf{u}$ and $q_2=(1 - f) \rho_2 \mathbf{u}$. We will need to
impose boundary conditions which match this definition. This is done
using the functions below. */

foreach_dimension()
static double boundary_q1ext_x (Point neighbor, Point point, scalar q1ext, void * data)
{
  return clamp(f[],0.,1.)*rho1*uext.x[];
}

foreach_dimension()
static double boundary_q2ext_x (Point neighbor, Point point, scalar q2ext, void * data)
{
  return (1. - clamp(f[],0.,1.))*rho2*uext.x[];
}

/**
Similarly, on trees we need prolongation functions which also follow
this definition. */

#if TREE
foreach_dimension()
static void prolongation_q1ext_x (Point point, scalar q1ext) {
  foreach_child()
    q1ext[] = clamp(f[],0.,1.)*rho1*uext.x[];
}

foreach_dimension()
static void prolongation_q2ext_x (Point point, scalar q2ext) {
  foreach_child()
    q2ext[] = (1. - clamp(f[],0.,1.))*rho2*uext.x[];
}
#endif

static scalar * tracers1_ext = NULL;

/**
We overload the *vof()* event to transport consistently the volume
fraction and the momentum of each phase. */

event phasechange (i++) {

  /**
  We allocate two temporary vector fields to store the two components
  of the momentum and set the boundary conditions and prolongation
  functions. The boundary conditions on $q_1$ and $q_2$ depend on the
  boundary conditions on $f$. */
  
  for (scalar s in {q1ext,q2ext}) {
    s.depends = list_add (s.depends, f);
    foreach_dimension()
      s.v.x.i = -1; // not a vector
  }
  for (int i = 0; i < nboundary; i++)
    foreach_dimension() {
      q1ext.x.boundary[i] = boundary_q1ext_x;
      q2ext.x.boundary[i] = boundary_q2ext_x;
    }
#if TREE
  foreach_dimension() {
    q1ext.x.prolongation = prolongation_q1ext_x;
    q2ext.x.prolongation = prolongation_q2ext_x;
  }
#endif

  /**
  We split the total momentum $q$ into its two components $q1ext$ and
  $q2ext$ associated with $f$ and $1 - f$ respectively. */

  foreach()
    foreach_dimension() {
      double fc = clamp(f[],0,1);
      q1ext.x[] = fc*rho1*uext.x[];
      q2ext.x[] = (1. - fc)*rho2*uext.x[];
    }

  /**
  Momentum $q2ext$ is associated with $1 - f$, so we set the *inverse*
  attribute to *true*. We use the same slope-limiting as for the
  velocity field. */

  foreach_dimension() {
    q2ext.x.inverse = true;
    q1ext.x.gradient = q2ext.x.gradient = uext.x.gradient;
  }

  /**
  We associate the transport of $q1ext$ and $q2ext$ with $f$ and transport
  all fields consistently using the VOF scheme. */

  tracers1_ext = f.tracers;
  //tracers1 = f.tracers;
  f.tracers = list_concat (tracers1_ext, (scalar *){q1ext, q2ext});
  //f.tracers = list_concat (tracers1, (scalar *){q1ext, q2ext});

#ifdef VARPROP
  for (scalar s in {q1ext, q2ext})
    s.conservative = false;
#endif
}

event tracer_advection (i++) {

  /**
  We remove the momentum fields from the tracers lists (to avoid
  memory problems), and we restore the vof tracers list. */

  free (f.tracers);
  //f.tracers = tracers1;
  f.tracers = tracers1_ext;

  /**
  We recover the advected velocity field using the total momentum and
  the density */

  foreach()
    foreach_dimension()
      uext.x[] = (q1ext.x[] + q2ext.x[])/rho(f[]);
}

