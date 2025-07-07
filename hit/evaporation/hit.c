/**
   Set-up of HIT without and with droplets/bubbles */

/*
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
//#include "sandbox/maxruntime.h"   
//#include "sandbox/perfs.h"   
//#include "sandbox/my_funcs.h" 
*/

#define F_ERR 1e-10
#define NGS 2
#define NLS 1
//#define FILTERED
#define SOLVE_TEMPERATURE
//#define USE_GSL 0
#define USE_SUNDIALS
#define FSOLVE_ABSTOL 1.e-3
#define USE_CLAPEYRON

#include "../../../basilisk-sandbox-ecipriano/src/navier-stokes/centered-evaporation.h"
#include "../../../basilisk-sandbox-ecipriano/src/navier-stokes/centered-doubled.h"
#include "two-phase.h"
//#include "../../../basilisk-sandbox-ecipriano/src/navier-stokes/conserving-evaporation.h"
#include "tension.h"
#include "../../../basilisk-sandbox-ecipriano/src/evaporation.h"
#include "../../../basilisk-sandbox-ecipriano/src/multicomponent.h"
#include "tag.h"
#include "view.h"
#include "sandbox/maxruntime.h"   
#include "sandbox/perfs.h"   
#include "sandbox/my_funcs.h"  

/**
   The code takes the following input parameters (preinitialized to 1): */

double RHO_RATIO = 1.0;
double MU_RATIO  = 1.0;
double L0D0_RATIO = 1.0;
double D0 = 1.0;
double sigma = 1.0;
double amp_force = 1.0;
int FORCED = 1;
int MAX_LEVEL = 1;
int MIN_LEVEL = 1;
double end_sim = 1.0;
double t_phy = 1.0;
int single_ph = 1;
//double mEvapVal = -0.10;

/**
   Parameters for evaporation: */

char* gas_species[NGS] = {"H2O", "AIR"};
char* liq_species[NLS] = {"H2O"};
char* inert_species[1] = {"AIR"};
double gas_start[NGS] = {0., 1.};
double liq_start[NLS] = {1.};
double inDmix1[NLS] = {0.};
double inDmix2[NGS] = {2.075e-5, 2.075e-5};
double inKeq[NLS] = {0.0};

///*
double Tboil[NLS] = {373.15};
double lambda1 = 0.60700;
double lambda2 = 0.02624;
double dhev = 2.44e6;
double cp1 = 4137.00;
double cp2 = 1040.00;
double TL0 = 288.15;
double TG0 = 303.15;
//*/

/**
   Define some output frequencies. */

double t_out_stat = 64.0;  // output frequency of global observables 
double t_out_dpst = 16.0;  // output frequency for dump restarting files 
double t_out_dpbk = 4.00;  // output frequency back-up restarting files
double t_out_mov  = 16.0;  // output frequency of movie
double t_out_tag  = 16.0;  // output frequency of tagging

/**
   For the restarting step. */

int counter_max = 2;
int counter = 0;

int main (int argc, char * argv[]) {

  /**
     Provide the inputs */

  //maxruntime (&argc, argv);

  if (argc > 1)
    RHO_RATIO = atof (argv[1]);
  if (argc > 2)
    MU_RATIO = atof(argv[2]);
  if (argc > 3)
    L0D0_RATIO = atof(argv[3]);
  if (argc > 4)
    D0 = atof(argv[4]);
  if (argc > 5)
    sigma = atof(argv[5]);
  if (argc > 6)
    amp_force = atof(argv[6]);
  if (argc > 7)
    FORCED = atof(argv[7]);
  if (argc > 8)
    MAX_LEVEL = atoi(argv[8]);
  if (argc > 9)
    MIN_LEVEL = atoi(argv[9]);
  if (argc > 10)
    end_sim  = atof(argv[10]);
  if (argc > 11)
    t_phy = atof(argv[11]);
  if (argc > 12)
    single_ph = atof(argv[12]);

  if (argc < 13) {

    fprintf(ferr, "Lack of command line arguments. Check! Need %d more arguments\n", 13-argc);
    return 1;

  }

  fprintf(stderr, "************************\n"), fflush (stderr);
  fprintf(stderr, "maximum runtime = %.10e seconds\n", _maxruntime), fflush (stderr);
  fprintf(stderr, "Check of input parameters only\n"), fflush (stderr);
  fprintf(stderr, " RHO_RATIO = %.10e\n ", RHO_RATIO), fflush (stderr);
  fprintf(stderr, " MU_RATIO  = %.10e\n ", MU_RATIO), fflush (stderr);
  fprintf(stderr, " L0D0_ratio  = %.10e\n ", L0D0_RATIO), fflush (stderr);
  fprintf(stderr, " D0 = %.10e\n ", D0), fflush (stderr);
  fprintf(stderr, " sigma  = %.10e\n ", sigma), fflush (stderr);
  fprintf(stderr, " amp_force = %.10e\n ", amp_force), fflush (stderr);
  fprintf(stderr, " FORCED = %d\n ", FORCED), fflush (stderr);
  fprintf(stderr, " MAXLEVEL = %d\n ", MAX_LEVEL), fflush (stderr);
  fprintf(stderr, " MINLEVEL = %d\n ", MIN_LEVEL), fflush (stderr);
  fprintf(stderr, " tend  = %.10e\n ", end_sim), fflush (stderr);
  fprintf(stderr, " t_phy = %.10e\n ", t_phy), fflush (stderr);
  fprintf(stderr, " single_ph = %d\n ", single_ph), fflush (stderr);
  fprintf(stderr, "************************\n"), fflush (stderr);

  /**
     Set the domain size + boundary conditions + origin of the domain. */

  L0 = L0D0_RATIO*D0;
  size (L0);
  foreach_dimension() {
    periodic (right);
  }
  origin (0.0, 0.0, 0.0);

  /**
     Set the initial grid. */

#if TREE
  N = 1 << (MAX_LEVEL-2);
#else
  N = 1 << MAX_LEVEL;
#endif

  /**
     Reference properties = gas. 
     Liquid properties computed from the density/viscosity ratio. */

  rho2 = 1.225;
  mu2  = 1.78e-5;
  rho1 = rho2*RHO_RATIO;
  mu1  = mu2*MU_RATIO;

  /**
     Set the surface tension.*/

  f.sigma = sigma;
  
  /**
     Set tolerance on flow divergence + maximum number of iterations. */

  //TOLERANCE = 1.0e-4;
  NITERMAX  = 20; 
   
  /**
     We finally run. */
										
  run();

}
										
/**
## Initial conditions */

void vorticity3D (vector u, vector omg, scalar omg_mod) {

  foreach() {
    omg.x[] = ( (u.z[0,1,0]-u.z[0,-1,0]) - (u.y[0,0,1]-u.y[0,0,-1]) )/(2.0*Delta);
    omg.y[] = ( (u.x[0,0,1]-u.x[0,0,-1]) - (u.z[1,0,0]-u.z[-1,0,0]) )/(2.0*Delta);
    omg.z[] = ( (u.y[1,0,0]-u.y[-1,0,0]) - (u.x[0,1,0]-u.x[0,-1,0]) )/(2.0*Delta);
    omg_mod[] = 0.;
    foreach_dimension() {
      omg_mod[] += sq(omg.x[]);
    }
    omg_mod[] = sqrt(omg_mod[]);
  }

}

# define POPEN(name, mode) fopen (name ".ppm", mode)

event init (i = 0) {

  /**
     Print relevant simulation parameters */

  fprintf(stderr, "************************\n"), fflush (stderr);
  fprintf(stderr, "A-posteriori check of simulation parameters\n"), fflush (stderr);
  fprintf(stderr, " RHO_RATIO = %.10e\n ", rho1/rho2), fflush (stderr);
  fprintf(stderr, " MU_RATIO  = %.10e\n ", mu1/mu2), fflush (stderr);
  fprintf(stderr, " L0D0_RATIO  = %.10e\n ", L0/D0), fflush (stderr);
  fprintf(stderr, " D0 = %.10e\n ", D0), fflush (stderr);
  fprintf(stderr, " sigma  = %.10e\n ", sigma), fflush (stderr);
  fprintf(stderr, " amp_force = %.10e\n ", amp_force), fflush (stderr);
  fprintf(stderr, " FORCED = %d\n ", FORCED), fflush (stderr);
  fprintf(stderr, " MAXLEVEL = %d\n ", MAX_LEVEL), fflush (stderr);
  fprintf(stderr, " MINLEVEL = %d\n ", MIN_LEVEL), fflush (stderr);
  fprintf(stderr, " tend  = %.10e\n ", end_sim), fflush (stderr);
  fprintf(stderr, " t_phy = %.10e\n ", t_phy), fflush (stderr);
  fprintf(stderr, " single_ph = %d\n ", single_ph), fflush (stderr);
  fprintf(stderr, "************************\n"), fflush (stderr);

  /**
  We set the molecular weights of the chemial species
  involved in the simulation (by default inMW=1). */

  inMW[0] = 18.0; inMW[1] = 29.0;

  /**
  We set the reference pressure. */

  Pref = 101325.0;

  /**
     Create directories for better organization of the output */
  
  fprintf(stderr, "Create directories if needed\n"), fflush (stderr);

  char comm[80];
  sprintf (comm, "mkdir -p restart_bk");
  system(comm);
  sprintf (comm, "mkdir -p tagging");
  system(comm);

  /**
     Start a new simulation or restore a previous one */

  if (single_ph == 1) {
    fprintf(stderr, "Single-phase simulation\n"), fflush (stderr);
    if (!restore (file = "restart.bin")) {
      fprintf(stderr, "Start a new simulation to generate sustained HIT\n"), fflush (stderr);
     
      double waveno = L0/(2.0*pi);
      double amp = 6.0;
      
      foreach() {
        f[] = 0.; // only gas at the beginning
    
        double xw = x/waveno;
        double yw = y/waveno;
        double zw = z/waveno;
    
        u.x[] = amp*( cos(yw) + sin(zw) );
        u.y[] = amp*( sin(xw) + cos(zw) );
        u.z[] = amp*( cos(xw) + sin(yw) );
    
      }
    }
    else {
      fprintf(stderr, "We restart from a single-phase precursor\n"), fflush (stderr);
    }
  }
  else {
    fprintf(stderr, "Start a two-phase simulation\n"), fflush (stderr);
    if (restore (file = "restart_sp.bin")) { // a precursor
      fprintf(stderr, "Restart from statistically stationary single-phase HIT\n"), fflush (stderr);
      fprintf(stderr, "We add the droplet\n"), fflush (stderr);
     
      refine( sq(2.0*0.5*D0)-sq(x-L0/2.0)-sq(y-L0/2.0)-sq(z-L0/2.0)>0 && level < MAX_LEVEL );
      fraction (f, -(sq(x-L0/2.0) + sq(y-L0/2.0) + sq(z-L0/2.0) - sq(0.5*D0)));
      
      foreach() {      
        foreach_dimension() {
          u.x[] = (1.0-f[])*u.x[]; // we set to 0 the velocity inside the droplet
        }
      }

      // print an image for checking
      vector omg[];
      scalar omg_mod[];
      vorticity3D(u, omg, omg_mod);
      
      coord ubar;
      foreach_dimension() {
        stats s = statsf(u.x);
        ubar.x = s.sum/s.volume;
      }
      
      double ke = 0., vd = 0., vol = 0.;
      foreach(reduction(+:ke) reduction(+:vd) reduction(+:vol)) {
        vol += dv();
        foreach_dimension() {
          ke += dv()*sq(u.x[] - ubar.x);
          vd += dv()*(sq(u.x[1] - u.x[-1]) +
            	  sq(u.x[0,1] - u.x[0,-1]) +
            	  sq(u.x[0,0,1] - u.x[0,0,-1]))/sq(2.*Delta);
        }
      }
      ke /= 2.*vol;
      vd *= mu2/vol; 
      double vel_rms = sqrt(2.0/3.0*ke); 
      double lambda  = vel_rms*sqrt(15.*(mu2/rho2)*(1.0/vd));
     
      foreach() {
        omg_mod[] /= vel_rms/lambda;
      }
      boundary({omg_mod});
      boundary({u.x,u.y,u.z});

      clear();
      view (fov = 44, camera = "iso", ty = .2,
            width = 1200, height = 1200, bg = {1,1,1}, samples = 4);
      box(false, lc = {1,1,1}, lw = 0.1);
      draw_vof("f");
      squares ("u.x", linear = true, map = cool_warm, alpha = 1.e-6, n = {1, 0, 0});
      squares ("u.y", linear = true, map = cool_warm, alpha = 1.e-6, n = {0, 1, 0});
      squares ("u.z", linear = true, map = cool_warm, alpha = 1.e-6, n = {0, 0, 1});
 
      char namedump[200];
      sprintf (namedump, "./3D_tw_%09d.ppm", i);
      save(namedump);
      //return 1;

    }
    else {
      fprintf(stderr, "We restart from a two-phase simulation\n"), fflush (stderr);
      restore (file = "restart.bin");

      // print an image for checking
      vector omg[];
      scalar omg_mod[];
      vorticity3D(u, omg, omg_mod);
      
      coord ubar;
      foreach_dimension() {
        stats s = statsf(u.x);
        ubar.x = s.sum/s.volume;
      }
      
      double ke = 0., vd = 0., vol = 0.;
      foreach(reduction(+:ke) reduction(+:vd) reduction(+:vol)) {
        vol += dv();
        foreach_dimension() {
          ke += dv()*sq(u.x[] - ubar.x);
          vd += dv()*(sq(u.x[1] - u.x[-1]) +
            	  sq(u.x[0,1] - u.x[0,-1]) +
            	  sq(u.x[0,0,1] - u.x[0,0,-1]))/sq(2.*Delta);
        }
      }
      ke /= 2.*vol;
      vd *= mu2/vol; 
      double vel_rms = sqrt(2.0/3.0*ke); 
      double lambda  = vel_rms*sqrt(15.*(mu2/rho2)*(1.0/vd));
      
      foreach() {
        omg_mod[] /= vel_rms/lambda;
      }
      
      boundary({omg_mod});
      boundary({u.x,u.y,u.z});
      clear();
      view (fov = 44, camera = "iso", ty = .2,
            width = 1200, height = 1200, bg = {1,1,1}, samples = 4);
      box(false, lc = {1,1,1}, lw = 0.1);
      draw_vof("f");
      squares ("u.x", linear = true, map = cool_warm, alpha = 1.e-6, n = {1, 0, 0});
      squares ("u.y", linear = true, map = cool_warm, alpha = 1.e-6, n = {0, 1, 0});
      squares ("u.z", linear = true, map = cool_warm, alpha = 1.e-6, n = {0, 0, 1});
 
      char namedump[200];
      sprintf (namedump, "./3D_tw_eva_%09d.ppm", i);
      save(namedump);
     
      /* 
      char dname[100];
      sprintf (dname, "dump_eva_%d.bin", i);
      dump (dname);
      return 1;
      */

    }
  }

}

/**
   We modify the acceleration event to sustain turbulence. */

event acceleration (i++) {

  face vector av = a;

  if (FORCED) { 
    coord ubar;
    foreach_dimension() {
      stats s = statsf(u.x);
      ubar.x = s.sum/s.volume;
    }
    foreach_face() {
      av.x[] += (1.-(f[] + f[-1])/2.)*amp_force*((u.x[] + u.x[-1])/2. - ubar.x);
    }
  }

}

/**
## Outputs */

/**
   We log some useful information. */

event log_simulation (i += 10) {

  int min_level = +100, max_level = -100;
  foreach(reduction(min:min_level) reduction(max:max_level)) {
    max_level = max(max_level,level);
    min_level = min(min_level,level);
  }

  if (pid() == 0) { 

    fflush(stderr);
    char name_1[80];
    sprintf (name_1,"log_simulation.out");
    FILE * log_sim = fopen(name_1,"a");
    fprintf (log_sim, "%8E %8E %8E %d %d\n", t, 1.0*i, dt, max_level, min_level);
    fclose(log_sim);

  }

}

/**
   We log the evolution of viscous dissipation, kinetic energy, and
   microscale Reynolds number. */

//event logfile (t = 0.0; t <= end_sim; t += t_phy/t_out_stat) {
event logfile (i += 10) {
//event logfile (i += 1) {

  coord ubar;
  foreach_dimension() {
    stats s = statsf(u.x);
    ubar.x = s.sum/s.volume;
  }

  double ke = 0., vd = 0., vol = 0.;
  foreach(reduction(+:ke) reduction(+:vd) reduction(+:vol)) {
    vol += dv();
    foreach_dimension() {
      ke += dv()*sq(u.x[] - ubar.x);
      vd += dv()*(sq(u.x[1] - u.x[-1]) +
		  sq(u.x[0,1] - u.x[0,-1]) +
		  sq(u.x[0,0,1] - u.x[0,0,-1]))/sq(2.*Delta);
    }
  }
  ke /= 2.*vol;
  vd *= mu2/vol; 
  double vel_rms   = sqrt(2.0/3.0*ke); 
  double lambda    = vel_rms*sqrt(15.*(mu2/rho2)*(1.0/vd));
  double len_ls    = cube(vel_rms)/vd;
  double re_lambda = rho2*vel_rms*lambda/mu2;
  double re_ls     = rho2*vel_rms*len_ls/mu2;
  double we_est    = 2.0*rho2*pow(vd,2.0/3.0)*pow(D0,5.0/3.0)/sigma; // estimated weber based on the input sigma

  if (pid() == 0) {

    fflush(stderr);

    FILE * fd = fopen("stats.out","a");

    fprintf (fd, "%8E %09d %8E %8E %8E %8E %8E %8E %8E %8E\n",
                  t, i, vd, ke, vel_rms, vd/(3.0*sq(vel_rms)), lambda, re_lambda, re_ls, we_est);
    fclose (fd);
    fflush (fd);

  }

  // We define some useful scalars 
  scalar YInt_H2O = YGIntList[0];
  scalar Y_H2O    = YList[0];

  // Compute volumetric quantities 
  double tot_vol0 = (pi/6.0)*cube(D0);
  double vol_l = 0.0; 
  double vol_g = 0.0; 
  double Yg_vol_avg = 0.0; 
  double Tg_vol_avg = 0.0; 
  double Tl_vol_avg = 0.0; 
  foreach(reduction(+:vol_l) reduction(+:vol_g) reduction(+:Yg_vol_avg)
	  reduction(+:Tg_vol_avg) reduction(+:Tl_vol_avg)) {
    vol_l += (0.0+f[])*dv();
    vol_g += (1.0-f[])*dv();
    Yg_vol_avg += (1.0-f[])*Y_H2O[]*dv();
    Tg_vol_avg += (1.0-f[])*T[]*dv();
    Tl_vol_avg += (0.0+f[])*T[]*dv();
  }
  Yg_vol_avg /= vol_g;
  Tg_vol_avg /= vol_g;
  Tl_vol_avg /= vol_l;

  // Compute surface quantities
  double thr = 1.0e-6;
  double mass_flux = 0;
  double tot_area  = 0;
  double Tint_avg1 = 0;
  double Tint_avg2 = 0;
  double Yint_avg1 = 0;
  double Yint_avg2 = 0;
  foreach(reduction(+:mass_flux) reduction(+:tot_area)
	  reduction(+:Tint_avg1) reduction(+:Tint_avg2)
	  reduction(+:Yint_avg1) reduction(+:Yint_avg2)) {
    if (f[] > thr && f[] < 1.0-thr) {
      coord n = interface_normal (point, f), p;
      double alpha = plane_alpha(f[], n);
      double ar = pow(Delta, dimension - 1)*plane_area_center(n, alpha, &p);
      mass_flux += mEvapTot[]*ar;
      tot_area  += ar;
      Tint_avg1 += TInt[]*ar;
      Tint_avg2 += T[]*ar;
      Yint_avg1 += YInt_H2O[]*ar;
      Yint_avg2 += Y_H2O[]*ar;
    }
  }
  tot_area += 1.0e-14; // to avoid division by 0

  Tint_avg1 /= tot_area;
  Tint_avg2 /= tot_area;
  Yint_avg1 /= tot_area;
  Yint_avg2 /= tot_area;

  double D32 = 6.0*vol_l/tot_area;
  double Deq = pow(6.0*vol_l/pi, 1.0/3.0);
  double sh_v1 = (-mass_flux/(rho2*tot_area*(Yint_avg1-Yg_vol_avg)))*D0/inDmix2[0];
  double sh_v2 = (-mass_flux/(rho2*tot_area*(Yint_avg1-Yg_vol_avg)))*D32/inDmix2[0];
  double sh_v3 = (-mass_flux/(rho2*tot_area*(Yint_avg1-Yg_vol_avg)))*Deq/inDmix2[0];

  if (pid() == 0) {

    fflush(stderr);
    char name_1[80];
    sprintf (name_1,"int_qtn.out");
    FILE * log_sim = fopen(name_1,"a");
    fprintf (log_sim, "%8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E %8E\n", 
		       t, 1.0*i, tot_area, vol_l/tot_vol0, -mass_flux, 
		       Yg_vol_avg, vol_g, 
		       Tg_vol_avg, Tl_vol_avg,
		       Tint_avg1, Tint_avg2,
		       Yint_avg1, Yint_avg2,
		       sh_v1, sh_v2, sh_v3);
    fclose(log_sim);

  }

}

event tagging (t = 0.0; t <= end_sim; t += t_phy/t_out_tag) {

  if (single_ph == 0) {

    char name_1[100], name_2[100], name_3[100];
    
    sprintf (name_1, "./tagging/num_bubbles_droplets_f1.out");
    sprintf (name_2, "./tagging/droplets_f1_%09d.out", i);
    sprintf (name_3, "./tagging/bubbles_f1_%09d.out", i);
    countDropsBubble (name_1, name_2, name_3, i, f);

  }
  
}

/*
event dump_pre (i = 178875) {
  
  char dname[100];
  sprintf (dname, "dump_eva_%d.bin", i);
  dump (dname);

}
*/

/** 
   Dump every fraction of t_phy */

event dumpstep (t = 0.0; t <= end_sim; t += t_phy/t_out_dpst) {
  
  if(counter < counter_max) {
    counter++;
  }
  else {
    counter = 1;
  }

  char dname[100];
  sprintf (dname, "dump_%d.bin", counter);
  dump (dname);

  /** 
     Add a symbolic link, log restarting info and size of the bin */

  if (pid () == 0) {

    char comm_2[80];
    sprintf (comm_2, "ln -sf dump_%d.bin restart.bin", counter);
    system(comm_2);

    fflush(stderr);
    char name_0[80];
    sprintf (name_0, "dump_%d.bin", counter);
    FILE * fp = fopen(name_0, "r");
    fseek(fp, 0L, SEEK_END);
    long int res = ftell(fp);
    fclose(fp);

    fflush(stderr);
    char name_1[80];
    sprintf (name_1,"log_restart.out");
    FILE * log_sim = fopen(name_1,"a");
    fprintf (log_sim, "%.10e %.10e %.10e %.10e\n", t, 1.0*i, 1.0*counter, 1.0*res);
    fclose(log_sim);

  }

}

event dump_backup (t = 0.0; t <= end_sim; t += t_phy/t_out_dpbk) {
 
  char dname[100];
  sprintf (dname, "./restart_bk/dump_%09d.bin", i);
  dump (dname);

}

/**
   We use adaptivity. */

#if TREE
event adapt (i++) {

  scalar H2O = YList[0];
  double femax = 1.00e-2;
  double uemax = 0.20*normf(u.x).avg;
  double semax = 0.20*normf(H2O).avg;
  double temax = 0.20*normf(T).avg;
  adapt_wavelet_leave_interface ({H2O,T,u.x,u.y,u.z}, {f},
                                 (double[]){semax, temax, uemax, uemax, uemax, femax}, MAX_LEVEL, MIN_LEVEL, 0);

  /*
  adapt_wavelet ((scalar *){f, u}, 
		 (double[]){femax, uemax, uemax, uemax}, MAX_LEVEL, MIN_LEVEL);
  */


}
#endif

event movie (t = 0.0; t <= end_sim; t += t_phy/t_out_mov) {

  vector omg[];
  scalar omg_mod[];
  vorticity3D(u, omg, omg_mod);

  coord ubar;
  foreach_dimension() {
    stats s = statsf(u.x);
    ubar.x = s.sum/s.volume;
  }

  double ke = 0., vd = 0., vol = 0.;
  foreach(reduction(+:ke) reduction(+:vd) reduction(+:vol)) {
    vol += dv();
    foreach_dimension() {
      ke += dv()*sq(u.x[] - ubar.x);
      vd += dv()*(sq(u.x[1] - u.x[-1]) +
		  sq(u.x[0,1] - u.x[0,-1]) +
		  sq(u.x[0,0,1] - u.x[0,0,-1]))/sq(2.*Delta);
    }
  }
  ke /= 2.*vol;
  vd *= mu2/vol; 
  double vel_rms  = sqrt(2.0/3.0*ke); 
  double lambda   = vel_rms*sqrt(15.*(mu2/rho2)*(1.0/vd));
 
  foreach() {
    omg_mod[] /= vel_rms/lambda;
  }

  boundary({omg_mod});
  boundary({u.x,u.y,u.z});
  clear();
  view (fov = 44, camera = "iso", ty = .2,
	width = 1200, height = 1200, bg = {1,1,1}, samples = 4);
  box(false, lc = {1,1,1}, lw = 0.1);
  draw_vof("f");
  //squares ("omg_mod", linear = true, map = cool_warm, alpha = 1.e-6, n = {1, 0, 0});
  //squares ("omg_mod", linear = true, map = cool_warm, alpha = 1.e-6, n = {0, 1, 0});
  //squares ("omg_mod", linear = true, map = cool_warm, alpha = 1.e-6, n = {0, 0, 1});
  squares ("u.x", linear = true, map = cool_warm, alpha = 1.e-6, n = {1, 0, 0});
  squares ("u.y", linear = true, map = cool_warm, alpha = 1.e-6, n = {0, 1, 0});
  squares ("u.z", linear = true, map = cool_warm, alpha = 1.e-6, n = {0, 0, 1});
  {
    static FILE * fp = POPEN ("movie", "a");
    save (fp = fp);
    fflush (fp);
  }

}

/**
## End 

   We want to run up end_sim physical time. */

event end (t = end_sim) {

  dump ("end.bin");

  if ( pid() == 0 ) {

    char comm[80];
    sprintf (comm, "ln -sf end.bin restart.bin");
    system(comm);

  }

  return 1; // exit

}
