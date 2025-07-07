/**
   We also want to count the drops and bubbles in the flow. */

void countDropsBubble (char * name_1, char * name_2, char * name_3, int istep, scalar c) {

  //scalar m1[]; // droplets
  //scalar m2[]; // bubbles
  scalar m1 = new scalar; // droplets
  scalar m2 = new scalar; // bubbles
  foreach() {
    m1[] = c[] > 0.5; // i.e. set m true if f[] is close to unity (droplets)
    m2[] = c[] < 0.5; // m true if f[] close to zero (bubbles)
  }
  int n1 = tag(m1);
  int n2 = tag(m2);

  /**
  Having counted the bubbles, now we find their size. This example
  is similar to the jet atomization problem. We are interested in
  the volumes and positions of each droplet/bubble.*/

  double v1[n1]; // droplet
  double b1x[n1], b1y[n1], b1z[n1];  // droplet
  double v2[n2]; // bubble
  double b2x[n2], b2y[n2], b2z[n2];  // bubble

  /**
  We initialize to zero */

  for (int j=0; j<n1; j++) {
    v1[j] = b1x[j] = b1y[j] = b1z[j] = 0.0;
  }
  for (int j=0; j<n2; j++) {
    v2[j] = b2x[j] = b2y[j] = b2z[j] = 0.0;
  }

  /**
  We proceed with calculation. */

  fprintf(stderr, "I enter the for loop for reduction\n");
  foreach(reduction(+:v1[:n1]) reduction(+:v2[:n2])
	  reduction(+:b1x[:n1]) reduction(+:b1y[:n1]) reduction(+:b1z[:n1])
	  reduction(+:b2x[:n2]) reduction(+:b2y[:n2]) reduction(+:b2z[:n2])) {

    // droplets
    if (m1[] > 0) {
      int j = m1[] - 1;
      v1[j]  += dv()*c[]; 
      b1x[j] += dv()*c[]*x; 
      b1y[j] += dv()*c[]*y; 
      b1z[j] += dv()*c[]*z; 
    }

    // bubbles
    if (m2[] > 0) {
      int j = m2[] - 1;
      v2[j]  += dv()*(1.0-c[]);
      b2x[j] += dv()*(1.0-c[])*x; 
      b2y[j] += dv()*(1.0-c[])*y; 
      b2z[j] += dv()*(1.0-c[])*z; 
    }

  }
  delete({m1,m2});
  fprintf(stderr, "I go out of the for loop for reduction\n");

  if (pid() == 0) {

    fflush(stderr);

    /**
    We first output the number of droplets and bubbles. */

    FILE * ftot = fopen(name_1,"a");
    fprintf (ftot, "%.10e %.10e %.10e %.10e\n", t, 1.0*istep, 1.0*(n1-1), 1.0*(n2-1)); // we remove the main region
    fclose(ftot);

    /**
    We output separately the volume and position of each droplet and bubble to file. */

    FILE * fdrop = fopen(name_2,"a");
    for (int j=0; j<n1; j++) {
      fprintf (fdrop, "%d %.10e %.10e %.10e %.10e\n",
                         j, v1[j], b1x[j]/v1[j], b1y[j]/v1[j], b1z[j]/v1[j]);
      //fprintf (fdrop, "%d %.10e\n", j, v1[j]);
    }
    fclose(fdrop);
    fflush(fdrop);

    FILE * fbubb = fopen(name_3,"a");
    for (int j=0; j<n2; j++) {
      fprintf (fbubb, "%d %.10e %.10e %.10e %.10e\n",
                        j, v2[j], b2x[j]/v2[j], b2y[j]/v2[j], b2z[j]/v2[j]);
      //fprintf (fbubb, "%d %.10e\n", j, v2[j]);
    }
    fclose(fbubb);
    fflush(fbubb);

  }

}

