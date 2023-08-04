/**
# Map Region

Map a region of the domain with a Heaviside function.
This function sets to 1 the cells containing the interface
and the gas phase (the liquid phase if inverse).
The number of layers defines how many layers of cells close to
the interface, and belonging to the liquid phase, must
be included (0 by default).

* *H*: Heaviside function to be filled
* *f*: Volume fraction field
* *nl*: Number of layers to be included
*/

#ifndef MAPREGION_H
#define MAPREGION_H

#ifndef F_ERR
# define F_ERR 1.e-10
#endif

struct _MapRegion {
  scalar H;       // heaviside function
  scalar f;       // vof field (f = 1 if liquid)
  int nl;         // number of additional layers (default 0, optional 1 or 2)
  int inverse;    // the vof field if = 1 if gas (default false)
  int interface;  // map just a narrow band around the interface (default false)
};

void mapregion (struct _MapRegion p) {
  scalar fc[];
  scalar H = p.H, f = p.f;
  int nl = p.nl ? p.nl : 0.;
  bool inverse = p.inverse ? true : false;
  bool interface = p.interface ? true : false;

  foreach()
    fc[] = (!inverse) ? f[] : 1. - f[];

  foreach() {
    if (interface)
      H[] = (fc[] > F_ERR && fc[] < 1.-F_ERR) ? 1. : 0.;
    else
      H[] = (fc[] < 1.-F_ERR) ? 1. : 0.;
    if (fc[] > 1.-F_ERR && nl > 0) {
      bool lightup = false;
      foreach_neighbor(nl) {
        if (fc[] > F_ERR && fc[] < 1.-F_ERR) {
          lightup = true;
          break;
        }
      }
      H[] = lightup ? 1. : 0.;
    }
  }
}

#endif