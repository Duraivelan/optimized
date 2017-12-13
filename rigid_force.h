#ifndef RIGID_FORCE_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define RIGID_FORCE_H

# include "defs.h"

void forceUpdate( vector<SubData>& particle,  double *p_energy , int* combine_now , int combine[][4], long long int* step ,const int NrParticles, const double Lx,  const double Ly,  const double Lz, double STRAIN );

#endif
