#ifndef BACKGROUND_PROFILES_H
#define BACKGROUND_PROFILES_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "focal-struct.h"
#include "focal.h"
#include "hdf5.h"
#include "grid_io.h"

void init_background_profiles(  gridConfiguration *gridCfg,  
                                beamAntennaConfiguration *beamCfg,
                                double n_e[NX/2][NY/2][NZ/2], 
                                double J_B0[NX][NY][NZ] );

int make_density_profile(   gridConfiguration *gridCfg,
                            beamAntennaConfiguration *beamCfg, 
                            double cntrl_para,  
                            double n_e[NX/2][NY/2][NZ/2] );

int make_B0_profile( gridConfiguration *gridCfg,
                     double cntrl_para, 
                     double J_B0[NX][NY][NZ] );

int inside_cylinder(int x, int y, int ant_x, int ant_y, int radius );

#endif  // BACKGROUND_PROFILES_H
