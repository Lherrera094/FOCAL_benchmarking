#ifndef GRID_IO_H
#define GRID_IO_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

#include "focal-struct.h"
#include "macros-grid.h"
#include "focal.h"
#include "hdf5.h"

//Write timetraces to file
int writeTimetraces2ascii( int T_end, double Period, 
                           char filename[], double **timetraces );

//#ifdef HDF5
int writeMyHDF_v4( int dim0, int dim1, int dim2, char filename[], char dataset[], double array_3D[dim0][dim1][dim2] );
int writeConfig2HDF( gridConfiguration *gridCfg, beamAntennaConfiguration *beamCfg, char filename[] );
int readMyHDF( int dim0, int dim1, int dim2, char filename[], char dataset[], double array_3D[dim0][dim1][dim2]);
//#endif

//#if defined(HDF5) && defined(DETECTOR_ANTENNA_1D)
int detAnt1D_write2hdf5(    int N_x, 
                            char filename[], char detAnt_groupName[], 
                            size_t detAnt_ypos, size_t detAnt_zpos,
                            double **detAnt_fields );

#endif  // GRID_IO_H
