#ifndef __SIMPLE2DCELL_CUH__
#define __SIMPLE2DCELL_CUH__

#include "std_include.h"
#include <cuda_runtime.h>
#include "periodicBoundaries.h"

/*!
 \file Simple2DCell.cuh
A file providing an interface to the relevant cuda calls for the Simple2DCell class
*/

/** @defgroup Simple2DCellKernels Simple2DCell Kernels
 * @{
 * \brief CUDA kernels and callers for the Simple2DCell class

 One might think that a "computeGeometry" function should be here, but this function depends
 too much on whether the primary degrees of freedom are cells or vertices
 */

//!Move degrees of freedom according to a set of displacements, and put them back in the unit cell
bool gpu_move_degrees_of_freedom(double2 *d_points,
                    double2 *d_disp,
                    int N,
                    periodicBoundaries &Box
                    );

//!The same as the above, but scale the displacements by a scalar (i.e., x[i] += scale*disp[i]
bool gpu_move_degrees_of_freedom(double2 *d_points,
                    double2 *d_disp,
                    double  scale,
                    int N,
                    periodicBoundaries &Box
                    );

//!The same as the above, but with substrate interactions
bool gpu_move_degrees_of_freedom_substrate(double2 *d_points,
                    			   double2 *d_disp,
					   double2 *d_ap,
					   double *d_sg,
					   double *d_ct,
                        		   double dx,
					   double deltax,
                        		   double *d_tau,
                    			   int N,
                    			   periodicBoundaries &Box
                    			   );

//!A utility function; set all copmonents of an integer array to value
bool gpu_set_integer_array(int *d_array,
                           int value,
                           int N
                           );

/** @} */ //end of group declaration

#endif
