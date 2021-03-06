#include <cuda_runtime.h>
#include "curand_kernel.h"
#include "selfPropelledParticleDynamics.cuh"

/** \file selfPropelledParticleDynamics.cu
    * Defines kernel callers and kernels for GPU calculations of simple active 2D cell models
*/

/*!
    \addtogroup simpleEquationOfMotionKernels
    @{
*/

/*!
Each thread calculates the displacement of an individual cell
*/
__global__ void spp_eom_integration_kernel(double2 *forces,
                                           double2 *velocities,
                                           double2 *displacements,
                                           double2 *motility,
                                           double *cellDirectors,
                                           curandState *RNGs,
                                           int N,
                                           double deltaT,
                                           int Timestep,
                                           double mu)
    {
    unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
    if (idx >=N)
        return;
    //update displacements
    displacements[idx].x = deltaT*(velocities[idx].x + mu*forces[idx].x);
    displacements[idx].y = deltaT*(velocities[idx].y + mu*forces[idx].y);

    //next, get an appropriate random angle displacement
    curandState_t randState;
    randState=RNGs[idx];
    double v0 = motility[idx].x;
    double Dr = motility[idx].y;
    double angleDiff = cur_norm(&randState)*sqrt(2.0*deltaT*Dr);
    RNGs[idx] = randState;
    //update director and velocity vector
    double currentTheta = cellDirectors[idx];
    velocities[idx].x = v0 * Cos(cellDirectors[idx]);
    velocities[idx].y = v0 * Sin(cellDirectors[idx]);
    if(velocities[idx].y != 0. && velocities[idx].x != 0.)
        {
        currentTheta = atan2(velocities[idx].y,velocities[idx].x);
        };
    cellDirectors[idx] = currentTheta + angleDiff;


    return;
    };

//!get the current timesteps vector of displacements into the displacement vector
bool gpu_spp_eom_integration(
                    double2 *forces,
                    double2 *velocities,
                    double2 *displacements,
                    double2 *motility,
                    double *cellDirectors,
                    curandState *RNGs,
                    int N,
                    double deltaT,
                    int Timestep,
                    double mu)
    {
    unsigned int block_size = 128;
    if (N < 128) block_size = 32;
    unsigned int nblocks  = N/block_size + 1;


    spp_eom_integration_kernel<<<nblocks,block_size>>>(
                                forces,velocities,displacements,motility,cellDirectors,
                                RNGs,
                                N,deltaT,Timestep,mu);
    HANDLE_ERROR(cudaGetLastError());
    return cudaSuccess;
    };

/** @} */ //end of group declaration

