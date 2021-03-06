#include <cuda_runtime.h>
#include "curand_kernel.h"
#include "selfPropelledCellVertexDynamics.cuh"

/** \file selfPropelledCellVertexDynamics.cu
    * Defines kernel callers and kernels for GPU calculations of simple active 2D cell models
*/

/*!
    \addtogroup simpleEquationOfMotionKernels
    @{
*/

/*!
  In this version of the active vertex model, the motility of a vertex is a straight average of the
  motility of the three adjacent cells
  */
__global__ void calculate_vertex_displacement_kernel(
                                        double2 *d_forces,
                                        double2 *d_displacements,
                                        double2 *motility,
                                        double  *d_cellDirectors,
                                        int      *d_vertexCellNeighbors,
                                        double  deltaT,
                                        double  mu,
                                        int      Nvertices)
    {
    unsigned int idx = blockDim.x * blockIdx.x + threadIdx.x;
    if (idx >= Nvertices)
        return;

    //the vertex motility is the average of the motility of the connected cells
    int vn1 = d_vertexCellNeighbors[3*idx];
    int vn2 = d_vertexCellNeighbors[3*idx+1];
    int vn3 = d_vertexCellNeighbors[3*idx+2];
    double v1 = motility[vn1].x;
    double v2 = motility[vn2].x;
    double v3 = motility[vn3].x;

    double directorx =
            (v1*Cos(d_cellDirectors[vn1])+v2*Cos(d_cellDirectors[vn2])+v3*Cos(d_cellDirectors[vn3]))/3.0;
    double directory =
            (v1*Sin(d_cellDirectors[vn1])+v2*Sin(d_cellDirectors[vn2])+v3*Sin(d_cellDirectors[vn3]))/3.0;
    //update positions from forces and motility
    d_displacements[idx].x = deltaT*(directorx + mu*d_forces[idx].x);
    d_displacements[idx].y = deltaT*(directory + mu*d_forces[idx].y);
    };

/*!
After the vertices have been moved, the directors of the cells have some noise.
  */
__global__ void rotate_directors_kernel(
                                        double  *d_cellDirectors,
                                        curandState *d_curandRNGs,
                                        double2 *motility,
                                        double  deltaT,
                                        int      Ncells)
    {
    unsigned int idx = blockDim.x * blockIdx.x + threadIdx.x;
    if (idx >= Ncells)
        return;

    //get the per-cell RNG, rotate the director, return the RNG
    curandState_t randState;
    randState=d_curandRNGs[idx];
    d_cellDirectors[idx] += cur_norm(&randState)*sqrt(2.0*deltaT*motility[idx].y);
    d_curandRNGs[idx] = randState;
    };



//!get the current timesteps vector of displacements into the displacement vector, rotate the cells
bool gpu_spp_cellVertex_eom_integration(
                    double2 *forces,
                    double2 *displacements,
                    double2 *motility,
                    double  *cellDirectors,
                    int      *vertexCellNeighbors,
                    curandState *RNGs,
                    int Nvertices,
                    int Ncells,
                    double deltaT,
                    int Timestep,
                    double mu)
    {
    unsigned int block_size = 128;
    if (Nvertices < 128) block_size = 32;
    unsigned int nblocks  = Nvertices/block_size + 1;

    //displace vertices
    calculate_vertex_displacement_kernel<<<nblocks,block_size>>>(forces,displacements,motility,
                                                         cellDirectors,vertexCellNeighbors,
                                                         deltaT,mu,Nvertices);
    HANDLE_ERROR(cudaGetLastError());

    //rotate cell directors
    if (Ncells < 128) block_size = 32;
    nblocks = Ncells/block_size + 1;
    rotate_directors_kernel<<<nblocks,block_size>>>(cellDirectors,RNGs,
                                                        motility,deltaT,Ncells);
    HANDLE_ERROR(cudaGetLastError());

    return cudaSuccess;
    };

/** @} */ //end of group declaration

