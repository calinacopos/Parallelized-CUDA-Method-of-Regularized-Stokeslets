// Stokes Flow with Drag Component kernel
// Last updated: 12/05/12

#ifndef _SLIP_KERNEL_H_
#define _SLIP_KERNEL_H_

#include <stdio.h>

// Thread block size
#define BLOCK_SIZE 16
#define width (16*16)
#define height (16*16)

const int threadsPerBlock = BLOCK_SIZE*BLOCK_SIZE;

// CUDA kernel
__global__ void SlipKernel(double* xd, double* yd, double* Fxd, double* Fyd, double* Vxd, double* Vyd, double* tmpx, double* tmpy, int N, double visc, double e, double pi)
{
        // shared memory arrays
        __shared__ double cache_x[BLOCK_SIZE][BLOCK_SIZE];
        __shared__ double cache_y[BLOCK_SIZE][BLOCK_SIZE];
        
        // cache thread ID
        int cx = threadIdx.x;
        int cy = threadIdx.y;

	// 2D thread ID
        int tx = blockIdx.x * blockDim.x + threadIdx.x;
        int ty = blockIdx.y * blockDim.y + threadIdx.y;

        // Each thread fills Stokeslet matrix 
        double rk = sqrt(powf(xd[tx]-xd[ty],2) + powf(yd[tx]-yd[ty],2));
        double sq = sqrtf(powf(rk,2) + powf(e,2));
        double p1 = (1.0f/(4.0f*visc*pi)) * (logf(sq+e)-(e*(sq+2*e))/(sq*(sq+e)));
        double p2 = (1.0f/(4.0f*visc*pi)) * (sq+2*e)/(sq*powf(sq+e,2));
        
        // Sub-Stokeslet matrix
        cache_x[cx][cy] = -1.0f*p1*Fxd[ty] + p2*(powf(xd[tx]-xd[ty],2)*Fxd[ty] + (xd[tx]-xd[ty])*(yd[tx]-yd[ty])*Fyd[ty]);
        cache_y[cx][cy] = -1.0f*p1*Fyd[ty] + p2*((xd[tx]-xd[ty])*(yd[tx]-yd[ty])*Fxd[ty] + powf(yd[tx]-yd[ty],2)*Fyd[ty]);        

        // Synchronize all threads in a block to ensure submatrix is computed and loaded
        __syncthreads();

        // Sum up rows of stokeslet matrix
	double vtx = 0.0f; double vty = 0.0f;
        if (cy == 0) {
		// only half the threads work (rest chill and go on for the ride)
		for(int k=0; k<BLOCK_SIZE; k++) {
		   vtx += cache_x[k][cx];
 		   vty += cache_y[k][cx];
		}
		// Write the partial results to global memory
 		tmpx[blockDim.x+(width/BLOCK_SIZE)*blockDim.y] = vtx;
		tmpy[blockDim.x+(width/BLOCK_SIZE)*blockDim.y] = vty;	
	}

	// Synchronize all results to ensure the temporary values are computed and loaded
	__syncthreads();

	double vfx = 0.0f; double vfy = 0.0f;
	if (ty == 0) {
		for(int k=0; k<(width/BLOCK_SIZE); k++) {
		   vfx += tmpx[blockDim.x*(width/BLOCK_SIZE) + k];
		   vfy += tmpy[blockDim.x*(width/BLOCK_SIZE) + k];
		}
		Vxd[tx] = vfx; 
		Vyd[tx] = vfy;
	}
}
#endif // #ifndef _SLIP_KERNEL_H
