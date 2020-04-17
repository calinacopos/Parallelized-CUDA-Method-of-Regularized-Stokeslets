#ifndef _TEST_KERNEL_H_
#define _TEST_KERNEL_H_

#include <stdio.h>

// Thread block size
#define BLOCK_SIZE 4
#define width 2*2
#define height 2*2

// CUDA kernel
__global__ void TestKernel(double* xd, double* yd, double* Fxd, double* Fyd, double* Gxd, double* Gyd, double* Vxd, double* Vyd, int N, double visc, double e, double pi)
{
        // 2D thread ID
        int tx = blockIdx.x * blockDim.x + threadIdx.x;
        int ty = blockIdx.y * blockDim.y + threadIdx.y;

        // Each thread fills Stokeslet matrix 
        double rk = sqrt(powf(xd[tx]-xd[ty],2) + powf(yd[tx]-yd[ty],2));
        double sq = sqrtf(powf(rk,2) + powf(e,2));
        double p1 = (1.0f/(4.0f*visc*pi)) * (logf(sq+e)-(e*(sq+2*e))/(sq*(sq+e)));
        double p2 = (1.0f/(4.0f*visc*pi)) * (sq+2*e)/(sq*powf(sq+e,2));
        Gxd[tx+N*N*ty] = -1.0f*p1*Fxd[ty] + p2*(powf(xd[tx]-xd[ty],2)*Fxd[ty] + (xd[tx]-xd[ty])*(yd[tx]-yd[ty])*Fyd[ty]);
        Gyd[tx+N*N*ty] = -1.0f*p1*Fyd[ty] + p2*((xd[tx]-xd[ty])*(yd[tx]-yd[ty])*Fxd[ty] + powf(yd[tx]-yd[ty],2)*Fyd[ty]);
	//Gxd[tx+N*N*ty] = 0.0f; Gyd[tx+N*N*ty] = 0.0f;        

        // Synchronize all threads to ensure the matrix is computed and loaded
        __syncthreads();

	//Vxd[tx] = 2.0f; Vyd[tx] = 2.0f;

        // Sum up rows of stokeslet matrix
	double vx = 0.0f; double vy = 0.0f;
        if (ty == 0) {
                for(int k=0; k<N*N; k++) {
                   vx += Gxd[tx*N*N+k];
                   vy += Gyd[tx*N*N+k];
                }
                // Write the array to device memory
                Vxd[tx] = vx; Vyd[tx] = vy;
        } 
}
#endif // #ifndef _TEST_KERNEL_H

