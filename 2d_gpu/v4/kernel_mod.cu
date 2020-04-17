// Stokes Flow with Drag Component kernel
// Last updated: 12/21/12

#ifndef _SLIP_KERNEL_H_
#define _SLIP_KERNEL_H_

// Thread block size
#define BLOCK_SIZE 16
#define width (128*128)
#define height (128*128)

// CUDA kernel
__global__ void SlipKernel(double* xd, double* yd, double* Fxd, double* Fyd, double* Vxd, double* Vyd, int N, double visc, double e, double pi)
{ 
        // block ID
        int bx = blockIdx.x;
        int by = blockIdx.y; // should be a constant (= 1)
 
        // cache thread ID
        int tx = threadIdx.x;
        int ty = threadIdx.y; // should also be a constant (= 1)

        // Index of the first sub-array to be filled and reduced
        int begin = bx * width;

        // Index of the last sub-array to be filled and reduced
        int end = begin + width - 1;

        // Stride used to iterate through the sub-arrays
        int stride = BLOCK_SIZE * BLOCK_SIZE;

        // Loop through all sub-arrays
        for (int i = begin; i <= end; i += stride) {
		// Declaration of shared memory arrays
		__shared__ float cache_x[BLOCK_SIZE*BLOCK_SIZE];
		__shared__ float cache_y[BLOCK_SIZE*BLOCK_SIZE];        	

		// Each thread fills Stokeslet matrix 
        	double rk = sqrt(powf(xd[i+tx]-xd[bx],2) + powf(yd[i+tx]-yd[bx],2));
        	double sq = sqrtf(powf(rk,2) + powf(e,2));
        	double p1 = (1.0f/(4.0f*visc*pi)) * (logf(sq+e)-(e*(sq+2*e))/(sq*(sq+e)));
        	double p2 = (1.0f/(4.0f*visc*pi)) * (sq+2*e)/(sq*powf(sq+e,2));
        
        	// Sub-Stokeslet matrix
        	cache_x[tx] = -1.0f*p1*Fxd[bx] + p2*(powf(xd[i+tx]-xd[bx],2)*Fxd[bx] + (xd[i+tx]-xd[bx])*(yd[i+tx]-yd[bx])*Fyd[bx]);
        	cache_y[tx] = -1.0f*p1*Fyd[bx] + p2*((xd[i+tx]-xd[bx])*(yd[i+tx]-yd[bx])*Fxd[bx] + powf(yd[i+tx]-yd[bx],2)*Fyd[bx]);        

        	// Synchronize all threads in a block to ensure submatrix is computed and loaded
        	__syncthreads();

		// Reduction
		int j = blockDim.x/2;
		// only half the threads work (rest chill and go on for the ride)
		while (j != 0) {
		   if (tx < j) {
		      cache_x[tx] += cache_x[tx+j];
		      cache_y[tx] += cache_y[tx+j];
		   }
		   __syncthreads();
		   j /= 2;
                }
		
		// Update velocity per stride
		if (tx == 0) {
		   Vxd[bx] += cache_x[0];
		   Vyd[bx] += cache_y[0];
		}

		// Synchronize to make sure computation has completed before loading new sub-arrays
		__syncthreads();
	}
	
}
#endif // #ifndef _SLIP_KERNEL_H
