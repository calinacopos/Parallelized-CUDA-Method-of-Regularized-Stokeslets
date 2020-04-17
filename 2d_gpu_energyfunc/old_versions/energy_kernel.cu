// Stokes Flow with Drag Component kernel
// Last updated: 02/14/13

#ifndef _SLIP_KERNEL_H_
#define _SLIP_KERNEL_H_

// Thread block size
#define MATRIX_SIZE 16
#define THREADS MATRIX_SIZE*MATRIX_SIZE
#define BLOCKS 16

#define BLOCK_AREA THREADS

__constant__ double INVERSE_QUARTER_PI =  M_1_PI/4.0f;
__constant__ double E = 0.046451612903224;
__constant__ double E_SQUARE= 0.002157752341310;

// CUDA kernel
__global__ void SlipKernel(double* xd, double* yd, double* Fxd, double* Fyd, double* Vxd, double* Vyd, double visc)
{ 
        // block ID
        int bx = blockIdx.x;
        // cache thread ID
        int tx = threadIdx.x;
	int idx = threadIdx.x + blockIdx.x * blockDim.x;

        // Index of the first sub-array to be filled and reduced
        int begin = bx * BLOCKS;
        // Index of the last sub-array to be filled and reduced
        int end = begin + BLOCKS - 1;
        // Stride used to iterate through the sub-arrays
        int stride = BLOCK_AREA;

	// Declaration of shared memory arrays
	__shared__ float cache_x[BLOCK_AREA];
	__shared__ float cache_y[BLOCK_AREA];        	

        // Loop through all sub-arrays
        for (int i = begin; i <= end; i += stride) {
		int at = tx >> 4; /* 4 right shifting is probably faster than division by 16 */
                int dueto = tx % 16;

		// Each thread fills Stokeslet matrix 
        	double rk = sqrt(powf(xd[at]-xd[dueto],2) + powf(yd[at]-yd[dueto],2));
        	double sq = sqrtf(powf(rk,2) + E_SQUARE);
        	double p1 = (INVERSE_QUARTER_PI/visc) * (logf(sq+E)-(E*(sq+2*E))/(sq*(sq+E)));
        	double p2 = (INVERSE_QUARTER_PI/visc) * (sq+2*E)/(sq*powf(sq+E,2));

        	// Sub-Stokeslet matrix
        	cache_x[tx] = -p1*Fxd[dueto] + p2*(powf(xd[at]-xd[dueto],2)*Fxd[dueto] + (xd[at]-xd[dueto])*(yd[at]-yd[dueto])*Fyd[dueto]);
        	cache_y[tx] = -p1*Fyd[dueto] + p2*((xd[at]-xd[dueto])*(yd[at]-yd[dueto])*Fxd[dueto] + powf(yd[at]-yd[dueto],2)*Fyd[dueto]);        

/* Uncomment/redefined for one of the testcase below */
#define TEST_REDUCTION2

#ifdef TEST_REDUCTION1
        	cache_x[tx] = (tx+tx/MATRIX_SIZE)%MATRIX_SIZE;
        	cache_y[tx] = (tx+tx/MATRIX_SIZE)%MATRIX_SIZE;
#endif
#ifdef TEST_REDUCTION2
        	cache_x[tx] = bx;
        	cache_y[tx] = bx;
#endif
#ifdef TEST_REDUCTION3
		// should be numbered 0 ... 15 on each row in each block
		cache_x[tx] = tx%MATRIX_SIZE;
		cache_y[tx] = tx%MATRIX_SIZE;
#endif

        	// Synchronize all threads in a block to ensure submatrix is computed and loaded
        	__syncthreads();

		//printf("DBG: thrd:%d blk:%d & stokeslet (%f, %f)\n", tx, bx, cache_x[tx], cache_y[tx]);

		// Reduction
		// only half the threads work (rest chill and go on for the ride)
		int j = blockDim.x/2;	// keeps track of active threads
		int k = MATRIX_SIZE/2;	// keeps track of which neighbor you add you value with & simulateounsly
					// many entries per row should be changed by this code
		while (j >= MATRIX_SIZE ) {
		   if ( (tx%MATRIX_SIZE) < k ) { // for each row we add your value + value of k away neighbor
		      cache_x[tx] = cache_x[tx] + cache_x[tx+k];
		      cache_y[tx] = cache_y[tx] + cache_y[tx+k];
		   }
		   j = j >> 1;
		   k = k >> 1;
		   __syncthreads();
                }

#if 0
for (i=0;i<MATRIX_SIZE;i++)
	printf("[%d] %d: %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", tx, i, cache_x[i*MATRIX_SIZE],cache_x[i*MATRIX_SIZE+1],cache_x[i*MATRIX_SIZE+2],cache_x[i*MATRIX_SIZE+3],cache_x[i*MATRIX_SIZE+4],cache_x[i*MATRIX_SIZE+5],cache_x[i*MATRIX_SIZE+6],cache_x[i*MATRIX_SIZE+7],cache_x[i*MATRIX_SIZE+8],cache_x[i*MATRIX_SIZE+9],cache_x[i*MATRIX_SIZE+10],cache_x[i*MATRIX_SIZE+11],cache_x[i*MATRIX_SIZE+12],cache_x[i*MATRIX_SIZE+13],cache_x[i*MATRIX_SIZE+14],cache_x[i*MATRIX_SIZE+15]);
#endif

		// Update velocity per stride
		// To write back we need 2 items:
		// index going from 0 to 15 for vectors Vxd, Vyd
		// based on the reduction above, entries 0, 16, ..., 240 have summed up the rows
		// only these entries should write back their results (notice they are divisible by 16)
		if ( tx%MATRIX_SIZE == 0 ) {
		   int index = tx/MATRIX_SIZE;
		   Vxd[index] = cache_x[tx];
		   Vyd[index] = cache_y[tx];

#ifdef TEST_REDUCTION2
		   if (cache_x[tx] != bx*BLOCKS) 
			printf("FAILED2: [%d:%d] %f\n", bx,tx,cache_x[tx]);
#endif
		}
	}
}
#endif // #ifndef _SLIP_KERNEL_H
