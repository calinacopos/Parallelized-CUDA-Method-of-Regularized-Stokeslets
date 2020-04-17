// Stokes Flow with Drag Component kernel
// Last updated: 02/14/13

#ifndef _SLIP_KERNEL_H_
#define _SLIP_KERNEL_H_

// Thread block size
#define POINTS 256
#define MATRIX_SIZE 16
#define THREADS 256
#define BLOCKS 256

__constant__ double INVERSE_QUARTER_PI =  M_1_PI/4.0f;

// CUDA kernel
__global__ void SlipKernel(double* xd, double* yd, double* Fxd, double* Fyd, double* Vxd, double* Vyd, double visc, double e, double esq)
{ 
        // block ID
        int bx = blockIdx.x;
        // cache thread ID
        int tx = threadIdx.x;
	int idx = threadIdx.x + blockIdx.x * THREADS;

	// Declaration of shared memory arrays
	__shared__ float cache_x[THREADS];
	__shared__ float cache_y[THREADS];        	

	int at = floor((float)idx/POINTS); /* 4 right shifting is probably faster than division by 16 */
        int dueto = idx % POINTS;

	// Each thread fills Stokeslet matrix 
        double rk = sqrt(powf(xd[at]-xd[dueto],2) + powf(yd[at]-yd[dueto],2));
        double sq = sqrtf(powf(rk,2) + esq);
        double p1 = (INVERSE_QUARTER_PI/visc) * (logf(sq+e)-(e*(sq+2*e))/(sq*(sq+e)));
        double p2 = (INVERSE_QUARTER_PI/visc) * (sq+2*e)/(sq*powf(sq+e,2));

        // Sub-Stokeslet matrix
        cache_x[tx] = -p1*Fxd[dueto] + p2*(powf(xd[at]-xd[dueto],2)*Fxd[dueto] + (xd[at]-xd[dueto])*(yd[at]-yd[dueto])*Fyd[dueto]);
        cache_y[tx] = -p1*Fyd[dueto] + p2*((xd[at]-xd[dueto])*(yd[at]-yd[dueto])*Fxd[dueto] + powf(yd[at]-yd[dueto],2)*Fyd[dueto]);        
//#define TEST_REDUCTION
#ifdef TEST_REDUCTION
        cache_x[tx] = 1.0; 
        cache_y[tx] = 1.0;
#endif

        // Synchronize all threads in a block to ensure submatrix is computed and loaded
        __syncthreads();

	//printf("DBG: thrd:%d block:%d & stokeslet (%f, %f)\n", tx, bx, cache_x[tx], cache_y[tx]);

	// Reduction
	// only half the threads work (rest chill and go on for the ride)
        int j = blockDim.x/2;   // keeps track of active threads
        int k = MATRIX_SIZE/2;  // keeps track of which neighbor you add you value with & simulateounsly
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
	if ( (idx%POINTS == 0) ) {
		for(int i=0; i<POINTS/MATRIX_SIZE; i++) {
		    Vxd[idx/POINTS] += cache_x[tx+i*MATRIX_SIZE];
		    Vyd[idx/POINTS] += cache_y[tx+i*MATRIX_SIZE];
		}
	}
}
#endif // #ifndef _SLIP_KERNEL_H
