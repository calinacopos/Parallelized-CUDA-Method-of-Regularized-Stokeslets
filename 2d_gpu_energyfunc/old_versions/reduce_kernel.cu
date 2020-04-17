

#ifndef _REDUCE_KERNEL_H_
#define _REDUCE_KERNEL_H_

// Thread block size
#define POINTS_F 2048
#define THREADS_F 1024
#define BLOCKS_F 8

__global__ void ReduceKernel(int dim, double* tmpx, double* tmpy, double* Vxd, double* Vyd) 
{
  // cache thread ID
  int tx = threadIdx.x;
  int idx = threadIdx.x + blockIdx.x * THREADS_F;

  for(int i=0; i<dim; i++) {
	Vxd[idx] += tmpx[idx*dim+i];
	Vyd[idx] += tmpy[idx*dim+i];
  }
}
#endif 
