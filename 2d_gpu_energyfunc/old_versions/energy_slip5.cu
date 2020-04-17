/* slip.cu
 * GPU Benchmark Immersed Boundary Unstructured Grid
 * Based on "The Method of Regularized Stokelets" by R.Cortez`
 * Elastic force computed using energy-based formulation by Devendran + Peskin
 * C.Copos 02/21/2012
 */

/* WHICH VERSION IS THIS? */
/* - velocity results are computed per thread block (256) exclusively in shared memory
 * - works only up to 2048 mesh points
 * - block number represents which node we compute velocity for
 * - thread number (tx) within a block represents who contributes to velocity calculation
 * - check against cpu & perfect agreement statically & dynamically for x-coord (2/21/2013)
 * - compared to v4 this version launches another kernel to do final reduction
 * - no time gain, though!
 */

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <cstdio>
#include <cstdlib>
#include <stdio.h>
#include <cuda.h>

#include "test_kernel5.cu"
#include "reduce_kernel.cu"

using namespace std;

const double visc = 1.0f;
const double drag = 1.0f;
const float def = 0.1;

// Lame constants
const double lambda = 1.0;
const double mu = 1.0; // shear modulus (G)

// time 
const double TMAX = 0.1f;
const double tstep = 0.0001f;/*0.0001f;*/
static int tstop = floor(TMAX/tstep);
//static int tstop = 1;

// curve constants
const int N = 2048;
const double ri = 1.2f; 
const double ds = ri/(N-1); 
const double e = 1.2f*ds; // e: parameter determining width of blobs or cutoffs

// vector structure
typedef struct vector{
  double x; // x-component
  double y; // y-component
} vector;

// 2x2 matrix
typedef struct matrix{
  double x1; // term (1,1)
  double y1; // term (1,2)
  double x2; // term (2,1)
  double y2; // term (2,2)
} matrix;

// vertex 
typedef struct vertex{
  vector ref; // reference coords
  vector def; // deformation coords 
  vector force; // force
  int exterior; // 1 if this is a boundary point and 0 if this is an interior point
} vertex;

// triangle
typedef struct triangle{
  int A;
  int B;
  int C;
  double def_area; // deformed area
  double ref_area; // reformed area
  matrix f1; // term 1 of the forcing calculation
  matrix f2_0_x; // x-component of term 2 of the forcing calculation for vertex 0 
  matrix f2_0_y; // y-component of term 2 of the forcing calculation for vertex 0  
  matrix f2_1_x; // x-component of term 2 of the forcing calculation for vertex 1
  matrix f2_1_y; // y-component of term 2 of the forcing calculation for vertex 1 
  matrix f2_2_x; // x-component of term 2 of the forcing calculation for vertex 2
  matrix f2_2_y; // y-component of term 2 of the forcing calculation for vertex 2
  double f3; // term 3 of the forcing calculation
} triangle;

// gpu timing
double gpuTime = 0.0f;

// Compute time difference in seconds
double diffclock(clock_t s, clock_t e) {
   double diffticks = s-e;
   double diffms = (diffticks)/CLOCKS_PER_SEC;
   return diffms;
}

// Set up preliminary info per triangle (i.e. reference area or term 3 and term 2) 
void ref_info(int Npts, triangle &tr, vertex nodes[]) {
  // term 3 (otherwise known as reference area)
  tr.ref_area = 0.5*fabs( (nodes[tr.B].ref.x-nodes[tr.A].ref.x)*(nodes[tr.C].ref.y-nodes[tr.A].ref.y) - (nodes[tr.C].ref.x-nodes[tr.A].ref.x)*(nodes[tr.B].ref.y-nodes[tr.A].ref.y) );
  tr.f3 = tr.ref_area;

  // determinant of S
  double detS;
  detS = (1.0)*((nodes[tr.B].ref.x-nodes[tr.A].ref.x)*(nodes[tr.C].ref.y-nodes[tr.A].ref.y) - (nodes[tr.C].ref.x-nodes[tr.A].ref.x)*(nodes[tr.B].ref.y-nodes[tr.A].ref.y));

  // term 2
  tr.f2_0_x.x1 = (1.0/detS)*(-1.0*nodes[tr.C].ref.y + nodes[tr.B].ref.y);
  tr.f2_0_x.y1 = (1.0/detS)*(nodes[tr.C].ref.x - nodes[tr.B].ref.x);
  tr.f2_0_x.x2 = 0.0;
  tr.f2_0_x.y2 = 0.0;

  tr.f2_0_y.x1 = 0.0;
  tr.f2_0_y.y1 = 0.0;
  tr.f2_0_y.x2 = (1.0/detS)*(-1.0*nodes[tr.C].ref.y + nodes[tr.B].ref.y);
  tr.f2_0_y.y2 = (1.0/detS)*(nodes[tr.C].ref.x - nodes[tr.B].ref.x);

  tr.f2_1_x.x1 = (1.0/detS)*(nodes[tr.C].ref.y - nodes[tr.A].ref.y);
  tr.f2_1_x.y1 = (1.0/detS)*(nodes[tr.A].ref.x - nodes[tr.C].ref.x);
  tr.f2_1_x.x2 = 0.0;
  tr.f2_1_x.y2 = 0.0;

  tr.f2_1_y.x1 = 0.0;
  tr.f2_1_y.y1 = 0.0;
  tr.f2_1_y.x2 = (1.0/detS)*(nodes[tr.C].ref.y - nodes[tr.A].ref.y);
  tr.f2_1_y.y2 = (1.0/detS)*(nodes[tr.A].ref.x - nodes[tr.C].ref.x);

  tr.f2_2_x.x1 = (1.0/detS)*(nodes[tr.A].ref.y - nodes[tr.B].ref.y);
  tr.f2_2_x.y1 = (1.0/detS)*(nodes[tr.B].ref.x - nodes[tr.A].ref.x);
  tr.f2_2_x.x2 = 0.0;
  tr.f2_2_x.y2 = 0.0;

  tr.f2_2_y.x1 = 0.0;
  tr.f2_2_y.y1 = 0.0;
  tr.f2_2_y.x2 = (1.0/detS)*(nodes[tr.A].ref.y - nodes[tr.B].ref.y);
  tr.f2_2_y.y2 = (1.0/detS)*(nodes[tr.B].ref.x - nodes[tr.A].ref.x);
}

// Set up deformation specific info per triangle (i.e. deformed area and term 1)
void def_info(int Npts, triangle &tr, vertex nodes[]) {
  // deformed area
  tr.def_area = 0.5*fabs((nodes[tr.B].def.x-nodes[tr.A].def.x)*(nodes[tr.C].def.y-nodes[tr.A].def.y) - (nodes[tr.B].def.y-nodes[tr.A].def.y)*(nodes[tr.C].def.x-nodes[tr.A].def.x) );

  // deformation gradient tensor
  matrix a;
  double detS;
  detS = (1.0)*((nodes[tr.B].ref.x-nodes[tr.A].ref.x)*(nodes[tr.C].ref.y-nodes[tr.A].ref.y) - (nodes[tr.C].ref.x-nodes[tr.A].ref.x)*(nodes[tr.B].ref.y-nodes[tr.A].ref.y));

  a.x1 = (1.0/detS)*( (nodes[tr.B].def.x-nodes[tr.A].def.x)*(nodes[tr.C].ref.y-nodes[tr.A].ref.y) + (nodes[tr.C].def.x-nodes[tr.A].def.x)*(nodes[tr.A].ref.y - nodes[tr.B].ref.y) );
  a.y1 = (1.0/detS)*( (nodes[tr.B].def.x-nodes[tr.A].def.x)*(nodes[tr.A].ref.x-nodes[tr.C].ref.x) + (nodes[tr.C].def.x-nodes[tr.A].def.x)*(nodes[tr.B].ref.x-nodes[tr.A].ref.x) );
  a.x2 = (1.0/detS)*( (nodes[tr.B].def.y-nodes[tr.A].def.y)*(nodes[tr.C].ref.y-nodes[tr.A].ref.y) + (nodes[tr.C].def.y-nodes[tr.A].def.y)*(nodes[tr.A].ref.y-nodes[tr.B].ref.y) ); 
  a.y2 = (1.0/detS)*( (nodes[tr.B].def.y-nodes[tr.A].def.y)*(nodes[tr.A].ref.x - nodes[tr.C].ref.x) + (nodes[tr.C].def.y-nodes[tr.A].def.y)*(nodes[tr.B].ref.x - nodes[tr.A].ref.x) ); 

  // inverse transpose of deformation gradient tensor (w/o outside normalizers i.e. determinants) 
  matrix ait;
  ait.x1 = a.y2;
  ait.y1 = (-1.0)*(a.x2);
  ait.x2 = (-1.0)*(a.y1);
  ait.y2 = a.x1;

  // Cauchy stress tensor
  matrix sigma;
  // Material displacement gradient tensor ( = deformation gradient tensor - I)
  matrix d; 
  d.x1 = a.x1 - 1.0;
  d.y1 = a.y1;
  d.x2 = a.x2;
  d.y2 = a.y2 - 1.0;
  sigma.x1 = lambda*(d.x1+d.y2) + 2.0*mu*d.x1;
  sigma.y1 = mu*(d.y1+d.x2);
  sigma.x2 = mu*(d.x2+d.y1);
  sigma.y2 = lambda*(d.x1+d.y2) + 2.0*mu*d.y2;

  // term 1 (otherwise known as 1st Piola-Kirchhoff tensor)
  tr.f1.x1 = ( sigma.x1*ait.x1 + sigma.y1*ait.x2 );
  tr.f1.y1 = ( sigma.x1*ait.y1 + sigma.y1*ait.y2 );
  tr.f1.x2 = ( sigma.x2*ait.x1 + sigma.y2*ait.x2 );
  tr.f1.y2 = ( sigma.x2*ait.y1 + sigma.y2*ait.y2 );
}

// Compute velocity vector for all points in the grid
void velocity(int Npts, int Ntris, vertex Nodes[], vector f[], double velx[], double vely[]) {
  int mem_size = Npts;  
  int dim = Npts/256; // number of thread blocks that make up a vector computation 

  // Allocate host memory for result (velocity) 
  // THIS IS UNNECESSARY & I SHOULD CHANGE THIS
  double *vx = (double*) malloc(mem_size*dim*sizeof(double));
  double *vy = (double*) malloc(mem_size*dim*sizeof(double));

  // Allocate and fill host memory for force
  double *fxh = (double*) malloc(mem_size*sizeof(double));
  double *fyh = (double*) malloc(mem_size*sizeof(double));
  for(int j=0; j<Npts; j++) {
	fxh[j] = f[j].x;
	fyh[j] = f[j].y;
  }

  // Allocate and fill host memory for position
  double *xh = (double*) malloc(mem_size*sizeof(double));
  double *yh = (double*) malloc(mem_size*sizeof(double));
  for(int j=0; j<Npts; j++) {  
  	xh[j] = Nodes[j].def.x; 
        yh[j] = Nodes[j].def.y;
  }

  // Allocate device memory for x, y, F, v, and G (where G is the Stokeslet matrix) 
  double *xd, *yd, *Fxd, *Fyd, *tmpx, *tmpy;
  cudaMalloc((void**) &xd, mem_size*sizeof(double));
  cudaMalloc((void**) &yd, mem_size*sizeof(double));
  cudaMalloc((void**) &Fxd, mem_size*sizeof(double));
  cudaMalloc((void**) &Fyd, mem_size*sizeof(double));
  cudaMalloc((void**) &tmpx, mem_size*dim*sizeof(double));
  cudaMalloc((void**) &tmpy, mem_size*dim*sizeof(double));

  // Initialize device memory to zero
  cudaMemset(tmpx, 0x0, mem_size*dim*sizeof(double));
  cudaMemset(tmpy, 0x0, mem_size*dim*sizeof(double));
  
  // Copy position and force arrays to allocated device memory locations
  cudaMemcpy(xd, xh, mem_size*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(yd, yh, mem_size*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(Fxd, fxh, mem_size*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(Fyd, fyh, mem_size*sizeof(double), cudaMemcpyHostToDevice);

  cudaEvent_t start, stop;
  float elapsedTime; 
  cudaEventCreate(&start);
  cudaEventRecord(start, 0);
 
  // Perform Stokeslet computation
  dim3 threads(THREADS, 1);
  dim3 grid(BLOCKS, 1);

  double esq = e*e;
  //printf("Number of threads per block: %d and number of blocks per grid: %d\n", THREADS, BLOCKS);
  SlipKernel<<< grid, threads >>>(xd, yd, Fxd, Fyd, tmpx, tmpy, visc, e, esq);

  cudaEventCreate(&stop);
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);

  // Report timing
  cudaEventElapsedTime(&elapsedTime, start, stop);
  gpuTime += elapsedTime;
 
  // Copy the result from device to host
  cudaMemcpy(vx, tmpx, mem_size*dim*sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(vy, tmpy, mem_size*dim*sizeof(double), cudaMemcpyDeviceToHost);

  /*for(int i=0; i<mem_size*dim; i++) {
  	printf("(vx, vy) = (%.16f,%.16f)\n", vx[i], vy[i]);
  }*/

  double *Vxd, *Vyd;
  cudaMalloc((void**) &Vxd, mem_size*sizeof(double));
  cudaMalloc((void**) &Vyd, mem_size*sizeof(double));
  cudaMemset(Vxd, 0x0, mem_size*sizeof(double));
  cudaMemset(Vyd, 0x0, mem_size*sizeof(double)); 

  cudaEvent_t start_f, stop_f;
  cudaEventCreate(&start_f);
  cudaEventRecord(start_f, 0);
  elapsedTime = 0.0f;

  // Perform final reduction
  dim3 threads_f(THREADS_F, 1);
  dim3 grid_f(BLOCKS_F, 1);

  ReduceKernel<<< grid_f, threads_f >>>(dim, tmpx, tmpy, Vxd, Vyd);

  cudaEventCreate(&stop_f);
  cudaEventRecord(stop_f, 0);
  cudaEventSynchronize(stop_f);

  cudaMemcpy(velx, Vxd, mem_size*sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(vely, Vyd, mem_size*sizeof(double), cudaMemcpyDeviceToHost);

  // Report timing
  cudaEventElapsedTime(&elapsedTime, start_f, stop_f);
  gpuTime += elapsedTime;

  /*for(int i=0; i<mem_size; i++) {
	printf("(vx, vy) = (%.16f, %.16f)\n", velx[i], vely[i]);
  }*/

  // Set velocity
  /*for(int j=0; j<Npts*dim; j+=dim) {
	int index = floor(j/dim);
	v[index].x = 0.0;
	v[index].y = 0.0;
	for(int i=0; i<dim; i++) {
	   v[index].x += vx[j+i];
	   v[index].y += vy[j+i];
	}
 	printf("final vel: (%.16f, %.16f)\n", v[index].x, v[index].y);
  }*/

  // Clean up
  free(xh); free(yh);
  free(fxh); free(fyh);
  free(vx); free(vy);
  cudaFree(xd); cudaFree(yd); 
  cudaFree(Fxd); cudaFree(Fyd); 
  cudaFree(tmpx); cudaFree(tmpy);
  cudaFree(Vxd); cudaFree(Vyd);
}

// Progression
void progress(int Npts, int Ntris, vertex Nodes[], triangle Tris[]) {
  vector pos_init[Npts];
  double vx[Npts];
  double vy[Npts];
  vector v[Npts];
  vector f[Npts];
  double ftime = 0.0f;

  // file handling
  ofstream f1, f2, f3;
  f1.open("initial_pos_conf.txt"); f2.open("final_pos_conf.txt"); f3.open("ref_pos_conf.txt");
 
  // print initial configuration (i.e. with initial deformation as described in startCurve() )
  for(int i=0; i<Npts; i++) {
        // zero the force
        Nodes[i].force.x = 0.0;
        Nodes[i].force.y = 0.0;
	pos_init[i].x = Nodes[i].def.x;
        pos_init[i].y = Nodes[i].def.y;
        f1 << pos_init[i].x << " " <<  pos_init[i].y << endl;
        f3 << Nodes[i].ref.x << " " << Nodes[i].ref.y << endl;
  }
  f1.close();
  f3.close();

  for(int t=0; t<tstop; t++) {
     clock_t fbegin = clock();
     float ref_Tarea = 0.0; float def_Tarea = 0.0;

     // CYCLE THROUGH TRIANGLES AND COMPUTE FORCES
     for(int j=0; j<Ntris; j++) {
        //printf("making a call for triangle: (%d, %d, %d)\n", Tris[j].A, Tris[j].B, Tris[j].C);
        def_info(Npts, Tris[j], Nodes);
        // vertex A     
        Nodes[Tris[j].A].force.x += -1.0*(Tris[j].ref_area)*(Tris[j].f1.x1*Tris[j].f2_0_x.x1 + Tris[j].f1.y1*Tris[j].f2_0_x.y1 + Tris[j].f1.x2*Tris[j].f2_0_x.x2 + Tris[j].f1.y2*Tris[j].f2_0_x.y2);
        Nodes[Tris[j].A].force.y += -1.0*(Tris[j].ref_area)*(Tris[j].f1.x1*Tris[j].f2_0_y.x1 + Tris[j].f1.y1*Tris[j].f2_0_y.y1 + Tris[j].f1.x2*Tris[j].f2_0_y.x2 + Tris[j].f1.y2*Tris[j].f2_0_y.y2);
	
	// vertex B
        Nodes[Tris[j].B].force.x += -1.0*(Tris[j].ref_area)*(Tris[j].f1.x1*Tris[j].f2_1_x.x1 + Tris[j].f1.y1*Tris[j].f2_1_x.y1 + Tris[j].f1.x2*Tris[j].f2_1_x.x2 + Tris[j].f1.y2*Tris[j].f2_1_x.y2);
        Nodes[Tris[j].B].force.y += -1.0*(Tris[j].ref_area)*(Tris[j].f1.x1*Tris[j].f2_1_y.x1 + Tris[j].f1.y1*Tris[j].f2_1_y.y1 + Tris[j].f1.x2*Tris[j].f2_1_y.x2 + Tris[j].f1.y2*Tris[j].f2_1_y.y2);

	// vertex C
        Nodes[Tris[j].C].force.x += -1.0*(Tris[j].ref_area)*(Tris[j].f1.x1*Tris[j].f2_2_x.x1 + Tris[j].f1.y1*Tris[j].f2_2_x.y1 + Tris[j].f1.x2*Tris[j].f2_2_x.x2 + Tris[j].f1.y2*Tris[j].f2_2_x.y2);
        Nodes[Tris[j].C].force.y += -1.0*(Tris[j].ref_area)*(Tris[j].f1.x1*Tris[j].f2_2_y.x1 + Tris[j].f1.y1*Tris[j].f2_2_y.y1 + Tris[j].f1.x2*Tris[j].f2_2_y.x2 + Tris[j].f1.y2*Tris[j].f2_2_y.y2);

	ref_Tarea += Tris[j].ref_area;
        def_Tarea += Tris[j].def_area;
     }
     clock_t fend = clock();

     for(int k=0; k<Npts; k++) {
        f[k].x = Nodes[k].force.x;
        f[k].y = Nodes[k].force.y;
        Nodes[k].force.x = 0.0;
        Nodes[k].force.y = 0.0;
     }

     // compute velocity fields
     velocity(Npts, Ntris, Nodes, f, vx, vy);

     // for each node in unstructured mesh
     for(int i=0; i<Npts; i++) {
        v[i].x = vx[i] + f[i].x/drag;
        v[i].y = vy[i] + f[i].y/drag;

        //printf("%f %f %f %f %f %f\n", Nodes[i].def.x, Nodes[i].def.y, v[i].x, v[i].y, f[i].x, f[i].y);

        Nodes[i].def.x = Nodes[i].def.x + tstep*v[i].x;
        Nodes[i].def.y = Nodes[i].def.y + tstep*v[i].y;
        if(t==tstop-1) { f2 << Nodes[i].def.x << " " << Nodes[i].def.y << endl; }
     }

     double fpart = diffclock(fend, fbegin);
     ftime = fpart + ftime;

  }
  f2.close();

  // compute final area
  printf("Total focing computation time (s): %.10f\n", ftime);
}

// Draw starting configuration
void startCurve() {
  
  // file handling
  ifstream f1;
  ifstream f2;
  ifstream f3;
  f1.open("data/UnitCirclePointsN2048.txt");
  f2.open("data/UnitCircleTrianglesN2048.txt");

  // determine length
  int Npoints = -1;
  int Ntris = -1;

  string c1;
  string c2;

  while( !f1.eof() ) {
        getline(f1, c1);
        Npoints++;
  }
  f1.close();
  while( !f2.eof() ) {
        getline(f2, c2);
        Ntris++;
  }
  f2.close();

  f1.open("data/UnitCirclePointsN2048.txt");
  f2.open("data/UnitCircleTrianglesN2048.txt");
  f3.open("data/UnitCircleBoundaryN2048.txt");

  vector Nodes[Npoints];
  triangle Tris[Ntris];
  int Boundary[Npoints];

  int counter = 0;
  double d1, d2;
  while(f1 >> d1 >> d2) {
        //printf("(%f, %f)\n", d1, d2);
        Nodes[counter].x = d1;
        Nodes[counter].y = d2;
        counter++;
  }
  f1.close();

  counter = 0;
  int i1, i2, i3;
  while(f2 >> i1 >> i2 >> i3) {
        Tris[counter].A = i1-1;
        Tris[counter].B = i2-1;
        Tris[counter].C = i3-1;
        //printf("[%d %d %d]\n", Tris[counter].A, Tris[counter].B, Tris[counter].C);
        counter++;
  }
  f2.close();

  counter = 0;
  int ext;
  // set all points to interior points
  for(int k=0; k<Npoints; k++) {
        Boundary[k] = 0;
  }
  while(f3 >> ext) {
        Boundary[ext-1] = 1;
        counter++;
  }
  f3.close();

  // output to array of vertices and array of triangles
  vertex Points[Npoints];

  for(int i=0; i<Npoints; i++) {
        Points[i].ref.x = Nodes[i].x;
        Points[i].ref.y = Nodes[i].y;
        Points[i].exterior = Boundary[i];

        // SPECIFY DEFORMATION HERE // Step 0: NO deformation
        //Points[i].def.x = Nodes[i].x;
        //Points[i].def.y = Nodes[i].y;

        // SPECIFY DEFORMATION HERE // Step 1: LINEAR deformation
        ///// expansion /////
        Points[i].def.x = (1.0 - def)*Nodes[i].x;
        Points[i].def.y = (1.0 - def)*Nodes[i].y;
        ///// shear /////
        //Points[i].def.x = Nodes[i].x + lambda*Nodes[i].y;
        //Points[i].def.y = Nodes[i].y;
        ///// vertical stretch /////
        //Points[i].def.x = Nodes[i].x;
        //Points[i].def.y = (1.0 + lambda)*Nodes[i].y;
        ///// uniaxial extension /////
        //Points[i].def.x = lambda*Nodes[i].x;
        //Points[i].def.y = (1.0/lambda)*Nodes[i].y; 

        // SPECIFY DEFORMATION HERE // Step 2: NONLINEAR deformation
        //Points[i].def.x = lambda*Nodes[i].x*Nodes[i].x;
        //Points[i].def.y = Points[i].def.y;
  }


  for(int j=0; j<Ntris; j++) {
        // find vertices
        int iA = Tris[j].A; // index of A vertex
        int iB = Tris[j].B; // index of B vertex
        int iC = Tris[j].C; // index of C vertex
        Points[iA].ref.x = Nodes[iA].x;
        Points[iA].ref.y = Nodes[iA].y;
        Points[iB].ref.x = Nodes[iB].x;
        Points[iB].ref.y = Nodes[iB].y;
        Points[iC].ref.x = Nodes[iC].x;
        Points[iC].ref.y = Nodes[iC].y;
   }

   for(int k=0; k<Ntris; k++) {
        // find forcing terms that remain constant with any deformation and timestep
        ref_info(Npoints, Tris[k], Points);
   }

   progress(Npoints, Ntris, Points, Tris);
}

// Main
int main(int argc, char **argv) {
  clock_t begin = clock();
  startCurve();
  clock_t end = clock();
  printf("GPU computation time (ms): %.10f \n", gpuTime);
  printf("Total computation time (s): %.10f\n", double(diffclock(end,begin)));

  return 0;
}



