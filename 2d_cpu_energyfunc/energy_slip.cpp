/* CPU Benchmark Immersed Boundary Unstructured Grid
 * Based on "The Method of Regularized Stokelets" by Cortez
 * Elastic force computed using energy-based formulation by Devendran + Peskin
 * C.Copos 04/10/2013
 */

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h> 
#include <time.h> 
#include <cstdio>
#include <cstdlib>
#include <stdio.h>

using namespace std;

const double PI = 3.14159265358979323846264338f;
const int N = 64;
const float visc = 1.0f;
const float flow = 1.0f;
const float drag = 1.0f;
const float def = 0.1; 

// Lame constants
const double lambda = 1.0; // bulk modulus
const double mu = 1.0; // shear modulus (G)

// time 
const float TMAX = 10.0f;
const float tstep = 0.01f;
//static int stop = floor(TMAX/tstep);
static int stop = 2000; /* RUNNING STATIC TESTS */

// curve constants
const float r0 = 1.0f;
const float ri = 1.2f; 
const double ds = ri/(N-1); const double e = 1.2f*ds; // e: parameter determining width of blobs or cutoffs

// vector
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
  //printf("Jacobian: %f\n", detS);

  a.x1 = (1.0/detS)*( (nodes[tr.B].def.x-nodes[tr.A].def.x)*(nodes[tr.C].ref.y-nodes[tr.A].ref.y) + (nodes[tr.C].def.x-nodes[tr.A].def.x)*(nodes[tr.A].ref.y - nodes[tr.B].ref.y) );
  a.y1 = (1.0/detS)*( (nodes[tr.B].def.x-nodes[tr.A].def.x)*(nodes[tr.A].ref.x-nodes[tr.C].ref.x) + (nodes[tr.C].def.x-nodes[tr.A].def.x)*(nodes[tr.B].ref.x-nodes[tr.A].ref.x) );
  a.x2 = (1.0/detS)*( (nodes[tr.B].def.y-nodes[tr.A].def.y)*(nodes[tr.C].ref.y-nodes[tr.A].ref.y) + (nodes[tr.C].def.y-nodes[tr.A].def.y)*(nodes[tr.A].ref.y-nodes[tr.B].ref.y) ); 
  a.y2 = (1.0/detS)*( (nodes[tr.B].def.y-nodes[tr.A].def.y)*(nodes[tr.A].ref.x - nodes[tr.C].ref.x) + (nodes[tr.C].def.y-nodes[tr.A].def.y)*(nodes[tr.B].ref.x - nodes[tr.A].ref.x) ); 
  //debug
  //printf("nodes[tr.A].ref.x = %f, nodes[tr.A].ref.y = %f\n", nodes[tr.A].ref.x, nodes[tr.A].ref.y);
  //printf("nodes[tr.B].ref.x = %f, nodes[tr.B].ref.y = %f\n", nodes[tr.B].ref.x, nodes[tr.B].ref.y);
  //printf("nodes[tr.C].ref.x = %f, nodes[tr.C].ref.y = %f\n", nodes[tr.C].ref.x, nodes[tr.C].ref.y);
  //double x_cent = (0.3333)*(nodes[tr.A].def.x + nodes[tr.B].def.x + nodes[tr.C].def.x);
  //double y_cent = (0.3333)*(nodes[tr.A].def.y + nodes[tr.B].def.y + nodes[tr.C].def.y);
  //printf("x-centroid: %.6f, y-centroid: %.6f\n", x_cent, y_cent);

  //printf("nodes[tr.A].def.x = %f, nodes[tr.A].def.y = %f\n", nodes[tr.A].def.x, nodes[tr.A].def.y);
  //printf("nodes[tr.B].def.x = %f, nodes[tr.B].def.y = %f\n", nodes[tr.B].def.x, nodes[tr.B].def.y);
  //printf("nodes[tr.C].def.x = %f, nodes[tr.C].def.y = %f\n", nodes[tr.C].def.x, nodes[tr.C].def.y);

  //printf("------\n");
  //printf("triangle: (%d, %d, %d), a.x1 = %f\n", tr.A, tr.B, tr.C, a.x1);
  //printf("triangle: (%d, %d, %d), a.y1 = %f\n", tr.A, tr.B, tr.C, a.y1); 
  //printf("triangle: (%d, %d, %d), a.x2 = %f\n", tr.A, tr.B, tr.C, a.x2); 
  //printf("triangle: (%d, %d, %d), a.y2 = %f\n", tr.A, tr.B, tr.C, a.y2);
  //printf("detS: %f & area: %f\n", detS, tr.ref_area);
  /*
  // inverse transpose of deformation gradient tensor (w/o outside normalizers i.e. determinants) 
  matrix ait;
  ait.x1 = a.y2;
  ait.y1 = (-1.0)*(a.x2);
  ait.x2 = (-1.0)*(a.y1);
  ait.y2 = a.x1;

  // Material displacement gradient tensor ( = deformation gradient tensor - I)
  matrix d; 
  d.x1 = a.x1 - 1.0;
  d.y1 = a.y1;
  d.x2 = a.x2;
  d.y2 = a.y2 - 1.0;

  // Cauchy stress tensor
  matrix sigma;
  sigma.x1 = 2*lambda*(d.x1+d.y2) + 4.0*mu*d.x1; // multiplied by 2 from derivative of little u (displ vector)
  sigma.y1 = mu*(d.y1+d.x2);
  sigma.x2 = mu*(d.x2+d.y1);
  sigma.y2 = lambda*(d.x1+d.y2) + 2.0*mu*d.y2;  
  */

  // term 1 (otherwise known as 1st Piola-Kirchhoff tensor)
  //tr.f1.x1 = ( sigma.x1*ait.x1 + sigma.y1*ait.x2 );
  //tr.f1.y1 = ( sigma.x1*ait.y1 + sigma.y1*ait.y2 );
  //tr.f1.x2 = ( sigma.x2*ait.x1 + sigma.y2*ait.x2 );
  //tr.f1.y2 = ( sigma.x2*ait.y1 + sigma.y2*ait.y2 ); 
  // new (4/08 pm)
  tr.f1.x1 = (a.x1-1)*(lambda+2*mu) + lambda*(a.y2-1);
  tr.f1.y2 = (a.y2-1)*(lambda+2*mu) + lambda*(a.x1-1);
  tr.f1.x2 = mu*(a.x2+a.y1);
  tr.f1.y1 = mu*(a.x2+a.y1);

  //printf("a.y2 = %f, a.x1 = %f, tr.f1.y2 = %f\n", a.y2, a.x1, tr.f1.y2);

  //printf("%f %f %f %f %f %f %f \n", nodes[tr.A].ref.x, nodes[tr.B].ref.x, nodes[tr.C].ref.x, tr.f1.x1, tr.f1.y1, tr.f1.x2, tr.f1.y2);

  //printf("-----\n");
  //printf("triangle: (%d, %d, %d), f1.x1 = %f\n", tr.A, tr.B, tr.C, tr.f1.x1);
  //printf("triangle: (%d, %d, %d), f1.y1 = %f\n", tr.A, tr.B, tr.C, tr.f1.y1);
  //printf("triangle: (%d, %d, %d), f1.x2 = %f\n", tr.A, tr.B, tr.C, tr.f1.x2);
  //printf("triangle: (%d, %d, %d), f1.y2 = %f\n", tr.A, tr.B, tr.C, tr.f1.y2);
}

// Compute pressure at point i due to all other points in the unstructured grid
double pressure(int i, int Npts, vertex Nodes[], vector f[]) {
  double p = 0.0; 
  double r = sqrt(powf(Nodes[i].def.x,2) + powf(Nodes[i].def.y,2));
  
  for (int ii=0; ii<Npts; ii++) {
	double pk = sqrt(powf(Nodes[ii].def.x,2) + powf(Nodes[ii].def.y,2));
	double rk = sqrt(powf(Nodes[i].def.x-Nodes[ii].def.x,2) + powf(Nodes[i].def.y-Nodes[ii].def.y,2));
        double sq = sqrtf(powf(rk,2)+powf(e,2));
	double h = (1.0f)/(2.0f*PI) * (powf(rk,2)+2.0f*powf(e,2)+e*sq)/((sq+e)*powf(sq,3));

	p += (f[ii].x*(Nodes[i].def.x-Nodes[ii].def.x) + f[ii].y*(Nodes[i].def.y-Nodes[ii].def.y)) * h;
  }
}

// Compute velocity at point i due to all other points in the unstructured grid
vector velocity(int i, int Npts, vertex Nodes[], vector f[]) {
  vector v;
  v.x = 0.0f ; v.y = 0.0f;
  double pvx = 0.0f; double pvy = 0.0f; // partial component velocities 

  for (int ii=0; ii<Npts; ii++) {
    double rk = sqrt(powf(Nodes[i].def.x-Nodes[ii].def.x,2) + powf(Nodes[i].def.y-Nodes[ii].def.y,2));
    double sq = sqrtf(powf(rk,2)+powf(e,2));
    pvx = 0.0f; pvy = 0.0f;
    double p1 = (1.0f/(4.0f*visc*PI)) * (logf(sq+e)-(e*(sq+2*e))/(sq*(sq+e)));
    double p2 = (1.0f/(4.0f*visc*PI)) * (sq+2*e)/(sq*powf(sq+e,2));

    pvx = -1.0f*p1*f[ii].x + p2*(powf(Nodes[i].def.x-Nodes[ii].def.x,2)*f[ii].x + (Nodes[i].def.x-Nodes[ii].def.x)*(Nodes[i].def.y-Nodes[ii].def.y)*f[ii].y);
    pvy = -1.0f*p1*f[ii].y + p2*((Nodes[i].def.x-Nodes[ii].def.x)*(Nodes[i].def.y-Nodes[ii].def.y)*f[ii].x + powf(Nodes[i].def.y-Nodes[ii].def.y,2)*f[ii].y);

    v.x += pvx;
    v.y += pvy;
 }

 return v;
}

// Progression
void progress(int Npts, int Ntris, vertex Nodes[], triangle Tris[]) {
  vector pos_init[Npts]; 
  vector vel;
  vector v[Npts];
  vector fpertr[Ntris];
  vector f[Npts];
  double p[Npts];
  double vtime = 0.0f, ftime = 0.0f;
 
   //printf("Here\n");
  // file handling
  ofstream f1, f2, f3, f4, f5;
  f1.open("initial_pos_conf.txt"); f2.open("final_pos_conf.txt"); f3.open("ref_pos_conf.txt");
  f4.open("pressure.txt");
  f5.open("all.txt");

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

  // CYCLE THROUGH POINTS AND ADD CONTRACTILE FORCES
  for(int k=0; k<Npts; k++) {
    Nodes[k].ref.x = (1.0-0.1)*Nodes[k].ref.x;
    Nodes[k].ref.y = (1.0-0.1)*Nodes[k].ref.y;
  }

  for(int t=0; t<stop; t++) {
     clock_t fbegin = clock();
     float ref_Tarea = 0.0; float def_Tarea = 0.0;
  
     // CYCLE THROUGH TRIANGLES AND COMPUTE FORCES
     for(int j=0; j<Ntris; j++) {
	def_info(Npts, Tris[j], Nodes);
	// vertex A	
 	Nodes[Tris[j].A].force.x += -1.0*(Tris[j].ref_area)*(Tris[j].f1.x1*Tris[j].f2_0_x.x1 + Tris[j].f1.y1*Tris[j].f2_0_x.y1 + Tris[j].f1.x2*Tris[j].f2_0_x.x2 + Tris[j].f1.y2*Tris[j].f2_0_x.y2);
	Nodes[Tris[j].A].force.y += -1.0*(Tris[j].ref_area)*(Tris[j].f1.x1*Tris[j].f2_0_y.x1 + Tris[j].f1.y1*Tris[j].f2_0_y.y1 + Tris[j].f1.x2*Tris[j].f2_0_y.x2 + Tris[j].f1.y2*Tris[j].f2_0_y.y2);
	//printf("%d: %f, %f\n", Tris[j].A, -1.0*(Tris[j].ref_area)*(Tris[j].f1.x1*Tris[j].f2_0_x.x1 + Tris[j].f1.y1*Tris[j].f2_0_x.y1 + Tris[j].f1.x2*Tris[j].f2_0_x.x2 + Tris[j].f1.y2*Tris[j].f2_0_x.y2),  -1.0*(Tris[j].ref_area)*(Tris[j].f1.x1*Tris[j].f2_0_y.x1 + Tris[j].f1.y1*Tris[j].f2_0_y.y1 + Tris[j].f1.x2*Tris[j].f2_0_y.x2 + Tris[j].f1.y2*Tris[j].f2_0_y.y2));
	
	//printf("%.5f %.5f %.5f %.5f\n ", Tris[j].f1.x1, Tris[j].f1.y1, Tris[j].f1.x2, Tris[j].f1.y2);

	// vertex B
	Nodes[Tris[j].B].force.x += -1.0*(Tris[j].ref_area)*(Tris[j].f1.x1*Tris[j].f2_1_x.x1 + Tris[j].f1.y1*Tris[j].f2_1_x.y1 + Tris[j].f1.x2*Tris[j].f2_1_x.x2 + Tris[j].f1.y2*Tris[j].f2_1_x.y2);
	Nodes[Tris[j].B].force.y += -1.0*(Tris[j].ref_area)*(Tris[j].f1.x1*Tris[j].f2_1_y.x1 + Tris[j].f1.y1*Tris[j].f2_1_y.y1 + Tris[j].f1.x2*Tris[j].f2_1_y.x2 + Tris[j].f1.y2*Tris[j].f2_1_y.y2);
	//printf("%d: %f, %f\n", Tris[j].B, -1.0*(Tris[j].ref_area)*(Tris[j].f1.x1*Tris[j].f2_1_x.x1 + Tris[j].f1.y1*Tris[j].f2_1_x.y1 + Tris[j].f1.x2*Tris[j].f2_1_x.x2 + Tris[j].f1.y2*Tris[j].f2_1_x.y2), -1.0*(Tris[j].ref_area)*(Tris[j].f1.x1*Tris[j].f2_1_y.x1 + Tris[j].f1.y1*Tris[j].f2_1_y.y1 + Tris[j].f1.x2*Tris[j].f2_1_y.x2 + Tris[j].f1.y2*Tris[j].f2_1_y.y2));

	// vertex C
	Nodes[Tris[j].C].force.x += -1.0*(Tris[j].ref_area)*(Tris[j].f1.x1*Tris[j].f2_2_x.x1 + Tris[j].f1.y1*Tris[j].f2_2_x.y1 + Tris[j].f1.x2*Tris[j].f2_2_x.x2 + Tris[j].f1.y2*Tris[j].f2_2_x.y2);
        Nodes[Tris[j].C].force.y += -1.0*(Tris[j].ref_area)*(Tris[j].f1.x1*Tris[j].f2_2_y.x1 + Tris[j].f1.y1*Tris[j].f2_2_y.y1 + Tris[j].f1.x2*Tris[j].f2_2_y.x2 + Tris[j].f1.y2*Tris[j].f2_2_y.y2);	
	//printf("%d: %f, %f\n", Tris[j].C, -1.0*(Tris[j].ref_area)*(Tris[j].f1.x1*Tris[j].f2_2_x.x1 + Tris[j].f1.y1*Tris[j].f2_2_x.y1 + Tris[j].f1.x2*Tris[j].f2_2_x.x2 + Tris[j].f1.y2*Tris[j].f2_2_x.y2), -1.0*(Tris[j].ref_area)*(Tris[j].f1.x1*Tris[j].f2_2_y.x1 + Tris[j].f1.y1*Tris[j].f2_2_y.y1 + Tris[j].f1.x2*Tris[j].f2_2_y.x2 + Tris[j].f1.y2*Tris[j].f2_2_y.y2));

	// debug
	//printf("TERM2\n");
	//printf("%d %d %d\n", Tris[j].A, Tris[j].B, Tris[j].C);
	//printf("t2_0.x: (%f %f %f %f)\n", Tris[j].f2_0_x.x1, Tris[j].f2_0_x.y1, Tris[j].f2_0_x.x2, Tris[j].f2_0_x.y2);
	//printf("t2_0.y: (%f %f %f %f)\n", Tris[j].f2_0_y.x1, Tris[j].f2_0_y.y1, Tris[j].f2_0_y.x2, Tris[j].f2_0_y.y2);
	//printf("t2_1.x: (%f %f %f %f)\n", Tris[j].f2_1_x.x1, Tris[j].f2_1_x.y1, Tris[j].f2_1_x.x2, Tris[j].f2_1_x.y2);
	//printf("t2_1.y: (%f %f %f %f)\n", Tris[j].f2_1_y.x1, Tris[j].f2_1_y.y1, Tris[j].f2_1_y.x2, Tris[j].f2_1_y.y2);
	//printf("t2_2.x: (%f %f %f %f)\n", Tris[j].f2_2_x.x1, Tris[j].f2_2_x.y1, Tris[j].f2_2_x.x2, Tris[j].f2_2_x.y2);
	//printf("t2_2.y: (%f %f %f %f)\n", Tris[j].f2_2_y.x1, Tris[j].f2_2_y.y1, Tris[j].f2_2_y.x2, Tris[j].f2_2_y.y2);
  	//printf("FORCES\n");
	//printf("f_0.x: %f\n", Nodes[Tris[j].A].force.x);
	//printf("f_0.y: %f\n", Nodes[Tris[j].A].force.y);
	//printf("f_1.x: %f\n", Nodes[Tris[j].B].force.x);
	//printf("f_1.y: %f\n", Nodes[Tris[j].B].force.y);
	//printf("f_2.x: %f\n", Nodes[Tris[j].C].force.x);
	//printf("f_2.y: %f\n", Nodes[Tris[j].C].force.y); 
	 
	ref_Tarea += Tris[j].ref_area;
	def_Tarea += Tris[j].def_area;
     } 
     // DONE CYCLING THROUGH TRIANGLES
     //printf("Reference total area: %f\n", ref_Tarea);
     //printf("Deformed total area: %f\n", def_Tarea); 
     clock_t fend = clock();

    for(int k=0; k<Npts; k++) {
		f[k].x = Nodes[k].force.x;
		f[k].y = Nodes[k].force.y;
   		Nodes[k].force.x = 0.0;	
		Nodes[k].force.y = 0.0; 
    }

    // for each node in unstructured mesh
    for(int i=0; i<Npts; i++) {		
		// compute pressure	
		p[i] = pressure(i, Npts, Nodes, f);
		f4 << Nodes[i].def.x << " " <<  Nodes[i].def.y << " " << p[i] << endl;	

		// compute velocity
		clock_t vbegin = clock();
		vel = velocity(i, Npts, Nodes, f);
		clock_t vend = clock();

		//printf("@ node %d (exterior: %d) the force is (%f, %f)\n", i, Nodes[i].exterior, f[i].x, f[i].y);	
		//if (Nodes[i].exterior == 0) { // interior point
		//   printf("%f %f %f %f\n", Nodes[i].ref.x, Nodes[i].ref.y, f[i].x, f[i].y);
		//}
		//printf(" expected force (%f, %f)\n", 2.0*lambda*def+4.0*mu*def, 0.0);
		v[i].x = vel.x + f[i].x/drag;
		v[i].y = vel.y + f[i].y/drag;

		//printf("@ node %d (exterior: %d) the velocity is (%f, %f)\n", i, Nodes[i].exterior, v[i].x, v[i].y);
		// vels
		//printf("(vx, vy) = (%.10f, %.10f)\n", vel.x, vel.y);
		//printf("%.16f %.16f %.16f %.16f %.16f %.16f\n", Nodes[i].def.x, Nodes[i].def.y, v[i].x, v[i].y, f[i].x, f[i].y);
		//printf("(fx, fy) = (%.10f, %.10f)\n", f[i].x, f[i].y);

		// timings
		double partv = diffclock(vend, vbegin);
		vtime = vtime + partv;
	}

	for(int i=0; i<Npts; i++) {
		Nodes[i].def.x = Nodes[i].def.x + tstep*v[i].x;
        Nodes[i].def.y = Nodes[i].def.y + tstep*v[i].y;
		if(t==stop-1) { f2 << Nodes[i].def.x << " " << Nodes[i].def.y << endl; }
    
		// print all to file
        f5 << t << " " << Nodes[i].def.x << " " << Nodes[i].def.y << " " << f[i].x << " " << f[i].y << " " << v[i].x << " " << v[i].y << endl;
     }

     double fpart = diffclock(fend, fbegin);
     ftime = fpart + ftime;

  }
  f2.close();
  f4.close();
  f5.close();

  // print out times
  //printf("Total forcing computation time (s): %.10f\n", ftime);
  //printf("Total velocity computation time (s): %.10f\n", vtime);
}

// Draw starting configuration
void startCurve() {
  
  // file handling
  ifstream f1;
  ifstream f2;
  ifstream f3;
  f1.open("UniformMesh/UnitCirclePointsN64.txt");
  f2.open("UniformMesh/UnitCircleTrianglesN64.txt");

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

  printf("Number of points: %d and number of triangles: %d\n", Npoints, Ntris);

  f1.open("UniformMesh/UnitCirclePointsN64.txt");
  f2.open("UniformMesh/UnitCircleTrianglesN64.txt");
  f3.open("UniformMesh/UnitCircleBoundaryN64.txt");

  vector Nodes[Npoints];
  triangle Tris[Ntris];
  int Boundary[Npoints];

  int counter = 0;
  double d1, d2;
  while(f1 >> d1 >> d2) {
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
    	//printf("%d %d %d\n",i1-1,i2-1,i3-1);
	    counter++;
  }
  f2.close();

  int nBoundary = 0;
  int ext;
  // set all points to interior points
  for(int k=0; k<Npoints; k++) {
	Boundary[k] = 0;
  }
  while(f3 >> ext) {
	Boundary[ext-1] = 1;
  	nBoundary++;
  } 
  f3.close();

  // output to array of vertices and array of triangles
  vertex Points[Npoints];

  for(int i=0; i<Npoints; i++) {
	Points[i].ref.x = Nodes[i].x;
	Points[i].ref.y = Nodes[i].y;
	Points[i].exterior = Boundary[i];
	
	// SPECIFY DEFORMATION HERE // Step 0: NO deformation
	Points[i].def.x = Nodes[i].x;
	Points[i].def.y = Nodes[i].y;
	
	// SPECIFY DEFORMATION HERE // Step 1: LINEAR deformation
	///// contraction /////
	//Points[i].def.x = (1.0 - def)*Nodes[i].x;
	//Points[i].def.y = (1.0 - def)*Nodes[i].y;
	///// expansion /////
  	//Points[i].def.x = (1.0 + def)*Nodes[i].x;
	//Points[i].def.y = (1.0 + def)*Nodes[i].y;
	///// shear /////
	//Points[i].def.x = Nodes[i].x + def*Nodes[i].y;
 	//Points[i].def.y = Nodes[i].y;
	///// vertical stretch /////
	//Points[i].def.x = Nodes[i].x;
	//Points[i].def.y = (1.0 + def)*Nodes[i].y;
	///// uniaxial extension /////
 	//Points[i].def.x = def*Nodes[i].x;
	//Points[i].def.y = (1.0/def)*Nodes[i].y; 
	
	// SPECIFY DEFORMATION HERE // Step 2: NONLINEAR deformation
	///// stretch /////
 	//Points[i].def.x = Nodes[i].x + def*Nodes[i].x*Nodes[i].x;
	//Points[i].def.y = Nodes[i].y;
	///// shear /////
	//Points[i].def.x = Nodes[i].x + def*Nodes[i].y*Nodes[i].y;
	//Points[i].def.y = Nodes[i].y;  	
  } 

  for(int j=0; j<Ntris; j++) {
 	// find vertices
	int iA = Tris[j].A; // index of A vertex
	int iB = Tris[j].B; // index of B vertex
	int iC = Tris[j].C; // index of C vertex
	 //printf("%d vertices %d %d %d\n",j,iA,iB,iC);  
	Points[iA].ref.x = Nodes[iA].x;
	Points[iA].ref.y = Nodes[iA].y;
	Points[iB].ref.x = Nodes[iB].x;
        Points[iB].ref.y = Nodes[iB].y;
	Points[iC].ref.x = Nodes[iC].x;
        Points[iC].ref.y = Nodes[iC].y;
   }

   // find focing terms that remain constant with any deformation and timestep
   for(int k=0; k<Ntris; k++) {	
	ref_info(Npoints, Tris[k], Points);
   }

   // progress points forward in time
   progress(Npoints, Ntris, Points, Tris);
}

// Main
int main(int argc, char **argv) {
  //printf("def = %f\n", def);
  clock_t begin = clock();
  startCurve();
  clock_t end = clock();
  //printf("Total computation time (s): %.10f\n", double(diffclock(end,begin)));

  return 0;
}

