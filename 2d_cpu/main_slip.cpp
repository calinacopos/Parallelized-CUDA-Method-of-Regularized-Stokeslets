/* main_slip.cpp 
 * GPU Benchmark Immersed Boundary Structured Grid
 * Based on "The Method of Regularized Stokelets" by R.Cortez
 * C.Copos 6/25/2012: switched to double precision
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

const double PI = 3.14159265358979323846264338;
const int N = 4; // # side square partitions 
const double visc = 1.0;
const double flow = 1.0;
const double drag = 0.0001;

// time 
const double TMAX = 0.1;
const double tstep = 0.0001;
static int stop = floor(TMAX/tstep);

// curve constants
const double r0 = 1.0;
const double ri = 1.2; 
const double ds = ri/(N-1); const double e = 1.2*ds; // e: parameter determining width of blobs or cutoffs

// vector structure
struct vector{
  double x; // x-component
  double y; // y-component
};

// spring constants
const double kl = 1.0*ds/*0.1*ds*/;
const double ks = 0.5*kl; // from victor camacho's elasticity paper

// Compute time difference in seconds
double diffclock(clock_t s, clock_t e) {
   double diffticks = s-e;
   double diffms = (diffticks)/CLOCKS_PER_SEC;
   return diffms;
}

// Compute average elastic force at point ij due to its nearest neighboring points
vector eforce(int i, int j, double x[N*N], double y[N*N]) {
  vector f; // final force
  double rl = r0/(N-1); // rest length for longitudinal springs
  double rd = sqrtf(2.0)*rl; // rest length for diagonal springs

  int i_0, i_1, i_2, j_0, j_1, j_2;
  i_0 = i-1; 
  i_1 = i; 
  i_2 = i+1; 
  j_0 = j-1; 
  j_1 = j;
  j_2 = j+1;

  double dlk_00, dlk_01, dlk_02, dlk_10, dlk_12, dlk_20, dlk_21, dlk_22;
  vector f_00, f_01, f_02, f_10, f_12, f_20, f_21, f_22;

  // distance between point (i,j) = (i_1,j_1) and points ...
    // top left corner
    if (i_1==0 && j_1==0) { 
	dlk_00 = 0.0; dlk_10 = 0.0; dlk_20 = 0.0; dlk_01 = 0.0; dlk_02 = 0.0;
	dlk_21 = sqrtf(powf(x[i_1+N*j_1]-x[i_2+N*j_1],2) + powf(y[i_1+N*j_1]-y[i_2+N*j_1],2));
	dlk_12 = sqrtf(powf(x[i_1+N*j_1]-x[i_1+N*j_2],2) + powf(y[i_1+N*j_1]-y[i_1+N*j_2],2));
	dlk_22 = sqrtf(powf(x[i_1+N*j_1]-x[i_2+N*j_2],2) + powf(y[i_1+N*j_1]-y[i_2+N*j_2],2));
    }
    // top right corner
    else if (i_1==(N-1) && j_1==0) {
    	dlk_00 = 0.0; dlk_10 = 0.0; dlk_20 = 0.0; dlk_21 = 0.0; dlk_22 = 0.0;
 	dlk_01 = sqrtf(powf(x[i_1+N*j_1]-x[i_0+N*j_1],2) + powf(y[i_1+N*j_1]-y[i_0+N*j_1],2));
	dlk_02 = sqrtf(powf(x[i_1+N*j_1]-x[i_0+N*j_2],2) + powf(y[i_1+N*j_1]-y[i_0+N*j_2],2));
	dlk_12 = sqrtf(powf(x[i_1+N*j_1]-x[i_1+N*j_2],2) + powf(y[i_1+N*j_1]-y[i_1+N*j_2],2));
    }
    // bottom left corner
    else if (i_1==0 && j_1==(N-1)) {
	dlk_00 = 0.0; dlk_01 = 0.0; dlk_02 = 0.0; dlk_12 = 0.0; dlk_22 = 0.0;
	dlk_10 = sqrtf(powf(x[i_1+N*j_1]-x[i_1+N*j_0],2) + powf(y[i_1+N*j_1]-y[i_1+N*j_0],2));
	dlk_20 = sqrtf(powf(x[i_1+N*j_1]-x[i_2+N*j_0],2) + powf(y[i_1+N*j_1]-y[i_2+N*j_0],2));
	dlk_21 = sqrtf(powf(x[i_1+N*j_1]-x[i_2+N*j_1],2) + powf(y[i_1+N*j_1]-y[i_2+N*j_1],2));
    }
    // bottom right corner
    else if (i_1==(N-1) && j_1==(N-1)) {
    	dlk_20 = 0.0; dlk_21 = 0.0; dlk_22 = 0.0; dlk_12 = 0.0; dlk_02 = 0.0;
	dlk_00 = sqrtf(powf(x[i_1+N*j_1]-x[i_0+N*j_0],2) + powf(y[i_1+N*j_1]-y[i_0+N*j_0],2));
	dlk_10 = sqrtf(powf(x[i_1+N*j_1]-x[i_1+N*j_0],2) + powf(y[i_1+N*j_1]-y[i_1+N*j_0],2));
	dlk_01 = sqrtf(powf(x[i_1+N*j_1]-x[i_0+N*j_1],2) + powf(y[i_1+N*j_1]-y[i_0+N*j_1],2));
    }
    // top edge
    else if (j_1==0) {
	dlk_00 = 0.0; dlk_10 = 0.0; dlk_20 = 0.0;
	dlk_01 = sqrtf(powf(x[i_1+N*j_1]-x[i_0+N*j_1],2) + powf(y[i_1+N*j_1]-y[i_0+N*j_1],2));
	dlk_21 = sqrtf(powf(x[i_1+N*j_1]-x[i_2+N*j_1],2) + powf(y[i_1+N*j_1]-y[i_2+N*j_1],2));
	dlk_02 = sqrtf(powf(x[i_1+N*j_1]-x[i_0+N*j_2],2) + powf(y[i_1+N*j_1]-y[i_0+N*j_2],2));
        dlk_12 = sqrtf(powf(x[i_1+N*j_1]-x[i_1+N*j_2],2) + powf(y[i_1+N*j_1]-y[i_1+N*j_2],2));
        dlk_22 = sqrtf(powf(x[i_1+N*j_1]-x[i_2+N*j_2],2) + powf(y[i_1+N*j_1]-y[i_2+N*j_2],2));
    }
    // right edge
    else if (i_1==(N-1)) {
	dlk_20 = 0.0; dlk_21 = 0.0; dlk_22 = 0.0;
	dlk_00 = sqrtf(powf(x[i_1+N*j_1]-x[i_0+N*j_0],2) + powf(y[i_1+N*j_1]-y[i_0+N*j_0],2));
        dlk_10 = sqrtf(powf(x[i_1+N*j_1]-x[i_1+N*j_0],2) + powf(y[i_1+N*j_1]-y[i_1+N*j_0],2));
	dlk_01 = sqrtf(powf(x[i_1+N*j_1]-x[i_0+N*j_1],2) + powf(y[i_1+N*j_1]-y[i_0+N*j_1],2));
	dlk_02 = sqrtf(powf(x[i_1+N*j_1]-x[i_0+N*j_2],2) + powf(y[i_1+N*j_1]-y[i_0+N*j_2],2));
        dlk_12 = sqrtf(powf(x[i_1+N*j_1]-x[i_1+N*j_2],2) + powf(y[i_1+N*j_1]-y[i_1+N*j_2],2));
    }
    // bottom edge
    else if (j_1==(N-1)) {
	dlk_02 = 0.0; dlk_12 = 0.0; dlk_22 = 0.0;
	dlk_00 = sqrtf(powf(x[i_1+N*j_1]-x[i_0+N*j_0],2) + powf(y[i_1+N*j_1]-y[i_0+N*j_0],2));
        dlk_10 = sqrtf(powf(x[i_1+N*j_1]-x[i_1+N*j_0],2) + powf(y[i_1+N*j_1]-y[i_1+N*j_0],2));
        dlk_20 = sqrtf(powf(x[i_1+N*j_1]-x[i_2+N*j_0],2) + powf(y[i_1+N*j_1]-y[i_2+N*j_0],2));
        dlk_01 = sqrtf(powf(x[i_1+N*j_1]-x[i_0+N*j_1],2) + powf(y[i_1+N*j_1]-y[i_0+N*j_1],2));
        dlk_21 = sqrtf(powf(x[i_1+N*j_1]-x[i_2+N*j_1],2) + powf(y[i_1+N*j_1]-y[i_2+N*j_1],2));
    }
    // left edge
    else if (i_1==0) {
	dlk_00 = 0.0; dlk_01 = 0.0; dlk_02 = 0.0;
	dlk_10 = sqrtf(powf(x[i_1+N*j_1]-x[i_1+N*j_0],2) + powf(y[i_1+N*j_1]-y[i_1+N*j_0],2));
        dlk_20 = sqrtf(powf(x[i_1+N*j_1]-x[i_2+N*j_0],2) + powf(y[i_1+N*j_1]-y[i_2+N*j_0],2));
	dlk_21 = sqrtf(powf(x[i_1+N*j_1]-x[i_2+N*j_1],2) + powf(y[i_1+N*j_1]-y[i_2+N*j_1],2));
	dlk_12 = sqrtf(powf(x[i_1+N*j_1]-x[i_1+N*j_2],2) + powf(y[i_1+N*j_1]-y[i_1+N*j_2],2));
        dlk_22 = sqrtf(powf(x[i_1+N*j_1]-x[i_2+N*j_2],2) + powf(y[i_1+N*j_1]-y[i_2+N*j_2],2));
    }
    // interior
    else {
 	dlk_00 = sqrtf(powf(x[i_1+N*j_1]-x[i_0+N*j_0],2) + powf(y[i_1+N*j_1]-y[i_0+N*j_0],2));
        dlk_10 = sqrtf(powf(x[i_1+N*j_1]-x[i_1+N*j_0],2) + powf(y[i_1+N*j_1]-y[i_1+N*j_0],2));
        dlk_20 = sqrtf(powf(x[i_1+N*j_1]-x[i_2+N*j_0],2) + powf(y[i_1+N*j_1]-y[i_2+N*j_0],2));
    	dlk_01 = sqrtf(powf(x[i_1+N*j_1]-x[i_0+N*j_1],2) + powf(y[i_1+N*j_1]-y[i_0+N*j_1],2));
        dlk_21 = sqrtf(powf(x[i_1+N*j_1]-x[i_2+N*j_1],2) + powf(y[i_1+N*j_1]-y[i_2+N*j_1],2));
        dlk_02 = sqrtf(powf(x[i_1+N*j_1]-x[i_0+N*j_2],2) + powf(y[i_1+N*j_1]-y[i_0+N*j_2],2));
   	dlk_12 = sqrtf(powf(x[i_1+N*j_1]-x[i_1+N*j_2],2) + powf(y[i_1+N*j_1]-y[i_1+N*j_2],2));
        dlk_22 = sqrtf(powf(x[i_1+N*j_1]-x[i_2+N*j_2],2) + powf(y[i_1+N*j_1]-y[i_2+N*j_2],2)); 
   }
  
  // finally, compute forces
  if (dlk_00 == 0.0) { f_00.x = 0.0; f_00.y = 0.0; }
  else {
	f_00.x = -ks*(dlk_00-rd)*((x[i_1+N*j_1] - x[i_0+N*j_0])/dlk_00); 
  	f_00.y = -ks*(dlk_00-rd)*((y[i_1+N*j_1] - y[i_0+N*j_0])/dlk_00);
  }

  if (dlk_10 == 0.0) { f_10.x = 0.0; f_10.y = 0.0; }
  else {
	f_10.x = -kl*(dlk_10-rl)*((x[i_1+N*j_1] - x[i_1+N*j_0])/dlk_10);
	f_10.y = kl*(dlk_10-rl)*((-y[i_1+N*j_1] + y[i_1+N*j_0])/dlk_10);
  }

  if (dlk_20 == 0.0) { f_20.x = 0.0; f_20.y = 0.0; }
  else {
	f_20.x = -ks*(dlk_20-rd)*((x[i_1+N*j_1] - x[i_2+N*j_0])/dlk_20);
        f_20.y = -ks*(dlk_20-rd)*((y[i_1+N*j_1] - y[i_2+N*j_0])/dlk_20);
  	
  }
  if (dlk_01 == 0.0) { f_01.x = 0.0; f_01.y = 0.0; }
  else {
	f_01.x = -kl*(dlk_01-rl)*((x[i_1+N*j_1] - x[i_0+N*j_1])/dlk_01);
        f_01.y = kl*(dlk_01-rl)*((-y[i_1+N*j_1] + y[i_0+N*j_1])/dlk_01);
  }
  if (dlk_21 == 0.0) { f_21.x = 0.0; f_21.y = 0.0; }
  else {
	f_21.x = -kl*(dlk_21-rl)*((x[i_1+N*j_1] - x[i_2+N*j_1])/dlk_21); 
	f_21.y = kl*(dlk_21-rl)*((-y[i_1+N*j_1] + y[i_2+N*j_1])/dlk_21);
  }
  if (dlk_02 == 0.0) { f_02.x = 0.0; f_02.y = 0.0; }
  else {
	f_02.x = -1.0*ks*(dlk_02-rd)*((x[i_1+N*j_1] - x[i_0+N*j_2])/dlk_02);
	f_02.y = -1.0*ks*(dlk_02-rd)*((y[i_1+N*j_1] - y[i_0+N*j_2])/dlk_02);	
  }
  if (dlk_12 == 0.0) { f_12.x = 0.0; f_12.y = 0.0; }
  else {
	f_12.x = -1.0*kl*(dlk_12-rl)*((x[i_1+N*j_1] - x[i_1+N*j_2])/dlk_12);
        f_12.y = 1.0*kl*(dlk_12-rl)*((-y[i_1+N*j_1] + y[i_1+N*j_2])/dlk_12);
  }
  if (dlk_22 == 0.0) { f_22.x = 0.0; f_22.y = 0.0; }
  else {
	f_22.x = -1.0*ks*(dlk_22-rd)*((x[i_1+N*j_1] - x[i_2+N*j_2])/dlk_22);
        f_22.y = -1.0*ks*(dlk_22-rd)*((y[i_1+N*j_1] - y[i_2+N*j_2])/dlk_22);	
  }

  // evaluate final force components
  f.x = f_00.x + f_10.x + f_20.x + f_01.x + f_21.x + f_02.x + f_12.x + f_22.x;
  f.y = f_00.y + f_10.y + f_20.y + f_01.y + f_21.y + f_02.y + f_12.y + f_22.y;

  // what's going on with the forces?
  //printf("%f %f %f %f\n", x[i][j], y[i][j], f.x, f.y);
  ///*
  //printf("Force @ position: (%f,%f) is: (%f,%f)\n", x[i][j], y[i][j], f.x, f.y);
  //printf("Force due to (0,0) neighbor is: (%f,%f)\n", f_00.x, f_00.y);
  //printf("Force due to (1,0) neighbor is: (%f,%f)\n", f_10.x, f_10.y);
  //printf("Force due to (2,0) neighbor is: (%f,%f)\n", f_20.x, f_20.y);
  //printf("Force due to (0,1) neighbor is: (%f,%f)\n", f_01.x, f_01.y);
  //printf("Force due to (2,1) neighbor is: (%f,%f)\n", f_21.x, f_21.y);
  //printf("Force due to (0,2) neighbor is: (%f,%f)\n", f_02.x, f_02.y);
  //printf("Force due to (1,2) neighbor is: (%f,%f)\n", f_12.x, f_12.y);
  //printf("Force due to (2,2) neighbor is: (%f,%f)\n", f_22.x, f_22.y);
  //*/
 
  return f;   
}

// Compute pressure at point ij due to all other points in the structured grid
/*double pressure(int i, int j, double x[N*N], double y[N*N]) {
  double p; vector f;
  double pp = 0.0f; // partial pressure
  double rk = 0.0f; double sq = 0.0f;

  for (int jj=0; jj<N; jj++) { // loop over all nodes in the grid
    for (int ii=0; ii<N; ii++) { 
  	f.x = eforce(ii, jj, x, y).x;
        f.y = eforce(ii, jj, x, y).y;

	double r = sqrtf(powf(x[i][j],2) + powf(y[i][j],2));
	double pk = sqrtf(powf(x[ii][jj],2)+powf(y[ii][jj],2));
        double theta = atan2f(y[i][j], x[i][j]);
        double thetaj = atan2f(y[ii][jj], x[ii][jj]);
        
	double dtheta;
	if (theta>PI) { dtheta = thetaj + 2*PI - theta; }
        else { dtheta = theta - thetaj; }

	// dealing with rounding off errors     
        if (dtheta > -0.000001 && dtheta < 0.000001) { rk = 0.0f; }
        else { rk = sqrtf(powf(r,2) + powf(pk,2) - 2*r*pk*cosf(dtheta)); } // rk^2 = |x-xk|^2 = r^2 + pk^2 - 2*r*pk*cos(theta-thetak) where r=|x|, pk=|xk| 
        sq = sqrtf(powf(rk,2)+powf(e,2));

	double h = (1.0f)/(2.0f*PI) * (powf(rk,2)+2.0f*powf(e,2)+e*sq)/((sq+e)*powf(sq,1.5f));
	
	pp = (f.x*(x[i][j]-x[ii][jj]) + f.y*(y[i][j]-y[ii][jj])) * h;
	p += pp;
    }
  }
 
  return p; 
}*/

// Compute velocity at point ij due to all other points in the structured grid
vector velocity(int i, int j, double x[N*N], double y[N*N]) {
  vector v; vector f;
  v.x = 0.0; v.y = 0.0; f.x = 0.0; f.y = 0.0;
  double pvx = 0.0; double pvy = 0.0; // partial component velocities 

  for (int jj=0; jj<N; jj++) {
    for (int ii=0; ii<N; ii++) { // loop over all nodes in the grid
	f.x = eforce(ii, jj, x, y).x;
	f.y = eforce(ii, jj, x, y).y;

	/* // OLD WAY OF CALCULATION Rk 
        double r = sqrt(powf(x[i][j],2) + powf(y[i][j],2));
        double pk = sqrt(powf(x[ii][jj],2) + powf(y[ii][jj],2));
	double theta = atan2f(y[i][j], x[i][j]);
	double thetaj = atan2(y[ii][jj], x[ii][jj]);
        double dtheta;
        if (theta>PI) { dtheta = thetaj + 2*PI - theta; }
        else { dtheta = theta - thetaj; }

        // dealing with rounding off errors     
        if (dtheta > -0.000001 && dtheta < 0.000001) { rk = 0.0f; }
        else { rk = sqrtf(powf(r,2) + powf(pk,2) - 2*r*pk*cosf(dtheta)); } // rk^2 = |x-xk|^2 = r^2 + pk^2 - 2*r*pk*cos(theta-thetak) where r=|x|, pk=|xk|
        */

        double rk = sqrt(powf(x[i+N*j]-x[ii+N*jj],2) + powf(y[i+N*j]-y[ii+N*jj],2));	 
    	double sq = sqrtf(powf(rk,2)+powf(e,2));
    	pvx = 0.0; pvy = 0.0;
    	double p1 = (1.0/(4.0*visc*PI)) * (logf(sq+e)-(e*(sq+2*e))/(sq*(sq+e)));
    	double p2 = (1.0/(4.0*visc*PI)) * (sq+2*e)/(sq*powf(sq+e,2));

    	pvx = -1.0*p1*f.x + p2*(powf(x[i+N*j]-x[ii+N*jj],2)*f.x + (x[i+N*j]-x[ii+N*jj])*(y[i+N*j]-y[ii+N*jj])*f.y);
    	pvy = -1.0*p1*f.y + p2*((x[i+N*j]-x[ii+N*jj])*(y[i+N*j]-y[ii+N*jj])*f.x + powf(y[i+N*j]-y[ii+N*jj],2)*f.y);
	//printf("stokeslet: (%.10f, %.10f)\n", pvx, pvy);

    	v.x += pvx;
    	v.y += pvy;
    }
 }

 return v;
}

// Compute polygonal area
//double area(double x[N*N], double y[N*N]) {
//   double A = 0.0f;
//   int k1; int k2; // k1 = k; k2 = k+1
//   double x_ext[4*N-4]; double y_ext[4*N-4]; // exterior points

   /* Check for shear that we obtain a perfectly rotate uniformly discretized grid
   // check horizontal sides
   for(int j=0; j<N; j++) {
     for(int i=0; i<N-1; i++) {
	double x_current = x[i][j]; double x_next = x[i+1][j];
	double y_current = y[i][j]; double y_next = y[i+1][j];
	double horiz_side = sqrtf(powf(x_current-x_next,2) + powf(y_current-y_next,2));
	printf("h(%d, %d) = %f\n", i, j, horiz_side);
     }
   }
   
   // check vertical sides
   for(int i=0; i<N; i++) {
     for(int j=0; j<N-1; j++) {
	double x_current = x[i][j]; double x_next = x[i][j+1];
     	double y_current = y[i][j]; double y_next = y[i][j+1];
	double vert_side = fabs(y_next - y_current);
	printf("v(%d, %d) = %f\n", i, j, vert_side);
     }
   }
   */

/*   // rearrange points in vector that contain exterior points
   for(int j=0; j<N; j++) {
     for(int i=0; i<N; i++) {
	if(j==0) { x_ext[i] = x[i][j]; y_ext[i] = y[i][j]; } // bottom
	else if((i==(N-1)) && (j!=0)) { x_ext[j+i] = x[i][j]; y_ext[j+i] = y[i][j]; } // right
	else if((j==(N-1)) && (i!=(N-1))) { x_ext[3*j-i] = x[i][j]; y_ext[3*j-i] = y[i][j]; } // bottom
	else if((i==0) && (j!=0) && (j!=(N-1))) { x_ext[4*(N-1)-j] = x[i][j]; y_ext[4*(N-1)-j] = y[i][j]; } // left
     } 
   }
  
   for(int k=0; k<(4*N-4); k++) {
        k1 = k;
        if(k1 == (4*N-5)) { k2 = 0; }
        else k2 = k+1;
        A += 0.5f * (x_ext[k1]*y_ext[k2]-x_ext[k2]*y_ext[k1]);
}

   return A;
}*/

// Progression
void progress(double x[N*N], double y[N*N]) {
  double x_i[N*N]; double y_i[N*N];
  double u[N*N]; double v[N*N];
  vector vel; vector ef; 
  double p;
  double vtime = 0.0, ftime = 0.0;

  // file handling
  ofstream f1, f2;
  f1.open("initial_conf.txt"); f2.open("final_conf.txt"); 

  // initial configuration
  for(int j=0; j<N; j++) {
        for(int i=0; i<N; i++) {
           x_i[i+N*j] = x[i+N*j]; y_i[i+N*j] = y[i+N*j]; // saving initial configuration for initial area calculation
           f1 << x_i[i+N*j] << " " << y_i[i+N*j] << endl;
	}
  }
 
  for(int t=0; t<stop; t++) {
     for(int j=0; j<N; j++) {
        for(int i=0; i<N; i++) {
	  //p = pressure(i, j, x, y);
	  clock_t vbegin = clock();
	  vel = velocity(i, j, x, y);
	  clock_t vend = clock();
	  clock_t fbegin = clock();
	  ef = eforce(i, j, x, y);
	  clock_t fend = clock();
	  
	  //printf("(vx,vy) = (%.10f, %.10f)\n", vel.x, vel.y);
	
	  u[i+N*j] = vel.x + ef.x/drag;
	  v[i+N*j] = vel.y + ef.y/drag;

	  // timings
	  double partv = diffclock(vend, vbegin);
	  double partf = diffclock(fend, fbegin);
	  vtime = vtime + partv; ftime = ftime + partf;

	  x[i+N*j] = x[i+N*j] + tstep*u[i+N*j];
          y[i+N*j] = y[i+N*j] + tstep*v[i+N*j];
	  if(t==stop-1) { f2 << x[i+N*j] << " " << y[i+N*j] << endl; }
	}
     }

  }

  // print out times
  printf("Total focing computation time (s): %.10f\n", ftime);
  printf("Total velocity computation time (s): %.10f\n", vtime);

  f1.close(); f2.close();
  
  // compute final area
  //printf("Starting area: %.16f, Final area: %.16f\n", fabs(area(x_i,y_i)), fabs(area(x,y)) );
}

// Draw starting configuration
void startCurve() {
  double xf[N*N]; double yf[N*N];

  ///*
  // stretch uniformly in each direction
  for(int j=0; j<N; j++) {
    for(int i=0; i<N; i++) {
 	double dx = ri/(N-1); double dy = ri/(N-1);
	xf[i+j*N] = -ri/2 + dx*i;
	yf[i+j*N] = ri/2 - dy*j;
    }
  }
  //*/
  
  /* 
  // stretch in y-direction only
  for(int j=0; j<N; j++) {
     for(int i=0; i<N; i++) {
	double dx = r0/(N-1); double dy = ri/(N-1);
	xf[i][j] = -r0/2 + dx*i;
	yf[i][j] = ri/2 - dy*j;
	glVertex2f(xf[i][j], yf[i][j]);
     }
  }
  */

  /* 
  // shear in both directions
  double lambda = 0.5f; // shear element
  for(int j=0; j<N; j++) { 
    for(int i=0; i<N; i++) {
        double dx = r0/(N-1); double dy = r0/(N-1);
        yf[i][j] = r0/2 - dy*j;
	xf[i][j] = -r0/2 + dx*i + lambda*yf[i][j];
        glVertex2f(xf[i][j], yf[i][j]);
    }
  }
  */	

  progress(xf, yf);
}


// Main
int main(int argc, char **argv) {
  clock_t begin = clock();
  startCurve();
  clock_t end = clock();
  printf("Total computation time (s): %.10f\n", double(diffclock(end,begin)));

  return 0;

}


