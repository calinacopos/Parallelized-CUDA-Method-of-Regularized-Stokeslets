/* main.cpp
 *
 * 
 * Unsteady motion of an interface immersed in viscous fluid
 * Based on "The Method of Regularized Stokelets" by R.Cortez` Example 3
 *
 * C.Copos 6/5/2012
 */

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cmath> 
#include <ctime> 
#include <cstdio>
#include <cstdlib>

using namespace std;

const double PI = 3.14159265358979323846264338;
const int N =  /*10*/ 200 /*100*/ /*200*/;
const double visc = 1.0;
const double flow = 1.0;

// time 
const int TMAX = 1;
const double tstep = 0.0001 /*0.0001f*/;
static int stop = TMAX/tstep;
static int incr = 10;

// curve constants
const double a = 0.05; 
const double r0 = 0.25; 
const double b = sqrt(pow(r0,2)-0.5*pow(a,2));
const double ds = /*2*PI*r0/N*/ /*0.1615965f*/ /*0.0323193*/ /*0.0161693*/ 0.0080858; 
const double e = 1.2*ds; // e: parameter determining width of blobs or cuttofs

// vector structure
struct vector{
  double x; // x-component
  double y; // y-component
};

// Evaluation helper functions
// Finding max distance from origin to interface
double findRMax(double x[N], double y[N]) {
  double rmax = 0.0;
  for (int i=0; i<N; i++) {
        double d = sqrt(pow(x[i],2)+pow(y[i],2));
        if (d > rmax) rmax = d;
  }
  return rmax;
}

// Finding min distance from origin to interface
double findRMin(double x[N], double y[N]) {
  double rmin = r0+a;
  for (int i=0; i<N; i++) {
        double d = sqrt(pow(x[i],2)+pow(y[i],2));
        if (d < rmin) rmin = d;
  }
  return rmin;
}

// Finding area
double findArea(double x[N], double y[N]) {
  double r = findRMax(x,y); 
  double A = PI*pow(r,2);
  return A;
}

vector force(int k, double x[N], double y[N]) {
   int p0, p1, p2;
   if(k==0) { p0=N-1; p1=k; p2=k+1;}
   else if (k==(N-1)) { p0=k-1; p1=k; p2=0;}
   else { p0=k-1; p1=k; p2=k+1;}

   vector c; // curvature
   vector f; // force
   double lkp, lkm, lk;
   double kappa; // magnitude of curvature  

   lkp = sqrt(pow(x[p2]-x[p1],2) + pow(y[p2]-y[p1],2));
   lkm = sqrt(pow(x[p1]-x[p0],2) + pow(y[p1]-y[p0],2));
   lk = 1.0;
   c.x = ((x[p2]-x[p1])/lkp-(x[p1]-x[p0])/lkm)/lk;
   c.y = ((y[p2]-y[p1])/lkp-(y[p1]-y[p0])/lkm)/lk;
   kappa = sqrt(pow(c.x,2)+pow(c.y,2));

   // alternative
   double lk1, lk2, lk3;
   lk1 = lkp;// k+1 to k
   lk2 = lkm;// k to k-1
   lk3 = sqrt(pow(x[p2]-x[p0],2) + pow(y[p2]-y[p0],2));// k+1 to k-1
   double s = 0.5*(lk1+lk2+lk3); //semi-perimeter
   double a = sqrt(s*(s-lk1)*(s-lk2)*(s-lk3)); //area
   double kappa_t1 = (4*a)/(fabs(lk1)*fabs(lk2)*fabs(lk3)); 
   
   vector c_t2; // curvature trial 2
   c_t2.x = 2*((x[p0]-x[p1])*lkp + (x[p2]-x[p1])*lkm)/(lkm*lkp*(lkm+lkp));
   c_t2.y = 2*((y[p0]-y[p1])*lkp + (y[p2]-y[p1])*lkm)/(lkm*lkp*(lkm+lkp));
   double kappa_t2 = sqrt(pow(c_t2.x,2)+pow(c_t2.y,2));

   f.x = (double)(flow*kappa_t2 - 1.0/r0)*(c_t2.x/kappa_t2);
   f.y = (double)(flow*kappa_t2 - 1.0/r0)*(c_t2.y/kappa_t2);

   // dealing with rounding-off errors
   if (f.x > -0.000001 && f.x < 0.000001) f.x = 0.0;
   if (f.y > -0.000001 && f.y < 0.000001) f.y = 0.0;

   return f;
}

// pressure
double pressure(int j, double x[N], double y[N]) {
  double pp = 0.0; // partial pressure
  double p; vector fk;
  for(int k=0; k<N; k++) {
	pp = 0.0;
  	if (x[k] > -0.000001 && x[k] < 0.000001) x[k] = 0.0;
        if (x[j] > -0.000001 && x[j] < 0.000001) x[j] = 0.0;
        if (y[k] > -0.000001 && y[k] < 0.000001) y[k] = 0.0;
        if (y[j] > -0.000001 && y[j] < 0.000001) y[j] = 0.0;
	fk.x = force(k, x, y).x;
	fk.y = force(k, x, y).y;
        
	double r = sqrt(pow(x[j],2)+pow(y[j],2)); // assuming center @ origin
        double pk = sqrt(pow(x[k],2)+pow(y[k],2));
        double theta = atan2f(y[j], x[j]);
        double thetak = atan2f(y[k], x[k]);
        double dtheta, rk;
        if (theta>PI) { dtheta = thetak + 2*PI - theta; }
        else { dtheta = theta - thetak; }

        // dealing with rounding off errors     
        if (dtheta > -0.000001 && dtheta < 0.000001) { rk = 0.0; }
        else { rk = sqrt(pow(r,2) + pow(pk,2) - 2*r*pk*cosf(dtheta)); } // rk^2 = |x-xk|^2 = r^2 + pk^2 - 2*r*pk*cos(theta-thetak) where r=|x|, pk=|xk| 

        double sq = sqrt(pow(rk,2)+pow(e,2));
	double h = (1/(2*PI)) * (pow(rk,2)+2*pow(e,2)+e*sq)/(pow(sq,3)*(sq+e));

	pp = (fk.x*(x[j]-x[k]) + fk.y*(y[j]-y[k])) * h;	
	p += pp; 
  }
}

// velocity
vector vel(int j, double x[N], double y[N]) {
  double pvx = 0.0; // partial x velocity
  double pvy = 0.0; // partial y velocity
  vector vel; vector fk;
  vel.x = 0.0; vel.y = 0.0;
  for(int k=0; k<N; k++) {
        //if (x[k] > -0.000001 && x[k] < 0.000001) x[k] = 0.0;
        //if (x[j] > -0.000001 && x[j] < 0.000001) x[j] = 0.0;
        //if (y[k] > -0.000001 && y[k] < 0.000001) y[k] = 0.0;
        //if (y[j] > -0.000001 && y[j] < 0.000001) y[j] = 0.0;
	fk.x = force(k, x, y).x;
	fk.y = force(k, x, y).y;

 	//printf("k = %d, x[%d] = %f, y[%d] = %f, f.x = %f, f.y = %f\n", k, k, x[k], k, y[k], fk.x, fk.y);

	// OLD WAY OF COMPUTING
	//double r = sqrt(pow(x[j],2)+pow(y[j],2)); // assuming center @ origin
	//double pk = sqrt(pow(x[k],2)+pow(y[k],2));
	//double theta = atan2f(y[j], x[j]);
        //double thetak = atan2f(y[k], x[k]);
  	//double dtheta, rk;
	//if (theta>PI) { dtheta = thetak + 2*PI - theta; }
	//else { dtheta = theta - thetak; }

        // dealing with rounding off errors	
	//if (dtheta > -0.000001 && dtheta < 0.000001) { rk = 0.0; }
	//else { rk = sqrt(pow(r,2) + pow(pk,2) - 2*r*pk*cosf(dtheta)); } // rk^2 = |x-xk|^2 = r^2 + pk^2 - 2*r*pk*cos(theta-thetak) where r=|x|, pk=|xk| 
	
	
	double rk = sqrt(powf(x[j]-x[k],2) + powf(y[j]-y[k],2));
	double sq = sqrt(pow(rk,2)+pow(e,2));	
	pvx = 0.0; pvy = 0.0;
	double g = (1.0/(4.0*visc*PI)) * (logf(sq+e)-(e*(sq+2*e))/(sq*(sq+e)));
	double h = (1.0/(4.0*visc*PI)) * (sq+2*e)/(sq*pow(sq+e,2));

	pvx = -g*fk.x + h*(pow(x[j]-x[k],2)*fk.x + (x[j]-x[k])*(y[j]-y[k])*fk.y);
	pvy = -g*fk.y + h*((x[j]-x[k])*(y[j]-y[k])*fk.x + pow(y[j]-y[k],2)*fk.y);

	vel.x += pvx;
	vel.y += pvy;
  }
  return vel;
}

void drawCurve(double x[N], double y[N]) {
  double ux[N]; 
  double uy[N];
  vector velocity; 
  double p[N];
  ofstream f1, f2, f3;
  f1.open("initial_conf.txt");
  f2.open("final_conf.txt");
  f3.open("radius.txt");

  for (int t=0; t<stop; t+=incr) {
	// report max and min radius over time
        f3 << findRMax(x,y) << " " << findRMin(x,y) << endl;

	// compute velocities
  	for (int j=0; j<N; j++) {
	   velocity = vel(j, x, y);
	   ux[j] = velocity.x;
           uy[j] = velocity.y;
	   p[j] = pressure(j, x, y);
	   if (t==0.0) f1 << x[j] << " " << y[j] << endl;
	}

	// compute positions
	for(int j=0; j<N; j++) {
	   x[j] = x[j] + tstep*ux[j];
           y[j] = y[j] + tstep*uy[j];
	}
  }

  for (int j=0; j<N; j++) {
	f2 << x[j] << " " << y[j] << endl;
  }

  // Verify the shape of the final solution and check convergence as t -> infty
  double Diff_r = fabs(findRMax(x,y)-findRMin(x,y));
  double Diff_area = fabs(findArea(x,y) - PI*pow(r0,2));
  printf("N= %d ds= %.16f Diff_r= %.16f Diff_area= %.16f \n", N, ds, Diff_r, Diff_area);
}

// Main
int main(int argc, char **argv) {
  double xi[N]; double yi[N];
  // placing evenly spaced points for initial curve at t = 0
  for (int i=0; i<N; i++) {
        double theta = (2*PI/N)*i;
        xi[i] = (b+a*cos(2*theta))*cos(theta);
        yi[i] = (b+a*cos(2*theta))*sin(theta);
        //xi[i] = r0*cos(theta); // DEBUG (circle)
        //yi[i] = r0*sin(theta); // DEBUG (circle)
  }

  drawCurve(xi,yi);

  return 0;
}


