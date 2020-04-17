/* main.cpp 
 * GPU Benchmark Immersed Boundary Structured Grid
 * Based on "The Method of Regularized Stokelets" by R.Cortez`
 * C.Copos 7/24/2012
 */

#include <GL/glut.h>
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
const int N = 32; // side of cube discretized
const float visc = 1.0f;
const float flow = 1.0f;

// time 
const int TMAX = 1000;
const float tstep = 0.1f /*0.0001f*/;
static double stop = TMAX/tstep;
static int incr = 1;

// curve constants
const float r0 = 1.0f;
const float ri = 1.5f; 
const double ds = ri/N; const double e = ds; // e: parameter determining width of blobs or cuttofs

// vector structure
struct vector{
  double x; // x-component
  double y; // y-component
};

// spring constants
const double kl = 0.1f;
const double kd = 0.5f*kl;

// OpenGL Initialization
void init();

// OpenGL Callback functions
void display(void);
void keyboard(unsigned char key);

// OpenGL Support Functions
void drawObject();

// Define the window position on screen
int window_x;
int window_y;

// Variable representing the window size and title
int window_width = 400;
int window_height = 400;

// Set OpenGL program initial state
void init(void) {
  glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
  glPointSize(2.0);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glEnable(GL_POINT_SMOOTH);
  glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
}

/* EVALUATION HELPER FUNCTIONS */
// Compute diagonal
double findDiag(double x[N][N], double y[N][N]) {
  double diag;
  for (int j=0; j<N; j++) {
	for (int i=0; i<N; i++) {
           double d = sqrtf(powf(x[i][j]-x[0][0],2)+powf(y[i][j]-y[0][0],2));
           if (d > diag || d == diag) diag = d;
  	}
  }
  return diag;
}

// Compute short side
double findMin(double x[N][N], double y[N][N]) {
  double rmin;
  double d1 = sqrtf(powf(x[N-1][0]-x[0][0],2) + powf(y[N-1][0]-y[0][0],2));
  double d2 = sqrtf(powf(x[0][N-1]-x[0][0],2) + powf(y[0][N-1]-y[0][0],2));

  if (d1 < d2) { rmin = d1; } 
  else { rmin = d2; }

  return rmin;
}

// Compute long side
double findMax(double x[N][N], double y[N][N]) {
  double rmax;
  double d1 = sqrtf(powf(x[N-1][0]-x[0][0],2) + powf(y[N-1][0]-y[0][0],2));
  double d2 = sqrtf(powf(x[0][N-1]-x[0][0],2) + powf(y[0][N-1]-y[0][0],2));

  if (d1 < d2) { rmax = d2; }
  else { rmax = d1; }

  return rmax;
}

double totalDiff(double x[N][N], double y[N][N]) {
  double xdiff = 0.0f; double ydiff = 0.0f;
  
  for (int j=0; j<N; j++) {
	double dy = 1.0f/(N-1);
	for (int i=0; i<N; i++) {
	   double dx = 1.0f/(N-1);
	   xdiff += fabs(x[i][j]-(-r0/2 + dx*i));
	   ydiff += fabs(y[i][j]-(-r0/2 + dy*j));
	}
  }

  return (sqrtf(powf(xdiff,2) + powf(ydiff,2)));
}
/* --------------------- */

// Compute force at point ij due to point kl
vector force(int i, int j, int k, int l, double x[N][N], double y[N][N]) {
  vector f; vector t; // tangent vector
  double rl = 1.0f/(N-1); // rest length for principal springs
  double rd = sqrtf(2)*rl; // rest length for diagonal springs
  double r; // rest lenght  
  double ks = 0.0f; // spring constant
  double d = sqrtf(powf(x[i][j]-x[k][l],2) + powf(y[i][j]-y[k][l],2));

  if (d < 0.00001f && d > -0.00001f) { f.x = 0.0f; f.y = 0.0f; } // make sure not dividing by zero
  else {
    // determine which spring constant to use
    if (i==k || j==l) { ks = kl; r = rl; }
    else { ks = kd; r = rd; }

    // evaluate tangent vector
    t.x = (x[i][j]-x[k][l])/d;
    t.y = (y[i][j]-y[k][l])/d;

    // evaluate helper functions
    f.x = ks*((d-r)/r)*t.x;
    f.y = ks*((d-r)/r)*t.y;
  }

  //printf("x[%d][%d] = %f, x[%d][%d] = %f, y[%d][%d] = %f, y[%d][%d] = %f, d = %f, f.x = %f, f.y = %f\n", i, j, x[i][j], k, l, x[k][l], i, j, y[i][j], k, l, y[k][l], d, f.x, f.y);

  return f;   
}

// Compute velocity at point ij due to all other points in the structured grid
vector velocity(int i, int j, double x[N][N], double y[N][N]) {
  vector v; vector f;
  v.x = 0.0f ; v.y = 0.0f;
  double pvx = 0.0f; double pvy = 0.0f; // partial component velocities  

  // determine neighbors
  int kstart; int kend; int lstart; int lend; int inner;
  if (i==0 && j==(N-1)) { 
     kstart = i; lstart = j-1;
     kend = i+1; lend = j;
  }
  else if (i==0 && j==0) {
     kstart = i; lstart = j;
     kend = i+1; lend = j+1;
  }
  else if (i==(N-1) && j==0) {
     kstart = i-1; lstart = j;
     kend = i; lend = j+1;
  }
  else if (i==(N-1) && j==(N-1)) {
     kstart = i-1; lstart = j-1;
     kend = i; lend = j;
  }
  else if (j==0) {
     kstart = i-1; lstart = j;
     kend = i+1; lend = j+1;
  }
  else if (i==0) {
     kstart = i; lstart = j-1;
     kend = i+1; lend = j+1;
  }
  else if (j==(N-1)) {
     kstart = i-1; lstart = j-1;
     kend = i+1; lend = j;
  }
  else if (i==(N-1)) {
     kstart = i-1; lstart = j-1;
     kend = i; lend = j+1;
  }
  else {
     kstart = i-1; lstart = j-1;
     kend = i+1; lend = j+1;
  }


  // Loop through neighbors and add up their contributions due to the elastic force on the velocity
  for (int l=lstart; l<=lend; l++) {
	for (int k=kstart; k<=kend; k++) {
	   double rk = 0.0f; double sq = 0.0f; double p = 0.0f; double r = 0.0f;
	   
	   if (k == i && l == j) { f.x = 0.0f; f.y = 0.0f; }
	   else { 
	     f.x = -1.0f * force(i, j, k, l, x, y).x;
	     f.y = -1.0f * force(i, j, k, l, x, y).y;
	 
	     //printf("k = %d, l = %d, x[%d][%d] = %f, y[%d][%d] = %f, f.x = %f, f.y = %f\n", k, l, i, j, x[i][j], i, j, y[i][j], f.x, f.y);

	     rk = sqrtf(powf(x[i][j]-x[k][l],2) + powf(y[i][j]-y[k][l],2));
	     sq = sqrtf(powf(rk,2)+powf(e,2));
	     pvx = 0.0f; pvy = 0.0f;
	     p = (1.0f/(4.0f*visc*PI)) * (logf(sq+e)-(e*(sq+2*e))/(sq*(sq+e)));
             r = (1.0f/(4.0f*visc*PI)) * (sq+2*e)/(sq*powf(sq+e,2));
	   }

	   if (inner == 1) { pvx = 0.0f; pvy = 0.0f; }
	   else {
             pvx = -1.0f*p*f.x + r*(powf(x[i][j]-x[k][l],2)*f.x + (x[i][j]-x[k][l])*(y[i][j]-y[k][l])*f.y);
             pvy = -1.0f*p*f.y + r*((x[i][j]-x[k][l])*(y[i][j]-y[k][l])*f.x + powf(y[i][j]-y[k][l],2)*f.y);

	     //printf("rk = %f, e = %f, sq = %f, p = %f, r = %f\n", rk, e, sq, p, r);
	     //printf("i = %d, j = %d, k = %d, l = %d, pvx = %f, pvy = %f\n", i, j, k, l, pvx, pvy);
	     //printf("k = %d, l = %d, x[%d][%d] = %f, y[%d][%d] = %f, f.x = %f, f.y = %f\n", k, l, i, j, x[i][j], i, j, y[i][j], f.x, f.y);

             v.x += pvx;
             v.y += pvy;
	   }
	}
 }

 return v;
}

// Progression
void progress(double x[N][N], double y[N][N]) {
  double u[N][N]; double v[N][N];
  vector vel;

  for (int t=0; t<stop; t+=incr) {
  	for(int j=0; j<N; j++) {
		for(int i=0; i<N; i++) {
		   vel = velocity(i, j, x, y);
		   u[i][j] = vel.x;
		   v[i][j] = vel.y;
		}
	}

	for(int j=0; j<N; j++) {
                for(int i=0; i<N; i++) {
 		   x[i][j] = x[i][j] + tstep*u[i][j];
        	   y[i][j] = y[i][j] + tstep*v[i][j];
        	   //printf("t = %d, x[%d][%d] = %f, y[%d][%d] = %f, u = %f, v = %f\n", t, i, j, x[i][j], i, j, y[i][j], u[i][j], v[i][j]);
  		}
	}
  }

  // display final stage only 
  glColor3f(0.0f, 0.0f, 1.0f); // blue
  glBegin(GL_POINTS);
  for (int j=0; j<N; j++) {
	for (int i=0; i<N; i++) {
           glVertex2f(x[i][j],y[i][j]);
	}
  }
  glEnd();
  glFlush();

  printf("FINAL Diag = %f, Smax= %f, Smin= %f, Diff= %f, Total diff = %f\n", findDiag(x,y), findMax(x,y), findMin(x,y), fabs(findMax(x,y)-findMin(x,y)), totalDiff(x,y));

}

// Draw initial configuration
void initialCurve() {
  double xi[N][N]; double yi[N][N];

  glColor3f(1.0f, 1.0f, 1.0f); // white
  glBegin(GL_POINTS);
  for (int j=0; j<N; j++) {
        double dy = ri/(N-1);
        for (int i=0; i<N; i++) {
          double dx = ri/(N-1);
          xi[i][j] = -ri/2 + dx*i;
          yi[i][j] = -ri/2 + dy*j;
          glVertex2f(xi[i][j], yi[i][j]);
        }
  }
  glEnd();
  glFlush();

  printf("INITIAL Diag = %f, Smax= %f, Smin= %f, Diff= %f\n", findDiag(xi,yi), findMax(xi,yi), findMin(xi,yi), fabs(findMax(xi,yi)-findMin(xi,yi)));

  progress(xi, yi);
}

// Draw contour
void markers() {
  glColor3f(1.0f, 1.0f, 1.0f); // white
  glBegin(GL_LINES);
  	glVertex2f(-0.8f, -0.8f);
	glVertex2f(-0.8f, 0.8f);
	glVertex2f(-0.8f, 0.8f);
	glVertex2f(0.8f, 0.8f);
	glVertex2f(0.8f, 0.8f);
	glVertex2f(0.8f, -0.8f);
	glVertex2f(0.8f, -0.8f);
	glVertex2f(-0.8f, -0.8f);
  glEnd();
  glFlush();
}

// This function is passed to the glutKeyboardFunc and is called whenever the user hits a key.
void keyboard(unsigned char key, int x, int y)
{
  if (key == 27) { // ESC (exit)
        exit(1);
  }

  glutPostRedisplay();
}

double diffclock(clock_t s, clock_t e) {
   double diffticks = s-e;
   double diffms = (diffticks*10.0f)/CLOCKS_PER_SEC;
   return diffms;
}

// This function is passed to glutDisplayFunc in order to display OpenGL contents on the window
void display(void) {
  clock_t begin = clock();
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  markers();
  initialCurve();

  clock_t end = clock();
  printf("Computation time (ms): %f\n", double(diffclock(end,begin)) );
  glutSwapBuffers();
}

// Main
int main(int argc, char **argv) {
  // GLUT Initialiation
  glutInit(&argc, argv);
  glutInitWindowSize(window_width, window_height);
  glutInitWindowPosition(window_x, window_y);
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
  glutCreateWindow("Immersed Structured Grid using Stokeslets");

  init();

  // Callback functions
  glutDisplayFunc(display);
  glutKeyboardFunc(keyboard);
  glutMainLoop();

  return 0;
}


