/* main.cpp 
 * GPU Benchmark Immersed Boundary Structured Grid
 * Based on "The Method of Regularized Stokelets" by R.Cortez`
 * C.Copos 10/12/2012
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
const int N = 32; // # side square partitions 
const float visc = 1.0f;
const float flow = 1.0f;

// time 
const float TMAX = 500.0f;
const float tstep = 0.1f /*0.0001f*/;
static double stop = TMAX/tstep;

// curve constants
const float r0 = 1.0f;
const float ri = 1.2f; 
const double ds = ri/(N-1); const double e = ds; // e: parameter determining width of blobs or cuttofs

// vector structure
struct vector{
  double x; // x-component
  double y; // y-component
};

// spring constants
const double kl = 1.0f;

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

// Compute average elastic force at point j due to points j+1 and j-1
vector force(int j, double x[4*N-4], double y[4*N-4]) {
  vector fp; vector tp; // tangent vector for j to j+1
  vector fm; vector tm; // tangent vector for j-1 to j
  vector f; // final force
  double r = r0/(N-1); // rest length for principal springs
  double ks = kl; // spring constant

  // determine correct neighbors
  int j_jm, j_jp; // j_jm = j-1; j_jp = j+1
  if (j==0) { j_jm = (4*N-5); j_jp = j + 1; }
  else if (j == (4*N-5)) { j_jm = 4*N - 6; j_jp = 0;}
  else { j_jm = j - 1; j_jp = j + 1; }

  double dlkp = sqrtf(powf(x[j_jp]-x[j],2) + powf(y[j_jp]-y[j],2)); // distance between nodes j and j+1
  double dlkm = sqrtf(powf(x[j]-x[j_jm],2) + powf(y[j]-y[j_jm],2)); // distance between nodes j and j-1
  double dlk = sqrtf(powf(x[j_jp]-x[j_jm],2) + powf(y[j_jp]-y[j_jm],2));  // distance between nodes j-1 and j+1

  if (dlk < 0.00001f && dlk > -0.00001f) { f.x = 0.0f; f.y = 0.0f; } // make sure not dividing by zero
  else {
    // evaluate tangent vectors
    tp.x = (x[j]-x[j_jp])/dlkp;
    tp.y = (y[j]-y[j_jp])/dlkp;
    tm.x = (x[j]-x[j_jm])/dlkm;
    tm.y = (y[j]-y[j_jm])/dlkm;

    // evaluate helper functions   
    fp.x = -1.0f*ks*(dlkp-r)*(tp.x);
    fp.y = -1.0f*ks*(dlkp-r)*(tp.y); 
    fm.x = -1.0f*ks*(dlkm-r)*(tm.x);
    fm.y = -1.0f*ks*(dlkm-r)*(tm.y);   

    // evaluate final force
    f.x = (double)(fp.x + fm.x);
    f.y = (double)(fp.y + fm.y); 
  }

  //printf("Force @ position: (%f,%f) is: (%f,%f)\n", x[j], y[j], f.x, f.y);
  //printf("Partial forces at above location are x(+,-) = (%f,%f), y(+,-) = (%f,%f)\n", fp.x, fm.x, fp.y, fm.y); 

  return f;   
}

// Compute pressure at point i due to all other points in the structured grid
double pressure(int i, double x[4*N-4], double y[4*N-4]) {
  double p; vector f;
  double pp = 0.0f; // partial pressure
  double rk = 0.0f; double sq = 0.0f;

  for (int j=0; j<(4*N-4); j++) { // loop over all nodes in the grid
  	f.x = force(j, x, y).x;
        f.y = force(j, x, y).y;

	double r = sqrtf(powf(x[i],2) + powf(y[i],2));
	double pk = sqrtf(powf(x[j],2)+powf(y[j],2));
        double theta = atan2f(y[i], x[i]);
        double thetaj = atan2f(y[j], x[j]);
        
	double dtheta;
	if (theta>PI) { dtheta = thetaj + 2*PI - theta; }
        else { dtheta = theta - thetaj; }

	// dealing with rounding off errors     
        if (dtheta > -0.000001 && dtheta < 0.000001) { rk = 0.0f; }
        else { rk = sqrtf(powf(r,2) + powf(pk,2) - 2*r*pk*cosf(dtheta)); } // rk^2 = |x-xk|^2 = r^2 + pk^2 - 2*r*pk*cos(theta-thetak) where r=|x|, pk=|xk| 
        sq = sqrtf(powf(rk,2)+powf(e,2));

	double h = (1.0f)/(2.0f*PI) * (powf(rk,2)+2.0f*powf(e,2)+e*sq)/((sq+e)*powf(sq,1.5f));
	
	pp = (f.x*(x[i]-x[j]) + f.y*(y[i]-y[j])) * h;
	p +=pp;
  }
 
  return p; 
}

// Compute velocity at point i due to all other points in the structured grid
vector velocity(int i, double x[4*N-4], double y[4*N-4]) {
  vector v; vector fj;
  v.x = 0.0f ; v.y = 0.0f;
  double pvx = 0.0f; double pvy = 0.0f; // partial component velocities 
  double rk = 0.0f; double sq = 0.0f;

  for (int j=0; j<(4*N-4); j++) { // loop over all nodes in the grid
	fj.x = force(j, x, y).x;
	fj.y = force(j, x, y).y;

        double r = sqrt(powf(x[i],2) + powf(y[i],2));
        double pk = sqrt(powf(x[j],2) + powf(y[j],2));
	double theta = atan2f(y[i],x[i]);
	double thetaj = atan2(y[j],x[j]);
        double dtheta;
        if (theta>PI) { dtheta = thetaj + 2*PI - theta; }
        else { dtheta = theta - thetaj; }

        // dealing with rounding off errors     
        if (dtheta > -0.000001 && dtheta < 0.000001) { rk = 0.0f; }
        else { rk = sqrtf(powf(r,2) + powf(pk,2) - 2*r*pk*cosf(dtheta)); } // rk^2 = |x-xk|^2 = r^2 + pk^2 - 2*r*pk*cos(theta-thetak) where r=|x|, pk=|xk|
	 
    	sq = sqrtf(powf(rk,2)+powf(e,2));
    	pvx = 0.0f; pvy = 0.0f;
    	float p1 = (1.0f/(4.0f*visc*PI)) * (logf(sq+e)-(e*(sq+2*e))/(sq*(sq+e)));
    	float p2 = (1.0f/(4.0f*visc*PI)) * (sq+2*e)/(sq*powf(sq+e,2));

    	pvx = -1.0f*p1*fj.x + p2*(powf(x[i]-x[j],2)*fj.x + (x[i]-x[j])*(y[i]-y[j])*fj.y);
    	pvy = -1.0f*p1*fj.y + p2*((x[i]-x[j])*(y[i]-y[j])*fj.x + powf(y[i]-y[j],2)*fj.y);

    	v.x += pvx;
    	v.y += pvy;
 }

 return v;
}

// Compute polygonal area
double area(double x[4*N-4], double y[4*N-4]) {
   double A = 0.0f;
   int i1; int i2; // i1 = i; i2 = i+1

   for(int i=0; i<(4*N-4); i++) {
	i1 = i;
   	if(i1 == (4*N-5)) { i2 = 0; }
   	else i2 = i+1;
  	A += 0.5f * (x[i1]*y[i2]-x[i2]*y[i1]);
   }   

   return A;
}

// Progression
void progress(double x[4*N-4], double y[4*N-4]) {
  double u[4*N-4]; double v[4*N-4];
  vector vel; 
  double p[4*N-4];

  for (int t=0; t<stop; t++) {
	for(int i=0; i<(4*N-4); i++) {
	  vel = velocity(i, x, y);
	  u[i] = vel.x;
	  v[i] = vel.y;
	  p[i] = pressure(i, x, y);
	}

  	for(int i=0; i<(4*N-4); i++) {
 	   x[i] = x[i] + tstep*u[i];
           y[i] = y[i] + tstep*v[i];
  	}
  }

  // display final stage only 
  glColor3f(0.0f, 0.0f, 1.0f); // blue
  glBegin(GL_POINTS);
  for (int i=0; i<(4*N-4); i++) {
     	glVertex2f(x[i],y[i]);
  }
  glEnd();
  glFlush();
 
  // compute final area
  printf("Starting area: %f, Final area: %f\n", powf(ri,2), fabs(area(x,y)) );
}

// Draw starting configuration
void startCurve() {
  double xf[4*N-4]; double yf[4*N-4];
  glColor3f(1.0f, 1.0f, 1.0f); // white
  glBegin(GL_POINTS);

  for (int k=0; k<(4*N-4); k++) {
        double dx = ri/(N-1);
        if (k>=0 && k<(N-1)) { // bottom side
          xf[k] = -ri/2 + dx*k;
          yf[k] = -ri/2;
        }
        else if (k>=(N-1) && k<(2*N-2)) {
          xf[k] = -ri/2 + dx*(N-1);
          yf[k] = -ri/2 + dx*(k-(N-1));
        }
        else if (k>=(2*N-2) && k<(3*N-3)) {
          xf[k] = ri/2 - dx*(k-(2*N-2));
          yf[k] = -ri/2 + dx*(N-1);
        }
        else {
          xf[k] = -ri/2;
          yf[k] = ri/2 - dx*(k-(3*N-3));
        }
	glVertex2f(xf[k], yf[k]);
	yf[k] = -1.0f*yf[k]; // changing signs
  }


  glEnd();
  glFlush();

  progress(xf, yf);
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
   double diffms = (diffticks)/CLOCKS_PER_SEC;
   return diffms;
}

// This function is passed to glutDisplayFunc in order to display OpenGL contents on the window
void display(void) {
  clock_t begin = clock();
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  markers();
  startCurve();

  clock_t end = clock();
  printf("Computation time (s): %f\n", double(diffclock(end,begin)) );
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


