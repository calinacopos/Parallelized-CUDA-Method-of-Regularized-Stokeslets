/* main.cpp 
 * Balloon with Fluid Immersed Boundary
 * Based on "The Method of Regularized Stokelets" by R.Cortez`
 * C.Copos 7/5/2012
 */

#include <GL/glut.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h> 
#include <ctime> 
#include <cstdio>
#include <cstdlib>

using namespace std;

const double PI = 3.14159265358979323846264338f;
const double RMAX = 0.297487f;
const double RMIN = 0.197487f;
const int N =  /*10*/ 50 /*100*/ /*200*/;
const double visc = 1.0f;
const double flow = 1.0f;

// time 
const int TMAX = 2;
const double tstep = 0.01f /*0.0001f*/;
static double stop = TMAX/tstep;
static int incr = 1;

// curve constants
const double a = 0.05f; const double r0 = 0.25f; const double b = sqrtf(powf(r0,2)-0.5*powf(a,2));
const double ds = /*2*PI*r0/N*/ /*0.1615965f*/ 0.0323193f /*0.0161693f*/ /*0.0080858f*/; const double e = 1.2f*ds; // e: parameter determining width of blobs or cuttofs

// vector structure
struct vector{
  double x; // x-component
  double y; // y-component
};

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

// Evaluation helper functions
// Finding max distance from origin to interface
double findRMax(double x[N], double y[N]) {
  double rmax;
  for (int i=0; i<N; i++) {
        double d = sqrtf(powf(x[i],2)+powf(y[i],2));
        if (d > rmax || d == rmax) rmax = d;
  }
  return rmax;
}

// Finding min distance from origin to interface
double findRMin(double x[N], double y[N]) {
  double rmin;
  for (int i=0; i<N; i++) {
        double d = sqrtf(powf(x[i],2)+powf(y[i],2));
        if (d < rmin || d == rmin) rmin = d;
  }
  return rmin;
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

   lkp = sqrtf(powf(x[p2]-x[p1],2) + powf(y[p2]-y[p1],2));
   lkm = sqrtf(powf(x[p1]-x[p0],2) + powf(y[p1]-y[p0],2));
   lk = 1.0f;
   c.x = ((x[p2]-x[p1])/lkp-(x[p1]-x[p0])/lkm)/lk;
   c.y = ((y[p2]-y[p1])/lkp-(y[p1]-y[p0])/lkm)/lk;
   kappa = sqrtf(powf(c.x,2)+powf(c.y,2));

   // alternative
   double lk1, lk2, lk3;
   lk1 = lkp;// k+1 to k
   lk2 = lkm;// k to k-1
   lk3 = sqrtf(powf(x[p2]-x[p0],2) + powf(y[p2]-y[p0],2));// k+1 to k-1
   double s = 0.5*(lk1+lk2+lk3); //semi-perimeter
   double a = sqrtf(s*(s-lk1)*(s-lk2)*(s-lk3)); //area
   double kappa_t1 = (4*a)/(fabs(lk1)*fabs(lk2)*fabs(lk3)); 
   
   vector c_t2; // curvature trial 2
   c_t2.x = 2*((x[p0]-x[p1])*lkp + (x[p2]-x[p1])*lkm)/(lkm*lkp*(lkm+lkp));
   c_t2.y = 2*((y[p0]-y[p1])*lkp + (y[p2]-y[p1])*lkm)/(lkm*lkp*(lkm+lkp));
   double kappa_t2 = sqrtf(powf(c_t2.x,2)+powf(c_t2.y,2));

   f.x = (double)(flow*kappa_t2 - 1.0f/r0)*(c_t2.x/kappa_t2);
   f.y = (double)(flow*kappa_t2 - 1.0f/r0)*(c_t2.y/kappa_t2);

   // dealing with rounding-off errors
   if (f.x > -0.000001 && f.x < 0.000001) f.x = 0.0f;
   if (f.y > -0.000001 && f.y < 0.000001) f.y = 0.0f;

   //printf("k = %d, lkp = %f, lkm = %f, kappa_t1 = %f, kappa_t2 = %f, c.x = %f, c.y = %f, f.x = %f, f.y = %f\n", k, lkp, lkm, kappa_t2, kappa_t1, c_t2.x, c_t2.y, f.x, f.y);

  return f;
}

// pressure
double pressure(int j, double x[N], double y[N]) {
  double pp = 0.0f; // partial pressure
  double p; vector fk;
  for(int k=0; k<N; k++) {
	pp = 0.0f;
  	if (x[k] > -0.000001 && x[k] < 0.000001) x[k] = 0.0f;
        if (x[j] > -0.000001 && x[j] < 0.000001) x[j] = 0.0f;
        if (y[k] > -0.000001 && y[k] < 0.000001) y[k] = 0.0f;
        if (y[j] > -0.000001 && y[j] < 0.000001) y[j] = 0.0f;
	fk.x = force(k, x, y).x;
	fk.y = force(k, x, y).y;
        
	double r = sqrtf(powf(x[j],2)+powf(y[j],2)); // assuming center @ origin
        double pk = sqrtf(powf(x[k],2)+powf(y[k],2));
        double theta = atan2f(y[j], x[j]);
        double thetak = atan2f(y[k], x[k]);
        double dtheta, rk;
        if (theta>PI) { dtheta = thetak + 2*PI - theta; }
        else { dtheta = theta - thetak; }

        // dealing with rounding off errors     
        if (dtheta > -0.000001 && dtheta < 0.000001) { rk = 0.0f; }
        else { rk = sqrtf(powf(r,2) + powf(pk,2) - 2*r*pk*cosf(dtheta)); } // rk^2 = |x-xk|^2 = r^2 + pk^2 - 2*r*pk*cos(theta-thetak) where r=|x|, pk=|xk| 

        double sq = sqrtf(powf(rk,2)+powf(e,2));
	double h = (1/(2*PI)) * (powf(rk,2)+2*powf(e,2)+e*sq)/(powf(sq,3)*(sq+e));

	pp = (fk.x*(x[j]-x[k]) + fk.y*(y[j]-y[k])) * h;	
	p += pp; 
  }
}

// velocity
vector vel(int j, double x[N], double y[N]) {
  double pvx = 0.0f; // partial x velocity
  double pvy = 0.0f; // partial y velocity
  vector vel; vector fk;
  vel.x = 0.0f; vel.y = 0.0f;
  for(int k=0; k<N; k++) {
        if (x[k] > -0.000001 && x[k] < 0.000001) x[k] = 0.0f;
        if (x[j] > -0.000001 && x[j] < 0.000001) x[j] = 0.0f;
        if (y[k] > -0.000001 && y[k] < 0.000001) y[k] = 0.0f;
        if (y[j] > -0.000001 && y[j] < 0.000001) y[j] = 0.0f;
	fk.x = force(k, x, y).x;
	fk.y = force(k, x, y).y;

 	//printf("k = %d, x[%d] = %f, y[%d] = %f, f.x = %f, f.y = %f\n", k, k, x[k], k, y[k], fk.x, fk.y);

	double r = sqrtf(powf(x[j],2)+powf(y[j],2)); // assuming center @ origin
	double pk = sqrtf(powf(x[k],2)+powf(y[k],2));
	double theta = atan2f(y[j], x[j]);
        double thetak = atan2f(y[k], x[k]);
  	double dtheta, rk;
	if (theta>PI) { dtheta = thetak + 2*PI - theta; }
	else { dtheta = theta - thetak; }

        // dealing with rounding off errors	
	if (dtheta > -0.000001 && dtheta < 0.000001) { rk = 0.0f; }
	else { rk = sqrtf(powf(r,2) + powf(pk,2) - 2*r*pk*cosf(dtheta)); } // rk^2 = |x-xk|^2 = r^2 + pk^2 - 2*r*pk*cos(theta-thetak) where r=|x|, pk=|xk| 
	
	double sq = sqrtf(powf(rk,2)+powf(e,2));	
	pvx = 0.0f; pvy = 0.0f;
	double g = (1.0f/(4*visc*PI)) * (logf(sq+e)-(e*(sq+2*e))/(sq*(sq+e)));
	double h = (1.0f/(4*visc*PI)) * (sq+2*e)/(sq*powf(sq+e,2));

	pvx = -g*fk.x + h*(powf(x[j]-x[k],2)*fk.x + (x[j]-x[k])*(y[j]-y[k])*fk.y);
	pvy = -g*fk.y + h*((x[j]-x[k])*(y[j]-y[k])*fk.x + powf(y[j]-y[k],2)*fk.y);

	vel.x += pvx;
	vel.y += pvy;
	//printf("k = %d, x[%d] = %f, y[%d] = %f, vel.x = %f, vel.y = %f\n", k, k, x[k], k, y[k], vel.x, vel.y);
  }
  return vel;
}

void drawCurve(double x[N], double y[N]) {
  double ux[N]; double uy[N];
  vector velocity; double p[N];
  ofstream file;
  file.open("pressure_output_t0.txt");

  for (int t=0; t<stop; t+=incr) {
	// compute velocities
  	for (int j=0; j<N; j++) {
	   velocity = vel(j, x, y);
	   ux[j] = velocity.x;
           uy[j] = velocity.y;
	   p[j] = pressure(j, x, y);
	   if (t==0) file << x[j] << " " << y[j] <<  " " << ux[j] << " " << uy[j] << " " << p[j] <<  endl;
	}

	// compute positions
	for(int j=0; j<N; j++) {
	   x[j] = x[j] + tstep*ux[j];
           y[j] = y[j] + tstep*uy[j];
	}
  }

  // print final configuration
  glColor3f(0.0f, 0.0f, 1.0f); // blue
  glBegin(GL_POINTS);
  for (int j=0; j<N; j++) {
	glVertex2f(x[j],y[j]);
  }
  glEnd();
  glFlush();

  // Verify the shape of the final solution and check convergence as t -> infty
  printf("FINAL Rmax= %f, Rmin= %f, Diff= %f\n", findRMax(x,y), findRMin(x,y), fabs(findRMax(x,y)-findRMin(x,y)));
}

// Draw initial curve (ellipse:  r(theta) = b + a*cos(2*theta) )
void initialCurve() {
  double xi[N]; double yi[N];
  // placing evenly spaced points for initial curve at t = 0
  glColor3f(1.0f, 1.0f, 1.0f); // white
  glBegin(GL_POINTS);
  for (int i=0; i<N; i++) {
        double theta = (2*PI/N)*i;
        xi[i] = (b+a*cos(2*theta))*cos(theta);
        yi[i] = (b+a*cos(2*theta))*sin(theta);
        //xi[i] = r0*cos(theta); // DEBUG (circle)
	//yi[i] = r0*sin(theta); // DEBUG (circle)
        glVertex2f(xi[i],yi[i]);
  }
  glEnd();
  glFlush();

  // draw curve through time
  drawCurve(xi, yi);
}

// Draw contour
void markers() {
  glColor3f(1.0f, 1.0f, 1.0f); // white
  glBegin(GL_LINES);
  	glVertex2f(-0.5f, -0.5f);
	glVertex2f(-0.5f, 0.5f);
	glVertex2f(-0.5f, 0.5f);
	glVertex2f(0.5f, 0.5f);
	glVertex2f(0.5f, 0.5f);
	glVertex2f(0.5f, -0.5f);
	glVertex2f(0.5f, -0.5f);
	glVertex2f(-0.5f, -0.5f);
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

// This function is passed to glutDisplayFunc in order to display OpenGL contents on the window
void display(void) {
  clock_t start, end;
  start = clock();

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  markers();
  initialCurve();

  end = clock();
  //printf("Computation time (seconds): %f\n", end-start);
  glutSwapBuffers();
}

void timer_function (int value) {
  static bool flag = true;
  static int count = 5;
  if (count == 10) { flag = false; }
  if (flag) { incr = 1; count++; }
  else { incr = -1; count--; }
  glutPostRedisplay();
  glutTimerFunc(200, timer_function, 0);
}

// Main
int main(int argc, char **argv) {
  // GLUT Initialiation
  glutInit(&argc, argv);
  glutInitWindowSize(window_width, window_height);
  glutInitWindowPosition(window_x, window_y);
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
  glutCreateWindow("Immersed Boundary using Stokelets");
 
  init();

  // Callback functions
  glutDisplayFunc(display);
  glutKeyboardFunc(keyboard);
  glutTimerFunc(1, timer_function, 0);
  glutMainLoop();

  return 0;
}


