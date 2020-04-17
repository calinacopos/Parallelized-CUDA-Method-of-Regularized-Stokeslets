/* main.cpp 

 * Normal Forces on a circle with Fluid Immersed Boundary (Example 4a)
 * Based on "The Method of Regularized Stokelets" by R.Cortez`
 * C.Copos 7/19/2012
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

const double PI = 3.14159265358979323846264338;
const int N =  /*10*/ /*50*/ 100 /*200*/;
const int M = 200; // number of discretized points on line (x, 3/10)
const double visc = 1.0;
const double flow = 1.0;

// curve constants
const double r0 = 1.0;
const double e = (2.0*PI)/(16.0*N); // e: parameter determining width of blobs or cuttofs

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

// Force computed using Stokeslets
vector force(int k, double x[N], double y[N]) {
   vector f; 
   double theta = atan2f(y[k], x[k]);
   double ds = 2.0*PI*r0/N;
 
   f.x = 2.0*sinf(3*theta)*cosf(theta)*ds;
   f.y = 2.0*sinf(3*theta)*sinf(theta)*ds;

   return f;
}

// Pressure computed using exact solution
double pressureE(int j, double xl[M], double yl[M]) {
    double p;
    double theta = atan2f(yl[j], xl[j]);
    double r = sqrtf(powf(xl[j],2) + powf(yl[j],2));  

    if ((xl[j] >= 0.0) && (xl[j] < sqrtf(91)/10)) {  // r < 1
	p = (-1.0)*powf(r,3)*sinf(3*theta);
    }
    else if ((xl[j] > (sqrtf(91)/10)) && (xl[j] < 2.0)) {  // r > 1
	p = powf(r,-3)*sinf(3*theta);
    }
    else { p = 0.0; }

    return p;
}

// Pressure computed using Stokeslets
double pressureS(int j, double x[N], double y[N], double xl[M], double yl[M]) {
  double p = 0.0; vector fk; fk.x = 0.0; fk.y = 0.0;
  for(int k=0; k<N; k++) {
	fk.x = force(k, x, y).x;
	fk.y = force(k, x, y).y;

	double rk = sqrtf(powf(xl[j] - x[k],2) + powf(yl[j] - y[k],2));
        double sq = sqrtf(powf(rk,2) + powf(e,2));
	double h = (1.0/(2.0*PI)) * (powf(rk,2) + 2.0*powf(e,2) + e*sq)/(powf(sq,3) * (sq+e));
	double pp = (fk.x*(xl[j]-x[k]) + fk.y*(yl[j]-y[k])) * h;	

	//printf("pp = %f, h = %f, rk = %f, sq = %f, e = %f, force = %f\n", pp, h, rk, sq, e, (fk.x*(xl[j]-x[k]) + fk.y*(yl[j]-y[k])));

	p = p + pp; 
  }

  return p;
}

// Velocity computer using exact solution
vector velE(int j, double xl[M], double yl[M]) {
  vector velE; vector f;
  velE.x = 0.0; velE.y = 0.0;
  double theta = atan2f(yl[j], xl[j]);
  double r = sqrtf(powf(xl[j],2) + powf(yl[j],2));

  if ((xl[j] >= 0.0) && (xl[j] < sqrtf(91)/10)) {  // r < 1
        velE.x = (3.0/8.0)*powf(r,2)*sinf(2*theta) + (1.0/16.0)*powf(r,4)*sinf(4*theta) - (1.0/4.0)*powf(r,4)*sinf(2*theta);
	velE.y = (3.0/8.0)*powf(r,2)*cosf(2*theta) - (1.0/16.0)*powf(r,4)*cosf(4*theta) - (1.0/4.0)*powf(r,4)*cosf(2*theta);
  }
  else if ((xl[j] >= (sqrtf(91)/10)) && (xl[j] < 2.0)) {  // r >= 1
    	velE.x = (1.0/8.0)*powf(r,-2)*sinf(2*theta) - (3.0/16.0)*powf(r,-4)*sinf(4*theta) + (1.0/4.0)*powf(r,-2)*sinf(4*theta);
	velE.y = (1.0/8.0)*powf(r,-2)*cosf(2*theta) + (3.0/16.0)*powf(r,-4)*cosf(4*theta) - (1.0/4.0)*powf(r,-2)*cosf(4*theta);
  }
  else { velE.x = 0.0; velE.y = 0.0; }

  return velE;
}

// Velocity computed using Stokeslets
vector velS(int j, double x[N], double y[N], double xl[M], double yl[M]) {
  double pvx = 0.0; // partial x velocity
  double pvy = 0.0; // partial y velocity
  vector vel; vector fk;
  vel.x = 0.0; vel.y = 0.0;
  for(int k=0; k<N; k++) {
	fk.x = force(k, x, y).x;
	fk.y = force(k, x, y).y;

	double r = sqrtf(powf(xl[j],2)+powf(yl[j],2)); // assuming center @ origin
	double pk = sqrtf(powf(x[k],2)+powf(y[k],2));
	double theta = atan2f(yl[j], xl[j]);
        double thetak = atan2f(y[k], x[k]);
  	double dtheta, rk;
	if (theta>PI) { dtheta = thetak + 2*PI - theta; }
	else { dtheta = theta - thetak; }

        /* // dealing with rounding off errors	
	if (dtheta > -0.000001 && dtheta < 0.000001) { rk = 0.0f; }
	else { rk = sqrtf(powf(r,2) + powf(pk,2) - 2*r*pk*cosf(dtheta)); } // rk^2 = |x-xk|^2 = r^2 + pk^2 - 2*r*pk*cos(theta-thetak) where r=|x|, pk=|xk| 
 	*/
	rk = sqrtf(powf(xl[j]-x[k],2) + powf(yl[j]-y[k],2));	

	double sq = sqrtf(powf(rk,2)+powf(e,2));	
	pvx = 0.0; pvy = 0.0;
	double g = (1.0/(4*visc*PI)) * (logf(sq+e)-(e*(sq+2*e))/(sq*(sq+e)));
	double h = (1.0/(4*visc*PI)) * (sq+2*e)/(sq*powf(sq+e,2));

	pvx = -1.0*g*fk.x + h*(powf(xl[j]-x[k],2)*fk.x + (xl[j]-x[k])*(yl[j]-y[k])*fk.y);
	pvy = -1.0*g*fk.y + h*((xl[j]-x[k])*(yl[j]-y[k])*fk.x + powf(yl[j]-y[k],2)*fk.y);

	vel.x += pvx;
	vel.y += pvy;
	//if(j==0) printf("k = %d, x[%d] = %f, y[%d] = %f, pvx = %f, pvy = %f, e = %f, h = %f, sq = %f\n", k, k, x[k], k, y[k], pvx, pvy, e, h, sq);
  }
  return vel;
}

void drawCurve(double x[N], double y[N], double xl[M], double yl[M]) {
  double uE[M]; double vE[M];
  double uS[M]; double vS[M];
  vector velocityE; vector velocityS; 
  double pS[M]; double pE[M];
  double errorP[M]; double errorU[M]; double errorV[M];  
  double eP = 0.0; double eU = 0.0; double eV = 0.0; // L2 norm errors

  ofstream file_1, file_2, file_3, file_4, file_5, file_6, file_7, file_8, file_9;
  file_1.open("pressure_exact.txt");
  file_2.open("pressure_stokeslets.txt");
  file_3.open("u_exact.txt");
  file_4.open("v_exact.txt");
  file_5.open("u_stokeslets.txt");
  file_6.open("v_stokeslets.txt");
  file_7.open("pressure_error.txt");
  file_8.open("u_error.txt");
  file_9.open("v_error.txt");

  for (int j=0; j<M; j++) {
	// exact velocity
	velocityE = velE(j, xl, yl);
	uE[j] = velocityE.x;
        vE[j] = velocityE.y;
	// stokeslets velocity
	velocityS = velS(j, x, y, xl, yl);
	uS[j] = velocityS.x;
	vS[j] = velocityS.y;
	// exact pressure
	pE[j] = pressureE(j, xl, yl); 
	// stokeslets pressure
	pS[j] = pressureS(j, x, y, xl, yl);

 	// compute discretized errors
	errorP[j] = pS[j] - pE[j]; 
	errorU[j] = uS[j] - uE[j];
   	errorV[j] = vS[j] - vE[j];

 	// compute L2 norm erros
	eP += powf(errorP[j],2);
	eU += powf(errorU[j],2);
	eV += powf(errorV[j],2); 

	// save to files
	file_1 << xl[j] << " " << pE[j] << endl; 
	file_2 << xl[j] << " " << pS[j] << endl; 
  	file_3 << xl[j] << " " << uE[j] << endl;
        file_4 << xl[j] << " " << vE[j] << endl;
	file_5 << xl[j] << " " << uS[j] << endl;
	file_6 << xl[j] << " " << vS[j] << endl;
  	file_7 << xl[j] << " " << errorP[j] << endl;
	file_8 << xl[j] << " " << errorU[j] << endl;
	file_9 << xl[j] << " " << errorV[j] << endl;
  }

  printf("M = %d, Error in p = %.16f, Error in u = %.16f, Error in v =%.16f\n", M, sqrtf(eP), sqrtf(eU), sqrtf(eV));  

  file_1.close(); file_2.close(); file_3.close(); file_4.close(); file_5.close(); file_6.close(); file_7.close(); file_8.close(); file_9.close();
}

// Draw initial curve
void initialCurve() {
  double xi[N]; double yi[N]; // circle radius 1
  double xl[N]; double yl[N]; // line at (x, 3/10)

  glColor3f(1.0f, 1.0f, 1.0f); // white
  glBegin(GL_POINTS);
  for (int i=0; i<N; i++) {
        double theta = (2*PI/N)*i;
        xi[i] = r0*cosf(theta);
	yi[i] = r0*sinf(theta); 
        glVertex2f(xi[i],yi[i]);
  }
  for (int i=0; i<M; i++) {
	xl[i] = (2.0f/M)*i;
	yl[i] = 0.3f;
	glVertex2f(xl[i],yl[i]);
  }
  glEnd();
  glFlush();

  // draw curve through time
  drawCurve(xi, yi, xl, yl);
}

// Draw contour
void markers() {
  glColor3f(1.0f, 1.0f, 1.0f); // white
  glBegin(GL_LINES);
  	glVertex2f(-1.0f, -1.0f);
	glVertex2f(-1.0f, 1.0f);
	glVertex2f(-1.0f, 1.0f);
	glVertex2f(1.0f, 1.0f);
	glVertex2f(1.0f, 1.0f);
	glVertex2f(1.0f, -1.0f);
	glVertex2f(1.0f, -1.0f);
	glVertex2f(-1.0f, -1.0f);
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
  glutMainLoop();

  return 0;
}


