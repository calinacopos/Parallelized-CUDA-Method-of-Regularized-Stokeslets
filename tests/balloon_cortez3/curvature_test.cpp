#include <math.h>
#include <cstdio>
#include <cstdlib>
#include <stdlib.h>
#include <iostream>

// vector structure
struct vector{
  float x; // x-component
  float y; // y-component
};


const float a[5][2]={{0.0f,0.0f}, {0.5f,-1.0f}, {0.75f,-0.1f}, {1.0f,-1.0f}, {2.0f,0.0f}};


int main() {
  for(int i=0; i<5; i++) {
	float lkm, lkp;
	int p0, p1, p2;
	vector c;
	if(i==0) { p0=4; p1=i; p2=i+1; }
	else if(i==4) { p0=i-1; p1=i; p2=0; }
	else { p0=i-1; p1=i; p2=i+1; }
	
	lkp = sqrtf(powf(a[p2][0]-a[p1][0],2) + powf(a[p2][1]-a[p1][1],2));
	lkm = sqrtf(powf(a[p1][0]-a[p0][0],2) + powf(a[p1][1]-a[p0][1],2));
	
	c.x = 2*((a[p0][0]-a[p1][0])*lkp + (a[p2][0]-a[p1][0])*lkm)/(lkm*lkp*(lkm+lkp));
   	c.y = 2*((a[p0][1]-a[p1][1])*lkp + (a[p2][1]-a[p1][1])*lkm)/(lkm*lkp*(lkm+lkp));
	float kappa = sqrtf(powf(c.x,2)+powf(c.y,2));

	printf("x[%d] = %f, y[%d] = %f, k = %f, c.x = %f, c.y = %f\n", i, a[i][0], i, a[i][1], kappa, c.x, c.y);	
  } 

  return 0;
}


