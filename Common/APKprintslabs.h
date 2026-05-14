#ifndef __APKprintslabs_h
#define __APKprintslabs_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "APKData_Structures.h"

int printslab(Base *);
int printslabs(Base *);//sushil
double FindReTheta(Zone *zone,ZoneBoundary *zBc,int indices[3], double *Theta, double *delta_star, double *delta);
double FindUavg(Zone *zone,ZoneBoundary *zBc,int indices[3]);

#endif

