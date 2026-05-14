#ifndef __APKinitialize_variable_h
#define __APKinitialize_variable_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "APKData_Structures.h"

double *initialize_variable(int i,Zone *z,double constant);
double *initializeFlux(int i,Zone *z);
double *initializeDiffFlux(int i,Zone *z);

#endif

