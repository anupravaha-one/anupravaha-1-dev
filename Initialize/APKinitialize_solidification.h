#ifndef __APKinitialize_solidification_h
#define __APKinitialize_solidification_h

#include<malloc.h>
#include<stdio.h>
#include "APKData_Structures.h"

#include"APKinitialize_variable.h"
#include"APKinitializeslab.h"
#include"APKinitializeVarCoefficient.h"
#include"APKinitializeSlabVarCoeff.h"
//#include "../Solidification/bisSolidification_updatevariable.h"
//#include "../Solidification/bisprint_enthalpy.h"
//int Solidification_updatevariable(Base *);
int initialize_solidification(Base *,int *);


/*#include<malloc.h>
#include<stdio.h>
#include "../Common/bisData_Structures.h"
int initialize_solidification(Base *,int *);
double* initializeVarCoefficient(int,Zone *,double,int);
int Solidification_updatevariable(Base *);
double *initializeslab(int,int,int *);
double *initialize_variable(int ,Zone *,double);
double *initializeSlabVarCoeff(int ,Zone *,double constant,int ,int ,int *);*/
#endif

