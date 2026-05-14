#ifndef __APKRefreshMassVelocities_h
#define __APKRefreshMassVelocities_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"

#include "APKExtractMassVelocities.h"
// #include"APKvelocitymassBc.h"

int RefreshMassVelocities(Base *base);
int RefreshMassVelocity(int vartag,Base *base);

#endif

