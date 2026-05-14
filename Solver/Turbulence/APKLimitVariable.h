#ifndef __APKLimitVariable_h
#define __APKLimitVariable_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKData_Structures.h"
#include"APK_File_Open.h"

int LimitVariable(Base *base,int vartag);
int LimitVariableZone(Base *base, Zone *zone, int vartag);

#endif

