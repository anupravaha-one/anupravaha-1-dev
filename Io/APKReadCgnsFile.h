#ifndef __APKReadCgnsFile_h
#define __APKReadCgnsFile_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "cgnslib.h"
#include "APKCGNSErrorCheck.h"

#include "APKData_Structures.h"

int read_2d_grid(int, Base *,int,int,FILE *);
Root ReadCgnsFile(char* );

#endif

