#ifndef __APKReadInitialData_h
#define __APKReadInitialData_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include "cgnslib.h"
#include "APKCGNSErrorCheck.h"

#include "APKData_Structures.h"

int ReadInitialData(char* CgnsFileName, Root CgnsRoot);

#endif

