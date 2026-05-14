#include "APKCorrection.h"

/**
This function performs coarse grid correction operation for 2D as well as 3D domain.
The prolongated error values are added to the solution variable.
Uh = Uh + Eh
*/

double Correction(Base *base,Zone *zone,int g)
{
	int i,j,k,NI,NJ,NK,*ptrSize;
	double *v,*e;
	
	if(g == 0)
	{
		v = zone->Flowsol->P;
		ptrSize = zone->Size ;
	}
	else
	{
		v = (zone->grid+g)->V;
		ptrSize = (zone->grid+g)->GSize;
	}
	NI=*(ptrSize+0);
	NJ=*(ptrSize+1);
	NK=*(ptrSize+2);
	
	e = (zone->grid+g)->E;
	
	for(k=0;k<NK-1;k++)
	for(j=0;j<NJ-1;j++)
	for(i=0;i<NI-1;i++)
	{
		*(v+(NI-1)*(NJ-1)*k+(NI-1)*j+i)=*(v+(NI-1)*(NJ-1)*k+(NI-1)*j+i)+*(e+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
	}
	return 0;
}

