#include "APKMfindZone2Variable.h"
/**
This function, findZone2Variable() is for transfering the values of a variable from the neighbouring zone sharing face cells into the
interface zone slabs. These values for a variable are required when the solver is solving at the interface.
*/

int MfindZone2Variable(int variableTag,Interface *facePtr,double *P,Zone *zone,int index2[3],int indexofSlabs)
{
	int k,j,i,NI,NJ,NK;
	int *sizePtrNeighbor;
	
	sizePtrNeighbor=zone->Size;
	i=index2[0]-1;
	j=index2[1]-1;
	k=index2[2]-1;
	
	NI=*(sizePtrNeighbor+0);
	NJ=*(sizePtrNeighbor+1);
	NK=*(sizePtrNeighbor+2);
	
//	*(slabPointer+indexofSlabs)= *(domainPointer+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
	
/*	if(i == 0)
	{
		*P = *((zone->grid)->DummyV[0]+(NJ-1)*k+j);
	}
	else if(i == NI-2)
	{
		*P = *((zone->grid)->DummyV[1]+(NJ-1)*k+j);
	}
	else if(j == 0)
	{
		*P = *((zone->grid)->DummyV[2]+(NI-1)*k+i);
	}
	else if(j == NJ-2)
	{
		*P = *((zone->grid)->DummyV[3]+(NI-1)*k+i);
	}
	else if(k == 0)
	{
		*P = *((zone->grid)->DummyV[4]+(NI-1)*j+i);
	}
	else if(k == NK-2)
	{
		*P = *((zone->grid)->DummyV[5]+(NI-1)*j+i);
	}*/
	
	return 0;
}
