#include "APKMfindZone2Variable2d.h"

int MfindZone2Variable2d(int variableTag,Interface *facePtr,double *P,Zone *zone,int index2[3],int indexofSlabs)
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
	
	*P = *((zone->grid)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
	
	return 0;
}
