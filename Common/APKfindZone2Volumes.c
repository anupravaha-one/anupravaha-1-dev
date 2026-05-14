#include "APKfindZone2Volumes.h"
/**
This function, findZone2Volumes() is used to get the neighbour zone cell volume. This is required in the solver when we are at the interface.
Zone pointer of the neighbour zone and the cell indices, for which volume is calculated are given as inputs.
*/

double findZone2Volumes(Zone *zone,int index2[3])
{
	int k,j,i,NI,NJ,NK;
	int *sizePtrNeighbor;
	/*printf("i m in findZone2Variable\n");*/
	sizePtrNeighbor=zone->Size;
	i=index2[0]-1;
	j=index2[1]-1;
	k=index2[2]-1;
	
	NI=*(sizePtrNeighbor+0);
	NJ=*(sizePtrNeighbor+1);
	NK=*(sizePtrNeighbor+2);
	return (*((zone->ZGeom)->Volume+(NI-1)*(NJ-1)*k+(NI-1)*j+i));
}
