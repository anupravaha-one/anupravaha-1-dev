#include "APKAdjoiningZonePointerXYZ.h"

/**
This function takes the Zone pointer and index(i,j,k) of a cell and returns the adresses of the x-coordinate, 
y-coordinate and z-coordinate of the vertex of the cell which corresponds to the minimum indices. 
*/
int AdjoiningZonePointerXYZ(Zone *zone,int index2[3],double *pointxyz[3])
{
	
	int NI,NJ,NK,i,j,k;
	
	NI=*(zone->Size+0);
	NJ=*(zone->Size+1);
	NK=*(zone->Size+2);
	//printf("size  of the neighbour zone is %4d%4d%4d\n",NI,NJ,NK);
	
	i=index2[0]-1;
	j=index2[1]-1;
	k=index2[2]-1;
	
	pointxyz[0]=((zone->CgnsGrid)->X+NI*NJ*k+NI*j+i);
	pointxyz[1]=((zone->CgnsGrid)->Y+NI*NJ*k+NI*j+i);
	pointxyz[2]=((zone->CgnsGrid)->Z+NI*NJ*k+NI*j+i);
	//printf("adress1=%u  address2=%u  adress3=    %u\n",pointxyz[0],pointxyz[1],pointxyz[2]);
	return 0;
}
