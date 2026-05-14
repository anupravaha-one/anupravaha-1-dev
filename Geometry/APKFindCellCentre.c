#include "APKFindCellCentre.h"

/**
This function is used to find the coordinates of the center of a cell.
 
As input arguments the vertex Size of the Zone and the addresses (pointer) of the coordinates of the vertex
with minimum indices in a cell are supplied.
Given the vertex size of the Zone and the pointer of the coordinates of the vertex with minimum 
indices, the coordinates of the other vertices can be retrieved. The coordinates of the cell center are 
calculated by taking the average of the respective coordinates of all the eight vertices. 
*/
int FindCellCentre(int *size,double *coordX,double *coordY,double *coordZ,double *centre)
{
	int NI,NJ,NK;   
	NI=*(size+0);
	NJ=*(size+1);
	NK=*(size+2);

	*(centre+0)=(*(coordX)+*(coordX+NI*NJ)+*(coordX+NI)+*(coordX+NI+NJ*NI)+*(coordX+1)+*(coordX+1+NI*NJ)+*(coordX+1+NI)+*(coordX+1+NI+NJ*NI))/8.0;
	*(centre+1)=(*(coordY)+*(coordY+NI*NJ)+*(coordY+NI)+*(coordY+NI+NJ*NI)+*(coordY+1)+*(coordY+1+NI*NJ)+*(coordY+1+NI)+*(coordY+1+NI+NJ*NI))/8.0;
	*(centre+2)=(*(coordZ)+*(coordZ+NI*NJ)+*(coordZ+NI)+*(coordZ+NI+NJ*NI)+*(coordZ+1)+*(coordZ+1+NI*NJ)+*(coordZ+1+NI)+*(coordZ+1+NI+NJ*NI))/8.0;
   
	return 0;   
}
