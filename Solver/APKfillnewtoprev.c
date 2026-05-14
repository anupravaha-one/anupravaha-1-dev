#include "APKfillnewtoprev.h"

/**
This function updates the old (previous iteration) values of the variable with the new 
(latest iteration) values.
*/

int fillnewtoprev(double *Address,double *PAddress,int *ptrSize)
{
	int i,j,k,index,NI,NJ,NK;
	NI=*(ptrSize+0); /*Storing current umass in prevflowsolution's umass (Note: umass or vmass or any var.)*/
	NJ=*(ptrSize+1);
	NK=*(ptrSize+2);
	
	for(k=0;k<NK-1;k++)
	for(j=0;j<NJ-1;j++)
	for(i=0;i<NI-1;i++)
	{
		index=(NJ-1)*(NI-1)*k+(NI-1)*j+i;
		*(PAddress+index) = *(Address+index);
	}
	return 0;
}
