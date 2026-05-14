#include "APKArrayDuplicator.h"

/** 
This fucction expands an array defined on cell centers to one of the size cell vertices by duplication of data.


The inputs to this function are number of vertices of bnoundary patch in either of the two directions out of 
i,j,k according to the direction of BC patch. Pointer to the variable array which is defined on all the face 
centers of the boundary cells. A local array is initialized then the input face centered size array is copied 
to the vertex size array duplicating the data required to fill the extra array points. The address of the local 
array is returned to the parent function. 
*/
double *ArrayDuplicator(int size1,int size2,double *solution)
{
	
	double *localBC;
	int j,k;
	localBC=(double*)malloc(size1*size2*sizeof(double));
	
	for(k=0;k<(size2);k++)
	{
		for(j=0;j<(size1);j++)
		{
			if(j==(size1-1) && k<(size2-1))	
			{
			*(localBC+(size1*k+j))=*(solution+(size1-1)*k+j-1);
			}
			else if(j<(size1-1) && k==(size2-1))
			{
			*(localBC+(size1)*k+j)=*(solution+(size1-1)*(k-1)+j);
			}
			else if(j==((size1)-1) && k==((size2)-1))
			{
			*(localBC+(size1)*k+j)=*(solution+((size1)-1)*(k-1)+j-1);
			}
			else 
			{
			*(localBC+(size1)*k+j)=*(solution+((size1)-1)*(k)+j);
			}	
		}
	}
	return(localBC);
}
