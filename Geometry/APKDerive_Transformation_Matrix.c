#include "APKDerive_Transformation_Matrix.h"

/** 
This function takes an integer and returns 1 if the integer is positive or zero and -1 if the integer is 
negative.
*/
int sgn(int x)
{

	int return_val = 1;
	if(x<0)  return_val = -1;
	
	return return_val;
}


/**
This function takes two integers and returns 1 if their magnitude is same otherwise returns zero.
*/

int del(int x, int y)
{
	int return_val = 0;
	if(abs(x) == abs(y)) return_val = 1;
	return return_val;
}

/**
This function derives the transformation matrix from the transformation rule stored in the Data structure. 
For details please refer to CGNS Recommended Practice AIAA R-101-2002 Page.No.66.
*/
int Derive_Transformation_Matrix(Interface *Intface) 
{
	/*printf("i'm in deriveTrans func \n");*/
	int i,j;
	Intface->Transformation_Matrix = (int *)malloc(9*sizeof(int));
	
	for(j=0;j<3;j++)
	{
		for(i=0;i<3;i++)
		{
			*(Intface->Transformation_Matrix+3*j+i) = sgn(*(Intface->Transform+i));
			*(Intface->Transformation_Matrix+3*j+i)*= del(*(Intface->Transform+i),j+1);
			/*printf("%4d",*(Intface->Transformation_Matrix+3*j+i));*/
			/*fprintf(k,"%4d",*(Intface->Transformation_Matrix+3*j+i));*/
		}
		/*printf("\n");*/
		/*fprintf(k,"\n");*/
	}
	return 0;
}




