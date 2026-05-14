#include "APKTransformCellIndicesCommon.h"

/**
This function transforms the the input cell indices into the index notation system of the neighboring Zone. 


The transformation of the index vector I1 (in present Zone's notation system) into index vector I2 (in the 
neighboring Zone's notation system) follows the rule I2=[T]*I1+A. Where T is the transformation (rotation) matrix 
and vector A is the translation vector.
*/
int TransformCellIndicesCommon(int index1[3],Interface *iFace,int index2[3])
{
	int *addConst,*trans,temp,i,j;
	addConst=iFace->addConstant;
	trans=iFace->Transformation_Matrix;
	
	for(j=0;j<3;j++)
	{
		temp=0;
		for(i=0;i<3;i++)
		temp=temp+(*(trans+j*3+i))*index1[i];
		
		index2[j]=temp+*(addConst+j);
		/*printf("%4d -> %4d\n",index1[j],index2[j]);*/
	}
	return 0;
}
