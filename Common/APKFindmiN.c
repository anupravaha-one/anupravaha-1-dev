#include "APKFindmiN.h"
/**
This function, FindmiN() is used to find minimum deltaX, deltaY and deltaZ in the domain. These values are used to find deltaT(time)
by using CFL condition.
*/
int FindmiN(GridCoordinate *zgrid,int *size,double *minVarX, double *minVarY, double *minVarZ)
{
	double *x,*y,*z;
	double comVarX,comVarY,comVarZ;
	int Ni,Nj,Nk,i,j,k;
	
	Ni=*(size+0);
	Nj=*(size+1);
	Nk=*(size+2);
	x=zgrid->X;
	y=zgrid->Y;
	z=zgrid->Z;
	
	
	*minVarX=*(x+1)-*(x+0);
	*minVarY=*(y+1)-*(y+0);
	*minVarZ=*(z+1)-*(z+0);
	
	for(k=0;k<Nk;k++)
	for(j=0;j<Nj;j++)
	for(i=0;i<Ni;i++)
	{
		comVarX=*(x+(Nj-1)*(Ni-1)*k+(Ni-1)*j+i+1)-*(x+(Nj-1)*(Ni-1)*k+(Ni-1)*j+i);
		
		if(*minVarX>comVarX)
		*minVarX=comVarX;
				
		comVarY=*(y+(Nj-1)*(Ni-1)*k+(Ni-1)*(j+1)+i)-*(y+(Nj-1)*(Ni-1)*k+(Ni-1)*j+i);
		
		if(*minVarY>comVarY)
		*minVarY=comVarY;
		
		comVarZ=*(z+(Nj-1)*(Ni-1)*(k+1)+(Ni-1)*j+i)-*(z+(Nj-1)*(Ni-1)*k+(Ni-1)*j+i);
		
		if(*minVarZ>comVarZ)
		*minVarZ=comVarZ;
	}
	
	return(0);
}

