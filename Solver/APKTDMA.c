#include "APKTDMA.h"
/**
The function TDMA() is used to implement the tridiagonal matrix algorithm.
*/
int TDMA(int nI,double *B,double *AP,double *AN,double *X)
{
	int i;
	double r;
	
	for(i=1;i<nI;i++)
	{
		r= *(AN+i*6+2)/(*(AP+i-1));
		*(AP+i)=*(AP+i)-r*(*(AN+(i-1)*6+3));
		*(B+i)=*(B+i)-r*(*(B+i-1));
	}
	*(X+nI-1)=*(B+nI-1)/(*(AP+nI-1));
	for(i=nI-2;i>=0;i--)
	{
		*(X+i)=(*(B+i)-(*(AN+i*6+3))*(*(X+i+1)))/(*(AP+i));
		//printf("X=%lf\n",*(X+i));
	}
	
	return (0);
}
