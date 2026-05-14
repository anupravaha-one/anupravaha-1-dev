#include "APKFindCellCentreArray.h"

int FindCellCentreArray(int *size,double *coordX,double *coordY,double *coordZ,double *centre)
{
	int NI,NJ,NK,i,j,k;   
	NI=*(size+0);
	NJ=*(size+1);
	NK=*(size+2);
	
	for(k=0;k<NK;k++)
	for(j=0;j<NJ;j++)
	for(i=0;i<NI;i++)
	{
	
	*(centre+(k*(NJ-1)*(NI-1)+(NI-1)*j+i)*3+0)=(*(coordX+k*(NJ-1)*(NI-1)+(NI-1)*j+i)+*(coordX+(k+1)*NJ*NI+NI*j+i)+*(coordX+k*NJ*NI+NI*(j+1)+i)+*(coordX+(k+1)*NJ*NI+NI*(j+1)+i)+*(coordX+k*NJ*NI+NI*j+i+1)+*(coordX+(k+1)*NJ*NI+NI*j+i+1)+*(coordX+k*NJ*NI+NI*(j+1)+i+1)+*(coordX+(k+1)*NJ*NI+NI*(j+1)+i+1))/8.0;
	*(centre+(k*(NJ-1)*(NI-1)+(NI-1)*j+i)*3+1)=(*(coordY+k*(NJ-1)*(NI-1)+(NI-1)*j+i)+*(coordY+(k+1)*NJ*NI+NI*j+i)+*(coordY+k*NJ*NI+NI*(j+1)+i)+*(coordY+(k+1)*NJ*NI+NI*(j+1)+i)+*(coordY+k*NJ*NI+NI*j+i+1)+*(coordY+(k+1)*NJ*NI+NI*j+i+1)+*(coordY+k*NJ*NI+NI*(j+1)+i+1)+*(coordY+(k+1)*NJ*NI+NI*(j+1)+i+1))/8.0;
	*(centre+(k*(NJ-1)*(NI-1)+(NI-1)*j+i)*3+2)=(*(coordZ+k*(NJ-1)*(NI-1)+(NI-1)*j+i)+*(coordZ+(k+1)*NJ*NI+NI*j+i)+*(coordZ+k*NJ*NI+NI*(j+1)+i)+*(coordZ+(k+1)*NJ*NI+NI*(j+1)+i)+*(coordZ+k*NJ*NI+NI*j+i+1)+*(coordZ+(k+1)*NJ*NI+NI*j+i+1)+*(coordZ+k*NJ*NI+NI*(j+1)+i+1)+*(coordZ+(k+1)*NJ*NI+NI*(j+1)+i+1))/8.0;
	
	}
	return 0;   
}
