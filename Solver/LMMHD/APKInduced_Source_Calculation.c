#include "APKInduced_Source_Calculation.h"

/**
This function is used for calculation of source term of each cell at once and storing the value in an array.
This value of source will be used for subsequent calculation of Induction equations.//raj
*/


int Induced_Source_Calculation(Base *base)
{
	int NI,NJ,NK,index[3],*ptrSize,indexP;
	int k,j,i,numZones,nZones;
	double *variable[3],grad[3];
	
	nZones = base->Number_of_zones;
	Zone *zone;
	
	for(numZones=0;numZones<nZones;numZones++)
	{
		zone = (base->CgnsZone+numZones);
		ptrSize= zone->Size;
		
		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);
		
		variable[0] = zone->Flowsol->LMMHD_Source[0];
		variable[1] = zone->Flowsol->LMMHD_Source[1];
		
		if(base->BaseParameters->flag2d==0)
			variable[2] = zone->Flowsol->LMMHD_Source[2];
		indexP=0;
		//flagrho = zone->Flowsol->flagRho;
		if(base->BaseParameters->flag2d==1)
		{
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				//rho = *(zone->Flowsol->Rho+indexP*flagrho);
				index[0]=k;index[1]=j;index[2]=i;
				
				InducedMagneticSource2d(zone,NI,NJ,NK,index,grad);
				
				*(variable[0]+indexP)=grad[0];
				
				
				*(variable[1]+indexP)=grad[1];
				
				indexP++;
				
			}
		}
		else if(base->BaseParameters->flag2d==0)
		{
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				//rho = *(zone->Flowsol->Rho+indexP*flagrho);
				index[0]=k;index[1]=j;index[2]=i;
				
				InducedMagneticSource(zone,NI,NJ,NK,index,grad);
				
				*(variable[0]+indexP)=grad[0];
				
				*(variable[1]+indexP)=grad[1];

				*(variable[2]+indexP)=grad[2];
				indexP++;
			}
		}
	}
	
	return(0);
}


