#include "APKPressure_Grad_Calculation.h"

/**
This function is used for calculation of pressure gradient of each cell at once and storing the value in an array.  This value of Pressure gradient will be used for subsequent calculation of Navier-Stokes equations.
*/


int Pressure_Grad_Calculation(Base *base)
{
	int NI,NJ,NK,index[3],*ptrSize,indexP;
	int k,j,i,numZones,nZones,flagrho;
	double *variable[3],rho,grad[3];
	
	nZones = base->Number_of_zones;
	Zone *zone;
	
	for(numZones=0;numZones<nZones;numZones++)
	{
		zone = (base->CgnsZone+numZones);
		ptrSize= zone->Size;
		
		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);
		
		variable[0] = zone->Flowsol->Pressure_Grad[0];
		variable[1] = zone->Flowsol->Pressure_Grad[1];
		
		if(base->BaseParameters->flag2d==0)
		variable[2] = zone->Flowsol->Pressure_Grad[2];
		indexP=0;
		flagrho = zone->Flowsol->flagRho;
		if(base->BaseParameters->flag2d==1)
		{
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{	
				if(base->flagMPFV==0)
				{
					rho = *(zone->Flowsol->Rho+indexP*flagrho);
				}
				else if(base->flagDroplet==1 && base->BaseParameters->flagVariableDia==1)
				{
					rho = 1.0;
				}
				else
				{
					rho = 1.0;
				}

				index[0]=k;index[1]=j;index[2]=i;
				CalcPressureGrad_Common(base,zone->Flowsol->P,zone->Flowsol->SlabPress, zone->ZGeom->Volume,zone->Flowsol->SlabVol,zone->ZGeom->SurfaceArea,index,ptrSize,grad);
				*(variable[0]+indexP)=grad[0]/rho;
				
				*(variable[1]+indexP)=grad[1]/rho;
				indexP++;
				
			}
		}
		else if(base->BaseParameters->flag2d==0)
		{
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				if(base->flagMPFV==0)
				{
					rho = *(zone->Flowsol->Rho+indexP*flagrho);
				}
				else if(base->flagDroplet==1 && base->BaseParameters->flagVariableDia==1)
				{
					rho = 1.0;
				}
				else
				{
					rho = 1.0;
				}

				index[0]=k;index[1]=j;index[2]=i;
				
				CalcPressureGrad_Common(base,zone->Flowsol->P,zone->Flowsol->SlabPress, zone->ZGeom->Volume,zone->Flowsol->SlabVol,zone->ZGeom->SurfaceArea,index,ptrSize,grad);
				*(variable[0]+indexP)=grad[0]/rho;
				*(variable[1]+indexP)=grad[1]/rho;
				*(variable[2]+indexP)=grad[2]/rho;
				indexP++;
			}
		}
	}
	return(0);
}


