#include "APKESP_Source_Calculation.h"

/**
This function is used for calculation of pressure gradient of each cell at once and storing the value in an array.
This value of Pressure gradient will be used for subsequent calculation of Navier-Stokes equations.
*/


int ESP_Source_Calculation(Base *base)
{
	int NI,NJ,NK,index[3],*ptrSize,indexP;
	int k,j,i,numZones,nZones,flagrho;
	double *variable[3],rho,grad[3],Scaling,*rhoion;
	
	nZones = base->Number_of_zones;
	Zone *zone;
	
	Scaling = base->BaseParameters->ScaleFactor;
	printf("Calculating ESP Momentum Source values in N-S Equations.....\n\n");
	for(numZones=0;numZones<nZones;numZones++)
	{
		zone = (base->CgnsZone+numZones);
		ptrSize= zone->Size;
		
		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);
		
		variable[0] = zone->Flowsol->ESP_Source[0];
		variable[1] = zone->Flowsol->ESP_Source[1];
		
		if(base->BaseParameters->flag2d==0)
			variable[2] = zone->Flowsol->ESP_Source[2];
		
		indexP=0;
		flagrho = zone->Flowsol->flagRho;
		
		rhoion = zone->Flowsol->RhoIon;
		if(base->BaseParameters->flag2d==1)
		{
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				rho = *(zone->Flowsol->Rho+indexP*flagrho);
				index[0]=k;index[1]=j;index[2]=i;
				CalcPressureGrad_Common(base,zone->Flowsol->Potential,zone->Flowsol->Slab_Potential, zone->ZGeom->Volume,zone->Flowsol->SlabVol,zone->ZGeom->SurfaceArea,index,ptrSize,grad);
				
				
				*(variable[0]+indexP)=grad[0]*Scaling*(*(rhoion+indexP))/rho;
				
				*(variable[1]+indexP)=grad[1]*Scaling*(*(rhoion+indexP))/rho;
				indexP++;
				
			}
		}
		else if(base->BaseParameters->flag2d==0)
		{
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				rho = *(zone->Flowsol->Rho+indexP*flagrho);
				index[0]=k;index[1]=j;index[2]=i;
				
				CalcPressureGrad_Common(base,zone->Flowsol->Potential,zone->Flowsol->Slab_Potential, zone->ZGeom->Volume,zone->Flowsol->SlabVol,zone->ZGeom->SurfaceArea,index,ptrSize,grad);
				*(variable[0]+indexP)=grad[0]*Scaling*(*(rhoion+indexP))/rho;
				*(variable[1]+indexP)=grad[1]*Scaling*(*(rhoion+indexP))/rho;
				*(variable[2]+indexP)=grad[2]*Scaling*(*(rhoion+indexP))/rho;
				indexP++;
			}
		}
	}
	return(0);
}


