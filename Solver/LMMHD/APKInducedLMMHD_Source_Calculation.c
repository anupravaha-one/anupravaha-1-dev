#include "APKInducedLMMHD_Source_Calculation.h"

/**
This function is used for calculation of pressure gradient of each cell at once and storing the value in an array.
This value of Pressure gradient will be used for subsequent calculation of Navier-Stokes equations.
*/


int InducedLMMHD_Source_Calculation(Base *base)
{
	int NI,NJ,NK,index[3],*ptrSize,indexP;
	int k,j,i,numZones,nZones,flagrho,flagpermeability;
	double *variable[3],rho,grad[3],magperm,Dkhi[3];
	double *volume;
	nZones = base->Number_of_zones;
	Zone *zone;
	
	
	//printf("Calculating LMMHD Momentum Source values in N-S Equations.....\n\n");
	for(numZones=0;numZones<nZones;numZones++)
	{
		zone = (base->CgnsZone+numZones);
		ptrSize= zone->Size;
		volume= zone->ZGeom->Volume;
		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);
		
		variable[0] = zone->Flowsol->LMMHD_Source[0];
		variable[1] = zone->Flowsol->LMMHD_Source[1];
		
		if(base->BaseParameters->flag2d==0)
			variable[2] = zone->Flowsol->LMMHD_Source[2];
		
		indexP=0;
		flagrho = zone->Flowsol->flagRho;
		
		
		flagpermeability=zone->Flowsol->flagpermeability;//rajendra

		
				
		
		if(base->BaseParameters->flag2d==1)
		{
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				rho = *(zone->Flowsol->Rho+indexP*flagrho);
				
				magperm = *(zone->Flowsol->Permeability+indexP*flagpermeability);
				
				index[0]=k;index[1]=j;index[2]=i;

				
				
				CalcKhiGrad_Common(base,zone->Flowsol->TotalMfieldX,zone->Flowsol->Slab_Total_MfieldX,zone->Flowsol->TotalMfieldY,zone->Flowsol->Slab_Total_MfieldY,zone->Flowsol->TotalMfieldZ,zone->Flowsol->Slab_Total_MfieldZ, zone->ZGeom->Volume,zone->Flowsol->SlabVol,zone->ZGeom->SurfaceArea,index,ptrSize,Dkhi);
 				
				
				
				InducedMomentumSource2d(zone,NI,NJ,NK,index,grad);
					


				*(variable[0]+indexP)=(1.0/(2.0*magperm*rho))*(2.0*grad[0] - Dkhi[0]);

				*(variable[1]+indexP)=(1.0/(2.0*magperm*rho))*(2.0*grad[1] - Dkhi[1]);
				

				
				
 				
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
				
				magperm = *(zone->Flowsol->Permeability+indexP*flagpermeability);
				index[0]=k;index[1]=j;index[2]=i;

				CalcKhiGrad_Common(base,zone->Flowsol->TotalMfieldX,zone->Flowsol->Slab_Total_MfieldX,zone->Flowsol->TotalMfieldY,zone->Flowsol->Slab_Total_MfieldY,zone->Flowsol->TotalMfieldZ,zone->Flowsol->Slab_Total_MfieldZ, zone->ZGeom->Volume,zone->Flowsol->SlabVol,zone->ZGeom->SurfaceArea,index,ptrSize,Dkhi);
				
 				
				
				
				
				InducedMomentumSource(zone,NI,NJ,NK,index,grad);
					


				*(variable[0]+indexP)=(1.0/(2.0*magperm*rho))*(2.0*grad[0] - Dkhi[0]);

				*(variable[1]+indexP)=(1.0/(2.0*magperm*rho))*(2.0*grad[1] - Dkhi[1]);
				
				*(variable[2]+indexP)=(1.0/(2.0*magperm*rho))*(2.0*grad[2] - Dkhi[2]);
				
 				
				indexP++;
				
			}
		}
	}
	return(0);
}


