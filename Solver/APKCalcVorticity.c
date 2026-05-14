#include "APKCalcVorticity.h"
/**
The function CalcVorticity() computes the three vorticity components from the velocities.
It loops over all the zones and finds the domain pointers and slab pointers of velocities and computes
vorticity components.
*/
int CalcVorticity(Base *base)
{
	Zone *zone;
	int i,j,k,cntZone,numZones,NI,NJ,NK;
	int index[3],indexP,*ptrSize;
	
	double dU_dY,dU_dZ,dV_dX,dV_dZ,dW_dX,dW_dY,*vorticityX,*vorticityY,*vorticityZ;
	double *velocity_U,*velocity_V,*velocity_W,*volume,*SurfArea,vol;
	double *slabVelocity_U[6],*slabVelocity_V[6],*slabVelocity_W[6],*slabVolume[6];
	
	numZones=(base->Number_of_zones);
	
	if(base->BaseParameters->flag2d==1)
	{
		for(cntZone=0;cntZone < numZones;cntZone++)
		{
			zone=(base->CgnsZone+cntZone);
			if(zone->flagSolidFluid == 0)
			{
				ptrSize= zone->Size;
	
				NI=*(ptrSize+0);
				NJ=*(ptrSize+1);
				NK=*(ptrSize+2);
	
				volume=MakeDomainArray(zone,5);
	
				vorticityX=zone->Flowsol->vorticity_X;
				vorticityY=zone->Flowsol->vorticity_Y;
				vorticityZ=zone->Flowsol->vorticity_Z;
	
				velocity_U=MakeDomainArray(zone,1);
				velocity_V=MakeDomainArray(zone,2);
	
				MakeSlabPtrArray(zone,slabVelocity_U,1);
				MakeSlabPtrArray(zone,slabVelocity_V,2);
	
				SurfArea=zone->ZGeom->SurfaceArea;
	
				MakeSlabPtrArray(zone,slabVolume,5);
	
				for(k=0;k<NK-1;k++)
				for(j=0;j<NJ-1;j++)
				for(i=0;i<NI-1;i++)
				{
					index[0]=k,index[1]=j,index[2]=i;
	
					indexP=k*(NI-1)*(NJ-1)+j*(NI-1)+i;
	
					dU_dY=CalcPressureGrad(velocity_U,slabVelocity_U,volume,slabVolume,SurfArea,index,ptrSize,2);
					dV_dX=CalcPressureGrad(velocity_V,slabVelocity_V,volume,slabVolume,SurfArea,index,ptrSize,1);
					
					vol = *(volume+indexP);
					
					*(vorticityX+indexP)=0.0;
					*(vorticityY+indexP)=0.0;
					*(vorticityZ+indexP)=(-1.0/vol)*(dV_dX-dU_dY);
				}
			}
		}
	}
	else
	{
		for(cntZone=0;cntZone < numZones;cntZone++)
		{
			zone=(base->CgnsZone+cntZone);
			if(zone->flagSolidFluid == 0)
			{
				ptrSize= zone->Size;
	
				NI=*(ptrSize+0);
				NJ=*(ptrSize+1);
				NK=*(ptrSize+2);
	
				volume=MakeDomainArray(zone,5);
	
				vorticityX=zone->Flowsol->vorticity_X;
				vorticityY=zone->Flowsol->vorticity_Y;
				vorticityZ=zone->Flowsol->vorticity_Z;
	
				velocity_U=MakeDomainArray(zone,1);
				velocity_V=MakeDomainArray(zone,2);
				velocity_W=MakeDomainArray(zone,3);
	
				MakeSlabPtrArray(zone,slabVelocity_U,1);
				MakeSlabPtrArray(zone,slabVelocity_V,2);
				MakeSlabPtrArray(zone,slabVelocity_W,3);
	
				SurfArea=zone->ZGeom->SurfaceArea;
	
				MakeSlabPtrArray(zone,slabVolume,5);
	
				for(k=0;k<NK-1;k++)
				for(j=0;j<NJ-1;j++)
				for(i=0;i<NI-1;i++)
				{
					index[0]=k,index[1]=j,index[2]=i;
	
					indexP=k*(NI-1)*(NJ-1)+j*(NI-1)+i;
	
					dU_dY=CalcPressureGrad(velocity_U,slabVelocity_U,volume,slabVolume,SurfArea,index,ptrSize,2);
					dU_dZ=CalcPressureGrad(velocity_U,slabVelocity_U,volume,slabVolume,SurfArea,index,ptrSize,3);
					dV_dX=CalcPressureGrad(velocity_V,slabVelocity_V,volume,slabVolume,SurfArea,index,ptrSize,1);
					dV_dZ=CalcPressureGrad(velocity_V,slabVelocity_V,volume,slabVolume,SurfArea,index,ptrSize,3);
					dW_dX=CalcPressureGrad(velocity_W,slabVelocity_W,volume,slabVolume,SurfArea,index,ptrSize,1);
					dW_dY=CalcPressureGrad(velocity_W,slabVelocity_W,volume,slabVolume,SurfArea,index,ptrSize,2);
	
					vol = *(volume+indexP);
					*(vorticityX+indexP)=(-1.0/vol)*(dW_dY-dV_dZ);
					*(vorticityY+indexP)=(-1.0/vol)*(dU_dZ-dW_dX);
					*(vorticityZ+indexP)=(-1.0/vol)*(dV_dX-dU_dY);
				}
			}
		}
	}
	return 0;
}


/**
The function CalcAbsVorticity() computes the absolute vorticity from the velocities. Even if the vorticity
is not asked to be calculated by the user this function calculates the absolute vorticity and stores it
in the data structure which will be used in SST turbulence model
*/

int CalcAbsVorticity(Base *base)
{
	Zone *zone;
	int i,j,k,cntZone,numZones,NI,NJ,NK;
	int index[3],indexP,*ptrSize;
	
	double dU_dY,dU_dZ,dV_dX,dV_dZ,dW_dX,dW_dY,VorticityX,VorticityY,VorticityZ,*Vorticity;
	double *velocity_U,*velocity_V,*velocity_W,*volume,*SurfArea,vol;
	double *slabVelocity_U[6],*slabVelocity_V[6],*slabVelocity_W[6],*slabVolume[6];
	
	numZones=(base->Number_of_zones);
	if(base->BaseParameters->flag2d==0)
	{
		for(cntZone=0;cntZone < numZones;cntZone++)
		{
			zone=(base->CgnsZone+cntZone);
			if(zone->flagSolidFluid == 0)
			{
				ptrSize= zone->Size;
		
				NI=*(ptrSize+0);
				NJ=*(ptrSize+1);
				NK=*(ptrSize+2);
		
				volume=MakeDomainArray(zone,5);
		
				velocity_U=MakeDomainArray(zone,1);
				velocity_V=MakeDomainArray(zone,2);
				velocity_W=MakeDomainArray(zone,3);
		
				MakeSlabPtrArray(zone,slabVelocity_U,1);
				MakeSlabPtrArray(zone,slabVelocity_V,2);
				MakeSlabPtrArray(zone,slabVelocity_W,3);
		
				SurfArea=zone->ZGeom->SurfaceArea;
		
				MakeSlabPtrArray(zone,slabVolume,5);
		
				Vorticity = zone->Flowsol->Vorticity;
		
				for(k=0;k<NK-1;k++)
				for(j=0;j<NJ-1;j++)
				for(i=0;i<NI-1;i++)
				{
					index[0]=k,index[1]=j,index[2]=i;
		
					indexP=k*(NI-1)*(NJ-1)+j*(NI-1)+i;
		
					dU_dY=CalcPressureGrad(velocity_U,slabVelocity_U,volume,slabVolume,SurfArea,index,ptrSize,2);
					dU_dZ=CalcPressureGrad(velocity_U,slabVelocity_U,volume,slabVolume,SurfArea,index,ptrSize,3);
					dV_dX=CalcPressureGrad(velocity_V,slabVelocity_V,volume,slabVolume,SurfArea,index,ptrSize,1);
					dV_dZ=CalcPressureGrad(velocity_V,slabVelocity_V,volume,slabVolume,SurfArea,index,ptrSize,3);
					dW_dX=CalcPressureGrad(velocity_W,slabVelocity_W,volume,slabVolume,SurfArea,index,ptrSize,1);
					dW_dY=CalcPressureGrad(velocity_W,slabVelocity_W,volume,slabVolume,SurfArea,index,ptrSize,2);
		
					vol = *(volume+indexP);
					VorticityX =(-1.0/vol)*(dW_dY-dV_dZ);
					VorticityY =(-1.0/vol)*(dU_dZ-dW_dX);
					VorticityZ =(-1.0/vol)*(dV_dX-dU_dY);
		
					*(Vorticity+indexP) = sqrt(VorticityX*VorticityX + VorticityY*VorticityY + VorticityZ*VorticityZ);
					//printf("Vorticity: \tX -> %e\tY -> %e\tZ -> %e\tAbs = %e\n",VorticityX,VorticityY,VorticityZ,*(Vorticity+indexP));
				}
			}
		}
	}
	else
	{
		for(cntZone=0;cntZone < numZones;cntZone++)
		{
			zone=(base->CgnsZone+cntZone);
			if(zone->flagSolidFluid == 0)
			{
				ptrSize= zone->Size;
		
				NI=*(ptrSize+0);
				NJ=*(ptrSize+1);
				NK=*(ptrSize+2);
		
				volume=MakeDomainArray(zone,5);
		
				velocity_U=MakeDomainArray(zone,1);
				velocity_V=MakeDomainArray(zone,2);
		
				MakeSlabPtrArray(zone,slabVelocity_U,1);
				MakeSlabPtrArray(zone,slabVelocity_V,2);
		
				SurfArea=zone->ZGeom->SurfaceArea;
		
				MakeSlabPtrArray(zone,slabVolume,5);
		
				Vorticity = zone->Flowsol->Vorticity;
				indexP=0;
				for(k=0;k<NK-1;k++)
				for(j=0;j<NJ-1;j++)
				for(i=0;i<NI-1;i++)
				{
					index[0]=k,index[1]=j,index[2]=i;
		
					dU_dY=CalcPressureGrad(velocity_U,slabVelocity_U,volume,slabVolume,SurfArea,index,ptrSize,2);
					dV_dX=CalcPressureGrad(velocity_V,slabVelocity_V,volume,slabVolume,SurfArea,index,ptrSize,1);					
					vol = *(volume+indexP);
					VorticityZ =(-1.0/vol)*(dV_dX-dU_dY);
		
					*(Vorticity+indexP) = sqrt(VorticityZ*VorticityZ);
					indexP++;
					//printf("Vorticity: \tX -> %e\tY -> %e\tZ -> %e\tAbs = %e\n",VorticityX,VorticityY,VorticityZ,*(Vorticity+indexP));
				}
			}
		}
	}
	return 0;
}


