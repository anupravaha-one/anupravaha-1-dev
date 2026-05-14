#include "APKElectricFieldXYZ.h"

/**
The function ElectricFieldXYZ() computes the three electric field components from
the electric potential. It loops over all the zones and finds the domain pointers
and slab pointers of potential and computes electric field components.
*/

int ElectricFieldXYZ(Base *base)
{
	Zone *zone;
	int i,j,k,cntZone,numZones,NI,NJ,NK;
	int index[3],indexP,*ptrSize;
	
	double Ey,Ez,Ex,*EfieldX,*EfieldY,*EfieldZ;
	double *potentialdomain,*volume,*SurfArea;
	double *potdomainSlab[6],*slabVolume[6];
	double p,q,r,Area; 

	
	numZones=(base->Number_of_zones);
	
	for(cntZone=0;cntZone < numZones ; cntZone++ )
	{

		zone=(base->CgnsZone+cntZone);
		ptrSize= zone->Size;
		
		
		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);

		volume=MakeDomainArray(zone,5);
		MakeSlabPtrArray(zone,slabVolume,5);
				
		potentialdomain=MakeDomainArray(zone,8);
		MakeSlabPtrArray(zone,potdomainSlab,8);
		
		EfieldX = zone->Flowsol->EfieldX;
		EfieldY = zone->Flowsol->EfieldY;
		EfieldZ = zone->Flowsol->EfieldZ;

		SurfArea=zone->ZGeom->SurfaceArea;

		

		for(k=0;k<NK-1;k++)
		for(j=0;j<NJ-1;j++)
		for(i=0;i<NI-1;i++)
		{
			index[0]=k,index[1]=j,index[2]=i;

			indexP=k*(NI-1)*(NJ-1)+j*(NI-1)+i;
			
			p=*(SurfArea+(indexP)*6*3+3*3+0);
			q=*(SurfArea+(indexP)*6*3+3*3+1);
			r=*(SurfArea+(indexP)*6*3+3*3+2);
			Area=sqrt(p*p+q*q+r*r);
			
			Ex = CalcPressureGrad(potentialdomain,potdomainSlab,volume,slabVolume,SurfArea,index,ptrSize,1);
			Ey = CalcPressureGrad(potentialdomain,potdomainSlab,volume,slabVolume,SurfArea,index,ptrSize,2);
			Ez = CalcPressureGrad(potentialdomain,potdomainSlab,volume,slabVolume,SurfArea,index,ptrSize,3);


			*(EfieldX+indexP) = Ex/(*(volume+indexP));
			*(EfieldY+indexP) = Ey/(*(volume+indexP));
			*(EfieldZ+indexP) = Ez/(*(volume+indexP));
			//*(EfieldZ+indexP)=0.0;

		}
		
	}
 
	return 0;
}
