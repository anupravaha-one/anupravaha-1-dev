#include "APKTURBCalcKWallTerm.h"

/**
CalcKWallTerm function calculates the Wall term which comes in the K equations for low reynolds number model(launder and sharma). 
(See Report for further details regarding wall term). This function calls TurbProduction() function which loops over the cells and 
according to 2D or 3D problem, calls grad2D() or grad3D() which actually calculates the wall term for k equation.
*/
int CalcKWallTerm(Base *base)
{
	int nZones;
	Zone *zone;
	int numZones;
	
	nZones=(base->Number_of_zones);
	
	for(numZones=0;numZones<nZones;numZones++)
	{
		zone = (base->CgnsZone+numZones);
		
		if(zone->flagSolidFluid==0)
		{
			KWall(base,zone,zone->Flowsol->Tcoeff->PhiK);
		}
	}
	
	return 0;
}

int KWall(Base *base,Zone *zone,double *Gk)
{
	int NI,NJ,NK,i,j,k,cellIndex,flagrho,flagviscosity;
	int *ptrSize,index[3];
	double *domainK,*slabK[6],*ptrVolume,Viscosity;
	double volume,rho;
	
	ptrSize=zone->Size;
	ptrVolume= zone->ZGeom->Volume;
	flagrho = zone->Flowsol->flagRho; 
	flagviscosity = zone->Flowsol->flagViscoeff;
	NI=*(ptrSize+0);
	NJ=*(ptrSize+1);
	NK=*(ptrSize+2);
	
	domainK = MakeDomainArray(zone,600);
	MakeSlabPtrArray(zone,slabK,600);
	if(base->BaseParameters->flag2d == 1)
	{
		for(k=0;k<NK-1;k++)
		for(j=0;j<NJ-1;j++)
		for(i=0;i<NI-1;i++)
		{
			index[0]=k,index[1]=j,index[2]=i;
			cellIndex = (NJ-1)*(NI-1)*k+(NI-1)*j+i;

			Viscosity = *(zone->Flowsol->Viscoeff+flagviscosity*cellIndex);
			rho = *(zone->PrevFlowsol->Rho+cellIndex*flagrho);
			volume = *(ptrVolume + cellIndex);

			*(Gk+cellIndex) = 2.0*(volume*Viscosity/rho) * Kgradient2D(base,zone,domainK,slabK,index);
		}
	}
	else
	{
		for(k=0;k<NK-1;k++)
		for(j=0;j<NJ-1;j++)
		for(i=0;i<NI-1;i++)
		{
			index[0]=k,index[1]=j,index[2]=i;
			cellIndex = (NJ-1)*(NI-1)*k+(NI-1)*j+i;

			Viscosity = *(zone->Flowsol->Viscoeff+flagviscosity*cellIndex);
			rho = *(zone->PrevFlowsol->Rho+cellIndex*flagrho);
			volume = *(ptrVolume + cellIndex);

			*(Gk+cellIndex) = 2.0*(volume*Viscosity/rho) * Kgradient3D(base,zone,domainK,slabK,index);
		}
	}
	
	return(0);
}

