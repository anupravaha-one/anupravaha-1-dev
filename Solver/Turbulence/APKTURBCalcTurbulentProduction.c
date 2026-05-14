#include "APKTURBCalcTurbulentProduction.h"

/**
CalcTurbulentProduction function calculates the production term which comes in the K and Epsilon equations. (See Report for further details
regarding Production term). This function calls TurbProduction() function which loops over the cells and according to
2D or 3D problem, calls grad2D() or grad3D() which actually calculates the production term.
*/
int CalcTurbulentProduction(Base *base)
{
	int nZones;
	Zone *zone;
	int numZones;
	
	nZones=(base->Number_of_zones);
	
	base->BaseParameters->Turb->TProductionMaxLimitCount = 0;
	
	for(numZones=0;numZones<nZones;numZones++)
	{
		zone = (base->CgnsZone+numZones);
		
		if(zone->flagSolidFluid==0)
		{
			TurbProduction(base,zone,zone->Flowsol->Tcoeff->TProduction);
		}
	}
	
	if(base->BaseParameters->Turb->TProductionMaxLimitCount != 0)//This should be done only for k-epsilon model
		printf("######### Max. Turbulent Production is limited in %d Cell(s) ########\n", base->BaseParameters->Turb->TProductionMaxLimitCount);
	
	return 0;
}
/**
TurbProduction() function loops over all the cells in a zone and according to
2D or 3D problem, calls grad2D() or grad3D() which actually calculates the turbulent production term.
*/
int TurbProduction(Base *base,Zone *zone,double *Gk)
{
	int i,j,k,NI,NJ,NK,cellIndex,flagrho;
	double *velocity_U,*velocity_V,*velocity_W;
	double *slabVelocity_U[6],*slabVelocity_V[6],*slabVelocity_W[6],*slabVolume[6];
	double *ptrVolume,volume,tViscosity,rho;
	int *ptrSize,index[3];
	
	ptrSize= zone->Size;
	ptrVolume= zone->ZGeom->Volume;
	flagrho = zone->Flowsol->flagRho;
	NI=*(ptrSize+0);
	NJ=*(ptrSize+1);
	NK=*(ptrSize+2);
	
	//Calculating Production term based on Current time velocity field
	if(base->flagMPFV==0)
	{
		velocity_U=MakeDomainArray(zone,1);
		velocity_V=MakeDomainArray(zone,2);
		velocity_W=MakeDomainArray(zone,3);
	
		MakeSlabPtrArray(zone,slabVelocity_U,1);
		MakeSlabPtrArray(zone,slabVelocity_V,2);
		MakeSlabPtrArray(zone,slabVelocity_W,3); 
	}
	else if(base->flagMPFV==1)
	{
		velocity_U=MakeDomainArray(zone,31);
		velocity_V=MakeDomainArray(zone,32);
		velocity_W=MakeDomainArray(zone,33);
	
		MakeSlabPtrArray(zone,slabVelocity_U,31);
		MakeSlabPtrArray(zone,slabVelocity_V,32);
		MakeSlabPtrArray(zone,slabVelocity_W,33); 
	}
	
// 	printf("\nDomain velocity--%lf",*(velocity_U+0));
// 	printf("\nThe slab velocity 0th--%lf",*(slabVelocity_U[0]+0));
// 	printf("\nThe slab velocity 1st--%lf",*(slabVelocity_U[1]+0));
// 	printf("\nThe slab velocity 2nd--%lf",*(slabVelocity_U[2]+0));
// 	printf("\nThe slab velocity 3rd--%lf",*(slabVelocity_U[3]+0));
// 	printf("\nThe slab velocity 4th--%lf",*(slabVelocity_U[4]+0));
// 	printf("\nThe slab velocity 5th--%lf",*(slabVelocity_U[5]+0)); 
	
	//Calculating Production term based on Past velocity field
	/*velocity_U=MakeDomainArray(zone,10);
	velocity_V=MakeDomainArray(zone,20);
	velocity_W=MakeDomainArray(zone,30);

	MakeSlabPtrArray(zone,slabVelocity_U,10);
	MakeSlabPtrArray(zone,slabVelocity_V,20);
	MakeSlabPtrArray(zone,slabVelocity_W,30); */
	
	MakeSlabPtrArray(zone,slabVolume,5);
	
	if(base->BaseParameters->flag2d == 1)
	{
		for(k=0;k<NK-1;k++)
		for(j=0;j<NJ-1;j++)
		for(i=0;i<NI-1;i++)
		{
			index[0]=k,index[1]=j,index[2]=i;
			cellIndex = (NJ-1)*(NI-1)*k+(NI-1)*j+i;

			tViscosity = *(zone->Flowsol->Tcoeff->TViscoeff+cellIndex);
			rho = *(zone->PrevFlowsol->Rho+cellIndex*flagrho);
			volume = *(ptrVolume + cellIndex);
			
			
			
			*(Gk+cellIndex) = (volume/rho) * grad2D(base,zone,velocity_U,slabVelocity_U,velocity_V,slabVelocity_V,index); //(volume*tViscosity/rho)
			//As Turbulent viscosity changes inside the dual time step loop tViscosity not included here in computation of Turbulent Production Term
// 			if(base->flagTurbModel == 1)
// 			{
// 				Eps = *(zone->Flowsol->Scal[1]+cellIndex);
// 				TProductionMaxLimit = fabs(20.0 * volume * Eps / tViscosity);//Limit Production outside the dual time loop after calculating the turbulence visocisty
// 				if(*(Gk+cellIndex) > TProductionMaxLimit)
// 				{
// 					*(Gk+cellIndex) = TProductionMaxLimit;
// 					base->BaseParameters->Turb->TProductionMaxLimitCount++;
// 				}
// 			}
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

			tViscosity = *(zone->Flowsol->Tcoeff->TViscoeff+cellIndex);
			rho = *(zone->PrevFlowsol->Rho+cellIndex*flagrho);
			volume = *(ptrVolume + cellIndex);
			
			
			
			*(Gk+cellIndex) = (volume/rho) * grad3D(base,zone,velocity_U,slabVelocity_U,velocity_V,slabVelocity_V,velocity_W,slabVelocity_W,index);
			
// 			if(base->flagTurbModel == 1)
// 			{
// 				Eps = *(zone->Flowsol->Scal[1]+cellIndex);
// 				TProductionMaxLimit = fabs(20.0 * volume * Eps / tViscosity);//Limit Production outside the dual time loop after calculating the turbulence visocisty
// 				if(*(Gk+cellIndex) > TProductionMaxLimit)
// 				{
// 					*(Gk+cellIndex) = TProductionMaxLimit;
// 					base->BaseParameters->Turb->TProductionMaxLimitCount++;
// 				}
// 			}
		}
	}
	
return(0);
}


