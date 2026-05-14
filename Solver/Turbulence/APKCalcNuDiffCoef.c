#include "APKCalcNuDiffCoef.h"
/**
This function calculates nu_eff = 1/sigma*(nu + nu_t); where 'sigma' is a constant in the Spalart-Allmaras turbulence model equal to a value of '2.0/3.0' 
*/
int CalcNuDiffCoef(Base *base)
{
	int NI,NJ,NK,numZones,cellindex,iface,zBc,index;
	int i,j,k,*ptrSize,flagrho,flagViscoeff,dir;
	int totalCells;//,NIcell,NJcell,NKcell;
	double constant,*domainPtr,*slabIptr,*slabBcptr,dummyViscosity=0.0,mu,nu,nu_t,rho;
	Interface *IfacePtr;
	ZoneBoundary *zBcPtr;
	Zone *zone;
	FILE *fp;

	fp=fopen("Settings/DummyViscocity.dat","r");
	fscanf(fp,"%lf", &dummyViscosity);
	fclose(fp);
	
	for(numZones=0; numZones < base->Number_of_zones; numZones++)
	{
		zone = base->CgnsZone+numZones;

		if(zone->flagSolidFluid==0)
		{
			ptrSize = zone->Size;
			NI=*(ptrSize+0);
			NJ=*(ptrSize+1);
			NK=*(ptrSize+2);

			constant = 1.5;//1.0/sigma;
			
			domainPtr = zone->Flowsol->Scalcoeff[0];
			
			flagrho = zone->Flowsol->flagRho;
			flagViscoeff = zone->Flowsol->flagViscoeff;
			rho = *(zone->PrevFlowsol->Rho);
			mu = *(zone->Flowsol->Viscoeff);
			
			cellindex = 0;
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				if(flagrho == 1)
				{
					rho = *(zone->PrevFlowsol->Rho+cellindex);
				}
				if(flagViscoeff == 1)
				{
					mu   = *(zone->Flowsol->Viscoeff+cellindex);
				}
				nu = mu/rho;
				nu_t = *(zone->PrevFlowsol->Scal[0]+cellindex);
				
				*(domainPtr+cellindex) = constant*(nu + nu_t) + dummyViscosity/rho;
				//printf("%d %d %d Mu_effective Domain = %e nu=%e nu_t=%e dummy=%e\n",k,j,i,*(domainPtr+cellindex),nu,nu_t,dummyViscosity);
				cellindex++;
			}
			
			for(iface=0;iface < zone->Number_of_Interfaces;iface++)
			{
				IfacePtr = zone->Iface+iface;
				dir = IfacePtr->Direction;
				
				slabIptr=IfacePtr->ZSlab->Scalcoeff[0];
				
				switch(dir)
				{
					case 1://Imax face
					case 2://Imin face
						totalCells = (NJ-1) * (NK-1);
						break;
					case 3://Jmax face
					case 4://Jmin face
						totalCells = (NI-1) * (NK-1);
						break;
					case 5://Kmax face
					case 6://Kmin face
						totalCells = (NI-1) * (NJ-1);
						break;
				}

				for(index=0; index < totalCells; index++)
				{
					if(flagrho == 1)
					{
						rho = *(IfacePtr->ZSlab->Rho+index);
					}
					if(flagViscoeff == 1)
					{
						mu   = *(IfacePtr->ZSlab->Viscoeff+index);
					}
					
					nu = mu/rho;
					nu_t = *(IfacePtr->ZSlab->ScalOld[0]+index);
						
					*(slabIptr+index) = constant*(nu + nu_t) + dummyViscosity/rho;
				}
			}
			
			for(zBc=0;zBc < zone->Number_of_BC_Patches;zBc++)
			{
				zBcPtr = zone->ZoneBc+zBc;
				dir = zBcPtr->Direction;
				
				slabBcptr=zBcPtr->ZSlab->Scalcoeff[0];
				
				switch(dir)
				{
					case 1://Imax face
					case 2://Imin face
						totalCells = (NJ-1) * (NK-1);
						break;
					case 3://Jmax face
					case 4://Jmin face
						totalCells = (NI-1) * (NK-1);
						break;
					case 5://Kmax face
					case 6://Kmin face
						totalCells = (NI-1) * (NJ-1);
						break;
				}

				for(index=0; index < totalCells; index++)
				{
					if(flagrho == 1)
					{
						rho = *(zBcPtr->ZSlab->Rho+index);
					}
					if(flagViscoeff == 1)
					{
						mu  = *(zBcPtr->ZSlab->Viscoeff+index);
					}
					
					nu = mu/rho;
					nu_t = *(zBcPtr->ZSlab->ScalOld[0]+index);
					
					*(slabBcptr+index) = constant*(nu + nu_t) + dummyViscosity/rho;
					
					//printf("dir=%d k=%d j=%d i=%d vartag=%d Mu_effective slab = %e nu=%e constant=%e nu_t=%e\n",dir,k,j,i,vartag,*(slabBcptr+index),nu,constant,nu_t);
				}
			}
		}
	}
	return 0;
}

