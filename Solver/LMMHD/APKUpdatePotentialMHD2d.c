#include "APKUpdatePotentialMHD2d.h"
/**
The function UpdatePotentialMHD2d() does the point by point Gauss Seidel Iterations for the PotentialMHD Poisson equation for 2D domain.
Note the diffusion coefficient for PotentialMHD is a constant and has a value of 1 always.
*/
double UpdatePotentialMHD2d(Base *base,Zone *zone)
{
	double residual=0.0,residZone=0.0,omega=1.0,*potentialMHDDomain,sumRHS,sumDiff,diffdia;
	double *RHS,potentialMHDDiff[6];
	int NI,NJ,NK,i,j,k,iFace,index[3],indexP,*ptrSize;
	double *slabVolume[6],*slabPotentialMHD[6];
	double *beta,*volumeDomain,*diffusionCoeffdomain;
	double delta_t,CoeffN[6];
	
	
	//printf("Inside updatePressure\n");
	delta_t = base->BaseParameters->TimeStep;
	potentialMHDDomain = zone->Flowsol->PotentialMHD;
	volumeDomain = zone->ZGeom->Volume;
	ptrSize = zone->Size ;
	beta=zone->ZGeom->Beta;
	
	NI=*(ptrSize+0);
	NJ=*(ptrSize+1);
	NK=*(ptrSize+2);
	
	MakeSlabPtrArray(zone,slabPotentialMHD,17);/*Called for making slabs of PotentialMHD*/
	MakeSlabPtrArray(zone,slabVolume,5);/*Called for making slabs of Volumes*/
	diffusionCoeffdomain = MakeDiffusionCoeffDomain(zone,17);
	//flagrho = zone->Flowsol->flagRho;
	
	for(k=0;k<NK-1;k++)
	for(j=0;j<NJ-1;j++)
	for(i=0;i<NI-1;i++)
	{
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		
		index[0]=k,index[1]=j,index[2]=i;
		
		//rho = *(zone->Flowsol->Rho+indexP*flagrho);
		
		diffdia = (-1) *Diffusion2dWithConstDiffCoeff(index,potentialMHDDomain,volumeDomain,ptrSize,slabPotentialMHD,slabVolume,diffusionCoeffdomain,beta,potentialMHDDiff,CoeffN);
		
		//printf("diffdia = %lf\n",diffdia);
		RHS=zone->Flowsol->PotentialMHDpoissonRHS+indexP*6;/*assign the pointer here*/
		
		sumDiff=0.0;
		sumRHS=0.0;
		for(iFace=0;iFace<4;iFace++)
		{
			sumDiff=sumDiff+potentialMHDDiff[iFace];
			sumRHS=sumRHS +(*(RHS+iFace)); 
		}
		residual= -sumRHS + sumDiff;
		//fprintf(test,"%d\t%d\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\n",i,j,sumRHS,*(RHS+0),*(RHS+1),*(RHS+2),*(RHS+3));
		//printf("i=%d j=%d k=%d sumMassFlux = %12.8elf sumDiff = %12.8elf\n",i,j,k,sumMassFlux,sumDiff);
		*(potentialMHDDomain+indexP) = *(potentialMHDDomain+indexP) + omega*residual/diffdia;
		
		residZone=residZone+residual*residual;
	}
	residZone=sqrt(residZone/((NI-1)*(NJ-1)*(NK-1)));
	//fclose(test);
	//exit(0);
	//printf("residZone = %lf\n",residZone);
	return(residZone);
}

