#include "APKUpdatePotential2d.h"

double UpdatePotential2d(Base *base,Zone *zone)
{
	double residual=0.0,residZone=0.0,*PotentialDomain,sumMassFlux,sumDiff,diffdia;
	double PotentialDiff[6];
	int NI,NJ,NK,i,j,k,iFace,index[3],indexP,*ptrSize;
	double *slabVolume[6],*slabPotential[6];
	double *beta,*volumeDomain,*diffusionCoeffdomain;
	double CoeffN[6];
	double *Residue,*CorrectionX,*DiagonalC,*NCoeff,DiagSource,source;
	
	//FILE *sav; 

	//sav=fopen("GS_Omega_Pot.dat","r");
	//fscanf(sav,"%lf",&omega);
	//fscanf(sav,"%lf",&delta_t);
	//fclose(sav);
	
	PotentialDomain = zone->Flowsol->Potential;
	volumeDomain = zone->ZGeom->Volume;
	ptrSize = zone->Size ;
	beta=zone->ZGeom->Beta;
	
	NI=*(ptrSize+0);
	NJ=*(ptrSize+1);
	NK=*(ptrSize+2);
	
	Residue=(double *)malloc((NJ-1)*sizeof(double));
	CorrectionX=(double *)malloc((NJ-1)*sizeof(double));
	DiagonalC=(double *)malloc((NJ-1)*sizeof(double));
	NCoeff=(double *)malloc(6*(NJ-1)*sizeof(double));
	
	MakeSlabPtrArray(zone,slabPotential,8);/*Called for making slabs of Pressure*/
	MakeSlabPtrArray(zone,slabVolume,5);/*Called for making slabs of Volumes*/
	diffusionCoeffdomain = MakeDiffusionCoeffDomain(zone,8);
	
	for(k=0;k<NK-1;k++)
	for(i=0;i<NI-1;i++)
	{
		for(j=0;j<NJ-1;j++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			diffdia = (-1.0)*Diffusion2dWithConstDiffCoeff(index,PotentialDomain,volumeDomain,ptrSize,slabPotential,slabVolume,diffusionCoeffdomain,beta,PotentialDiff,CoeffN);
			
			DiagSource=CalcSource(&source,8,zone,base,index);//Source Term
			source = -1.0*source;
			sumDiff=0.0;
			sumMassFlux=0.0;
			for(iFace=0;iFace<4;iFace++)
			{
				sumDiff=sumDiff+PotentialDiff[iFace];
				*(NCoeff+j*6+iFace) = CoeffN[iFace];
			}
			residual= sumDiff-source ;
			
			*(Residue+j)=residual;
			
			*(DiagonalC+j)=diffdia;
		}
		TDMA(NJ-1,Residue,DiagonalC,NCoeff,CorrectionX);
		for(j=0;j<NJ-1;j++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			*(PotentialDomain+indexP) = *(PotentialDomain+indexP) + *(CorrectionX+j);
			residZone=residZone+(*(CorrectionX+j))*(*(CorrectionX+j));
		}
	}
	residZone=sqrt(residZone/((NI-1)*(NJ-1)*(NK-1)));
	
	//printf("\nResidPotential = %e\n",residZone);
	
	free(Residue);
	free(CorrectionX);
	free(DiagonalC);
	free(NCoeff);
	return(residZone);
}

