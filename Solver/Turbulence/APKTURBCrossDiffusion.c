#include "APKTURBCrossDiffusion.h"

int CrossDiffusion(Base *base)
{
	Zone *zone;
	int i,j,k,cntZone,numZones,NI,NJ,NK;
	int index[3],indexP,*ptrSize;
	int flagViscoeff,flagrho;
	double *KOld,*OmegaOld,*volume,*SurfArea;
	double *slabKOld[6],*slabOmegaOld[6],*slabVolume[6];
	double DK[3],DOmega[3],kOld,omegaOld,SigmaOmega2,F1=0.0,CrossDiff=0.0;
	double vol,rho,mu,nu,y,CDkw,arg1;

	numZones=(base->Number_of_zones);

	for(cntZone=0;cntZone < numZones;cntZone++)
	{
		zone=(base->CgnsZone+cntZone);
		flagrho = zone->Flowsol->flagRho;
		flagViscoeff = zone->Flowsol->flagViscoeff;

		if(zone->flagSolidFluid == 0)
		{
			ptrSize= zone->Size;

			NI=*(ptrSize+0);
			NJ=*(ptrSize+1);
			NK=*(ptrSize+2);

			zone->Flowsol->CrossDiff = (double *)malloc((NI-1)*(NJ-1)*(NK-1)*sizeof(double));
			if(zone->Flowsol->CrossDiff == NULL)
			{
				printf("Memory Allocation Failed in CrossDiffusion()\n");
				exit(420);
			}


			KOld	= zone->PrevFlowsol->Scal[0];
			OmegaOld= zone->PrevFlowsol->Scal[1];

			MakeSlabPtrArray(zone,slabKOld,6000);
			MakeSlabPtrArray(zone,slabOmegaOld,6010);

			SurfArea=zone->ZGeom->SurfaceArea;

			volume=zone->ZGeom->Volume;
			MakeSlabPtrArray(zone,slabVolume,5);

			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				index[0]=k,index[1]=j,index[2]=i;
				indexP=k*(NI-1)*(NJ-1)+j*(NI-1)+i;

				kOld = *(KOld + indexP);
				omegaOld = *(OmegaOld + indexP);
				SigmaOmega2 = zone->Flowsol->Tcoeff->SigmaOmega2;
				vol = *(zone->ZGeom->Volume+indexP);
				rho = *(zone->PrevFlowsol->Rho+indexP*flagrho);
				mu = *(zone->Flowsol->Viscoeff+indexP*flagViscoeff);
				nu = mu/rho;

				CalcGradVar(base,zone,KOld,slabKOld,index,DK);
				CalcGradVar(base,zone,OmegaOld,slabOmegaOld,index,DOmega);
				CrossDiff = DK[0]*DOmega[0] + DK[1]*DOmega[1] + DK[2]*DOmega[2];

				CDkw = max(2.0*rho*SigmaOmega2/omegaOld*CrossDiff , pow(10.0,-20));
				y= *(zone->ZGeom->minWallDistance+indexP);
				arg1 = min(max((sqrt(kOld)/(0.09*omegaOld*y)),(500.0*nu/(y*y*omegaOld))) , 4.0*rho*SigmaOmega2*kOld/(CDkw*y*y));
				F1 = tanh(arg1*arg1*arg1*arg1);
				*(zone->Flowsol->Tcoeff->F1 + indexP) = F1;

				*(zone->Flowsol->CrossDiff+indexP) = 2.0 *(1.0-F1)*SigmaOmega2*1.0/omegaOld * CrossDiff * vol;
				//Here itself calculate all the other SST parameters
			}
		}
	}

	return(110);
}


