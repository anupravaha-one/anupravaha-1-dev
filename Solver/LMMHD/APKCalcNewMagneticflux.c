#include"APKCalcNewMagneticflux.h"

int CalcNewMagneticflux2d(Base *base)
{
	int numzones,cntZ,*ptrSize,NI,NJ,NK;
	int i,j,k,index[3],indexP,iface;
	Zone *zone;
	double psiDiff[6],diffdia,*psiDomain,*volumeDomain,*beta;
	double *slabPsi[6],*slabVolume[6],CoeffN[6];
	double delta_t,*diffusionCoeffdomain;
	
	delta_t = base->BaseParameters->TimeStep;
	numzones= base->Number_of_zones;
	
	for(cntZ=0;cntZ<numzones;cntZ++)
	{
		zone = (base->CgnsZone+cntZ);
		if(zone->flagSolidFluid==0)
		{
			ptrSize = zone->Size;
			NI=*(ptrSize+0);
			NJ=*(ptrSize+1);
			NK=*(ptrSize+2);
			beta = zone->ZGeom->Beta;
			psiDomain=MakeDomainArray(zone,59); //variableDomain stores variable pointer
			volumeDomain=MakeDomainArray(zone,5);

			MakeSlabPtrArray(zone,slabPsi,59);/*Called for making slabs of Psi*/
			MakeSlabPtrArray(zone,slabVolume,5);/*Called for making slabs of Volumes*/
			diffusionCoeffdomain = MakeDiffusionCoeffDomain(zone,59);
			indexP=0;
			for(k=0;k<NK-1;k++)   /*Storing current solution in Previous solution*/
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				//indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;

				index[0]=k,index[1]=j,index[2]=i;
				
				diffdia=Diffusion2dWithConstDiffCoeff_opti(index,psiDomain,volumeDomain,ptrSize,slabPsi,slabVolume,diffusionCoeffdomain,beta,psiDiff,CoeffN);
				for(iface=0;iface<4;iface++)
					*(zone->Flowsol->Magneticflux+indexP*6+iface)=*(zone->Flowsol->Magneticflux+indexP*6+iface) - psiDiff[iface];
				indexP++;
			}
		}
	}
	
	return 0;
}
int CalcNewMagneticflux(Base *base)
{
	int numzones,cntZ,*ptrSize,NI,NJ,NK;
	int i,j,k,index[3],indexP,iface;
	Zone *zone;
	double psiDiff[6],diffdia,*psiDomain,*volumeDomain,*beta;
	double *slabPsi[6],*slabVolume[6];
	double delta_t,*diffusionCoeffdomain;
	
	delta_t = base->BaseParameters->TimeStep;
	numzones= base->Number_of_zones;
	
	for(cntZ=0;cntZ<numzones;cntZ++)
	{
		zone = (base->CgnsZone+cntZ);
		if(zone->flagSolidFluid==0)
		{
			ptrSize = zone->Size;
			NI=*(ptrSize+0);
			NJ=*(ptrSize+1);
			NK=*(ptrSize+2);
			beta = zone->ZGeom->Beta;
			psiDomain=MakeDomainArray(zone,59); //variableDomain stores variable pointer
			volumeDomain=MakeDomainArray(zone,5);

			MakeSlabPtrArray(zone,slabPsi,59);/*Called for making slabs of Psi*/
			MakeSlabPtrArray(zone,slabVolume,5);/*Called for making slabs of Volumes*/
			diffusionCoeffdomain = MakeDiffusionCoeffDomain(zone,59);
			indexP=0;
			for(k=0;k<NK-1;k++)   /*Storing current solution in Previous solution*/
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				//indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;

				index[0]=k,index[1]=j,index[2]=i;
				
				diffdia=Diffusion3dWithConstDiffCoeff_opti(index,psiDomain,volumeDomain,ptrSize,slabPsi,slabVolume,diffusionCoeffdomain,beta,psiDiff);
				for(iface=0;iface<6;iface++)
					*(zone->Flowsol->Magneticflux+indexP*6+iface)=*(zone->Flowsol->Magneticflux+indexP*6+iface) - psiDiff[iface];
				indexP++;
			}
		}
	}
	return 0;
}

