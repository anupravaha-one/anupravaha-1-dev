#include "APKElectricFluxXYZ.h"
/**
This function ElectricFluxXYZ() is used for calculating the electric flux
(-1.0*Grad(Phi)) which is used for convective terms in UpdateVarCharge().
The electric flux also includes the ion mobility (b) which is merely a constant
value.
*/
int ElectricFluxXYZ(Base *base)
{
	int nZones,numZones;
	int k,j,i,face, index[3];
	double *domainVolume,*Potential;
	double *slabVolume[6],*slabPot[6];	
	double *SurfArea,*F, *beta;
	double CoeffN[6], Fdiff[6];
	Zone* zone;
	nZones=(base->Number_of_zones);
	//FILE *fp1;
	//fp1 = fopen("APKElectricFlux.txt","w");
	int NI,NJ,NK,*ptrSize;
	
	
	for(numZones=0;numZones<nZones;numZones++)
	{
		zone=(base->CgnsZone+numZones);
		
		ptrSize=zone->Size;

		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);

		F = zone->Flowsol->Electricflux;

		SurfArea=zone->ZGeom->SurfaceArea;

		Potential=MakeDomainArray(zone,8);
		domainVolume=MakeDomainArray(zone,5);

		MakeSlabPtrArray(zone,slabPot,8);
		MakeSlabPtrArray(zone,slabVolume,5);
		beta= zone->ZGeom->Beta;

		for(k=0;k<NK-1;k++)
		for(j=0;j<NJ-1;j++)
		for(i=0;i<NI-1;i++)
		{
			

			index[0]=k,index[1]=j,index[2]=i;
			
			//printf("\nIon Mobility = %e\n", *(zone->Flowsol->MobIon));

			Diffusion2dWithConstDiffCoeff(index,Potential,domainVolume,ptrSize,slabPot,slabVolume,(zone->Flowsol->MobIon),beta,Fdiff,CoeffN);
			

			
			face = 0;  // west face
				*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = -1.0*Fdiff[0];
				//fprintf(fp1,"=%.8lf\t",(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face))/(*(zone->Flowsol->MobIon+0)));
			face = 1;  // east face
				*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = -1.0*Fdiff[1];
				//fprintf(fp1,"=%.8lf\t",(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face))/(*(zone->Flowsol->MobIon+0)));
			face = 2;  // south face
				*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = -1.0*Fdiff[2];
				//fprintf(fp1,"=%.8lf\t",(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face))/(*(zone->Flowsol->MobIon+0)));
			face = 3;  //north face
				*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = -1.0*Fdiff[3];
				//fprintf(fp1,"=%.8lf\t",(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face))/(*(zone->Flowsol->MobIon+0)));
			face = 4;  // bottom face
				*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = -1.0*Fdiff[4];
				//fprintf(fp1,"=%.8lf\t",(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face))/(*(zone->Flowsol->MobIon+0)));
			face = 5;  // top face
				*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = -1.0*Fdiff[5]; 
				//fprintf(fp1,"=%.8lf\n",(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face))/(*(zone->Flowsol->MobIon+0)));
		}/*End of looping over NI,NJ,NK*/
	}/*End of looping over the zones */
	//fclose(fp1);
	return 0;
}//End of Main
