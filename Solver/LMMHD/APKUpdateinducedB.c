#include "APKUpdateinducedB.h"

int UpdateinducedB(Base *base)
{
	int nZones,numZones;
	int k,j,i;
	int Index,indexB,indexT,indexE,indexW,indexS,indexN,indexP;

	double *domainVolume,*domainPsi;
	double *slabVolume[6],*slabPsi[6];
	double *SurfArea;
	Zone* zone;
	nZones=(base->Number_of_zones);
	
	int NI,NJ,NK,*ptrSize;
	double psiw,psie,psin,psis,psit,psib;
	double psi_P,psi_E,psi_W,psi_S,psi_N,psi_T,psi_B;
	double vol_P,vol_E,vol_W,vol_N,vol_S,vol_T,vol_B;
	
	for(numZones=0;numZones<nZones;numZones++)
	{
		zone=(base->CgnsZone+numZones);
		
		ptrSize=zone->Size;

		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);

		SurfArea=zone->ZGeom->SurfaceArea;

		domainPsi=MakeDomainArray(zone,59);
		domainVolume=MakeDomainArray(zone,5);

		MakeSlabPtrArray(zone,slabPsi,59);
		MakeSlabPtrArray(zone,slabVolume,5);

		for(k=0;k<NK-1;k++)
		for(j=0;j<NJ-1;j++)
		for(i=0;i<NI-1;i++)
		{
			indexP=(NJ-1)*(NI-1)*k+(NI-1)*j+i;

			indexE=indexP+1;//(NJ-1)*(NI-1)*k+(NI-1)*j+(i+1);
			indexW=indexP-1;//(NJ-1)*(NI-1)*k+(NI-1)*j+(i-1);
			indexN=indexP+(NI-1);//(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+i;
			indexS=indexP-(NI-1);//(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i;
			indexT=indexP+(NI-1)*(NJ-1);//(NJ-1)*(NI-1)*(k+1)+(NI-1)*j+i;
			indexB=indexP-(NI-1)*(NJ-1);//(NJ-1)*(NI-1)*(k-1)+(NI-1)*j+i;
// 			indexE=(NJ-1)*(NI-1)*k+(NI-1)*j+(i+1);
// 			indexW=(NJ-1)*(NI-1)*k+(NI-1)*j+(i-1);
// 			indexN=(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+i;
// 			indexS=(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i;
// 			indexT=(NJ-1)*(NI-1)*(k+1)+(NI-1)*j+i;
// 			indexB=(NJ-1)*(NI-1)*(k-1)+(NI-1)*j+i;

			psi_P=*(domainPsi+indexP);
			vol_P=*(domainVolume+indexP);

			if(i==(NI-2))
			{
				psi_E=*(slabPsi[0]+(NJ-1)*k+j);
				vol_E=*(slabVolume[0]+(NJ-1)*k+j);
			}
			else
			{
				psi_E=*(domainPsi+indexE);
				vol_E=*(domainVolume+indexE);
			}

			if(i==0)
			{
				psi_W =*(slabPsi[1]+(NJ-1)*k+j);
				vol_W =*(slabVolume[1]+(NJ-1)*k+j);
			}
			else
			{
				psi_W=*(domainPsi+indexW);
				vol_W=*(domainVolume+indexW);
			}

			if(j==(NJ-2))
			{
				psi_N=*(slabPsi[2]+(NI-1)*k+i);
				vol_N=*(slabVolume[2]+(NI-1)*k+i);
			}
			else
			{
				psi_N=*(domainPsi+indexN);
				vol_N=*(domainVolume+indexN);
			}

			if(j==0)
			{
				psi_S=*(slabPsi[3]+(NI-1)*k+i);
				vol_S=*(slabVolume[3]+(NI-1)*k+i);
			}
			else
			{
				psi_S=*(domainPsi+indexS);
				vol_S=*(domainVolume+indexS);
			}

			if(k==(NK-2))
			{
				psi_T=*(slabPsi[4]+(NI-1)*j+i);
				vol_T=*(slabVolume[4]+(NI-1)*j+i);
			}
			else
			{
				psi_T=*(domainPsi+indexT);
				vol_T=*(domainVolume+indexT);
			}

			if(k==0)
			{
				psi_B=*(slabPsi[5]+(NI-1)*j+i);
				vol_B=*(slabVolume[5]+(NI-1)*j+i);
			}
			else
			{
				psi_B=*(domainPsi+indexB);
				vol_B=*(domainVolume+indexB);
			}


			psie=((psi_E*vol_P+psi_P*vol_E)/(vol_P+vol_E));
			psiw=((psi_W*vol_P+psi_P*vol_W)/(vol_P+vol_W));
			psis=((psi_S*vol_P+psi_P*vol_S)/(vol_P+vol_S));
			psin=((psi_N*vol_P+psi_P*vol_N)/(vol_P+vol_N));
			psit=((psi_T*vol_P+psi_P*vol_T)/(vol_P+vol_T));
			psib=((psi_B*vol_P+psi_P*vol_B)/(vol_P+vol_B));

			Index=((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6*3;
			*((zone->Flowsol)->TotalMfieldX+((NJ-1)*(NI-1)*k+(NI-1)*j+i)) = *((zone->Flowsol)->TotalMfieldX+((NJ-1)*(NI-1)*k+(NI-1)*j+i))- ((psiw*(*(SurfArea+Index+3*0+0))+psie*(*(SurfArea+Index+3*1+0))+psis*(*(SurfArea+Index+3*2+0))+psin*(*(SurfArea+Index+3*3+0))+psib*(*(SurfArea+Index+3*4+0))+psit*(*(SurfArea+Index+3*5+0)))/vol_P);
			
			*((zone->Flowsol)->TotalMfieldY+((NJ-1)*(NI-1)*k+(NI-1)*j+i)) = *((zone->Flowsol)->TotalMfieldY+((NJ-1)*(NI-1)*k+(NI-1)*j+i))- ((psiw*(*(SurfArea+Index+3*0+1))+psie*(*(SurfArea+Index+3*1+1))+psis*(*(SurfArea+Index+3*2+1))+psin*(*(SurfArea+Index+3*3+1))+psib*(*(SurfArea+Index+3*4+1))+psit*(*(SurfArea+Index+3*5+1)))/vol_P);

			*((zone->Flowsol)->TotalMfieldZ+((NJ-1)*(NI-1)*k+(NI-1)*j+i)) = *((zone->Flowsol)->TotalMfieldZ+((NJ-1)*(NI-1)*k+(NI-1)*j+i))- ((psiw*(*(SurfArea+Index+3*0+2))+psie*(*(SurfArea+Index+3*1+2))+psis*(*(SurfArea+Index+3*2+2))+psin*(*(SurfArea+Index+3*3+2))+psib*(*(SurfArea+Index+3*4+2))+psit*(*(SurfArea+Index+3*5+2)))/vol_P);
		}/*End of looping over NI,NJ,NK*/
		ExtractVariables(56,base,base->CgnsZone+numZones);
		ExtractVariables(57,base,base->CgnsZone+numZones);
		ExtractVariables(58,base,base->CgnsZone+numZones);
	}/*End of looping over the zones */
	return 0;
}
