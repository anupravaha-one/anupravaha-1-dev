#include "APKCalcFaceValue.h"
/**
These functions are used to compute the face value for the variables. Added by Mrunalini
*/
int CalcFaceValue2d(Base *base,int variableTag)
{
	int i,j,k,NI,NJ,NK,*ptrSize;
	Zone *zone;
	int cntZ,numZones;
	int indexB,indexT,indexE,indexW,indexS,indexN,indexP;

	double *domainVolume,*domainVariable;
	double *slabVolume[6],*slabVariable[6];	
		
	int K_Increment,J_Increment,I_Increment;

	double vare,varw,varn,vars;
	double var_P,var_E,var_W,var_N,var_S;
	double vol_P,vol_E,vol_W,vol_N,vol_S;
	
	double *F;
	int iFace;

	FILE *fp;
	fp = fopen("SOLUTIONS/rhof.dat","a");
	numZones = base->Number_of_zones;
	
	for(cntZ=0;cntZ<numZones;cntZ++)
	{
		zone = (base->CgnsZone+cntZ);	
		ptrSize=zone->Size;

		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);
		if(variableTag==710)
		{
			F = zone->PrevFlowsol->RhoG_face;
		}
		if(variableTag==650)
			F = zone->PrevFlowsol->VarthetaD_face;
		
		domainVariable= MakeDomainArray(zone,variableTag);
		domainVolume=MakeDomainArray(zone,5);

		MakeSlabPtrArray(zone,slabVariable,variableTag);
		MakeSlabPtrArray(zone,slabVolume,5);

		K_Increment = (NI-1)*(NJ-1);
		J_Increment = (NI-1);
		I_Increment = 1;

		indexP=0;
		
		for(k=0;k<NK-1;k++)
		for(j=0;j<NJ-1;j++)
		for(i=0;i<NI-1;i++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;

			indexE=indexP+1;
			indexW=indexP-1;
			indexN=indexP+J_Increment;
			indexS=indexP-J_Increment;
			indexT=indexP+K_Increment;
			indexB=indexP-K_Increment;

			var_P=*(domainVariable+indexP);
			
			vol_P=*(domainVolume+indexP);
		
			if(i==(NI-2))
			{
				var_E=*(slabVariable[0]+(NJ-1)*k+j);					
				vol_E=*(slabVolume[0]+(NJ-1)*k+j);
			}
			else
			{
				var_E=*(domainVariable+indexE);		
				vol_E=*(domainVolume+indexE);
			}
			
			if(i==0)
			{
				var_W =*(slabVariable[1]+(NJ-1)*k+j);			
				vol_W =*(slabVolume[1]+(NJ-1)*k+j);
			}
			else
			{
				var_W=*(domainVariable+indexW);			
				vol_W=*(domainVolume+indexW);
			}
		
			if(j==(NJ-2))
			{
				var_N=*(slabVariable[2]+(NI-1)*k+i);			
				vol_N=*(slabVolume[2]+(NI-1)*k+i);
			}
			else
			{
				var_N=*(domainVariable+indexN);		
				vol_N=*(domainVolume+indexN);
			}
			if(j==0)
			{
				var_S=*(slabVariable[3]+(NI-1)*k+i);			
				vol_S=*(slabVolume[3]+(NI-1)*k+i);
			}
			else
			{
				var_S=*(domainVariable+indexS);			
				vol_S=*(domainVolume+indexS);				
			}
		
			vare=((var_E*vol_P+var_P*vol_E)/(vol_P+vol_E));
			varw=((var_W*vol_P+var_P*vol_W)/(vol_P+vol_W));
			vars=((var_S*vol_P+var_P*vol_S)/(vol_P+vol_S));
			varn=((var_N*vol_P+var_P*vol_N)/(vol_P+vol_N));
			
			*(F+(indexP*6)+0)   = varw;
			*(F+(indexP*6)+1) = vare;	
			*(F+(indexP*6)+2)  = vars;			
			*(F+(indexP*6)+3)  = varn;		
			*(F+(indexP*6)+4) = 0.0;			
			*(F+(indexP*6)+5)  =0.0;	

			if(variableTag==710)	
			{
				for(iFace =0; iFace<6;iFace++ )	
				{
					fprintf(fp,"\n%d,%d --%d \t %lf",i,j,iFace, *(F+(indexP*6)+iFace));
				}
			}
		
		}//End of looping over NI,NJ,NK
		fprintf(fp,"\n\n");
	}
	fclose(fp);
	return 0;
}

/**
These functions are used to compute the face value for the variables.
*/
int CalcFaceValue(Base *base,int variableTag)
{
	int i,j,k,NI,NJ,NK,*ptrSize;
	Zone *zone;
	int cntZ,numZones;
	int indexB,indexT,indexE,indexW,indexS,indexN,indexP;

	double *domainVolume,*domainVariable;
	double *slabVolume[6],*slabVariable[6];	
		
	int K_Increment,J_Increment,I_Increment;

	double vare,varw,varn,vars,varb,vart;
	double var_P,var_E,var_W,var_N,var_S,var_B,var_T;
	double vol_P,vol_E,vol_W,vol_N,vol_S,vol_B,vol_T;
	
	double *F;
	int face_index;

	numZones = base->Number_of_zones;

	for(cntZ=0;cntZ<numZones;cntZ++)
	{
		zone = (base->CgnsZone+cntZ);	
		ptrSize=zone->Size;

		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);
		if(variableTag==7100)
		{
			F = zone->PrevFlowsol->RhoG_face;
		}
		if(variableTag==650)
			F = zone->PrevFlowsol->VarthetaD_face;
		
		domainVariable= MakeDomainArray(zone,variableTag);
		domainVolume=MakeDomainArray(zone,5);

		MakeSlabPtrArray(zone,slabVariable,variableTag);
		MakeSlabPtrArray(zone,slabVolume,5);

		K_Increment = (NI-1)*(NJ-1);
		J_Increment = (NI-1);
		I_Increment = 1;

		indexP=0;
	
		for(k=0;k<NK-1;k++)
		for(j=0;j<NJ-1;j++)
		for(i=0;i<NI-1;i++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;

			indexE=indexP+1;
			indexW=indexP-1;
			indexN=indexP+J_Increment;
			indexS=indexP-J_Increment;
			indexT=indexP+K_Increment;
			indexB=indexP-K_Increment;

			var_P=*(domainVariable+indexP);
			
			vol_P=*(domainVolume+indexP);
		
			if(i==(NI-2))
			{
				var_E=*(slabVariable[0]+(NJ-1)*k+j);					
				vol_E=*(slabVolume[0]+(NJ-1)*k+j);
			}
			else
			{
				var_E=*(domainVariable+indexE);		
				vol_E=*(domainVolume+indexE);
			}
			
			if(i==0)
			{
				var_W =*(slabVariable[1]+(NJ-1)*k+j);			
				vol_W =*(slabVolume[1]+(NJ-1)*k+j);
			}
			else
			{
				var_W=*(domainVariable+indexW);			
				vol_W=*(domainVolume+indexW);
			}
		
			if(j==(NJ-2))
			{
				var_N=*(slabVariable[2]+(NI-1)*k+i);			
				vol_N=*(slabVolume[2]+(NI-1)*k+i);
			}
			else
			{
				var_N=*(domainVariable+indexN);		
				vol_N=*(domainVolume+indexN);
			}
			if(j==0)
			{
				var_S=*(slabVariable[3]+(NI-1)*k+i);			
				vol_S=*(slabVolume[3]+(NI-1)*k+i);
			}
			else
			{
				var_S=*(domainVariable+indexS);			
				vol_S=*(domainVolume+indexS);				
			}
			if(k==(NK-2))
			{
				var_T=*(slabVariable[4]+(NI-1)*j+i);			
				vol_T=*(slabVolume[4]+(NI-1)*j+i);
			}
			else
			{
				var_T=*(domainVariable+indexT);
				vol_T=*(domainVolume+indexT);
			}
			if(k==0)
			{
				var_B=*(slabVariable[5]+(NI-1)*j+i);
				vol_B=*(slabVolume[5]+(NI-1)*j+i);
			}
			else
			{
				var_B=*(domainVariable+indexB);
				vol_B=*(domainVolume+indexB);			
			}

			vare=((var_E*vol_P+var_P*vol_E)/(vol_P+vol_E));
			varw=((var_W*vol_P+var_P*vol_W)/(vol_P+vol_W));
			vars=((var_S*vol_P+var_P*vol_S)/(vol_P+vol_S));
			varn=((var_N*vol_P+var_P*vol_N)/(vol_P+vol_N));
			vart=((var_T*vol_P+var_P*vol_T)/(vol_P+vol_T));
			varb=((var_B*vol_P+var_P*vol_B)/(vol_P+vol_B));

			*(F+(indexP*6)+0)   = varw;
			*(F+(indexP*6)+1) = vare;	
			*(F+(indexP*6)+2)  = vars;			
			*(F+(indexP*6)+3)  = varn;		
			*(F+(indexP*6)+4) = varb;			
			*(F+(indexP*6)+5)  =vart;			

		}//End of looping over NI,NJ,NK
	}
	return 0;
}

