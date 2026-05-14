#include "APKCalcUcrossB.h"
int CalcUcrossB(Base *base)
{
	int nZones,numZones;
	int k,j,i,indexP,NI,NJ,NK,*ptrSize;
	double *domainU,*domainV,*domainW,*slabU[6],*slabV[6],*slabW[6];
	double *domainB_X,*domainB_Y,*domainB_Z,*slabB_X[6],*slabB_Y[6],*slabB_Z[6];
	
	double *domainUcrossB_X,*domainUcrossB_Y,*domainUcrossB_Z,*slabUcrossB_X[6],*slabUcrossB_Y[6],*slabUcrossB_Z[6];
	Zone* zone;
	nZones=(base->Number_of_zones);
	for(numZones=0;numZones<nZones;numZones++)
	{
		zone=(base->CgnsZone+numZones);
		
		ptrSize=zone->Size;

		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);
		
		domainU=MakeDomainArray(zone,1);
		domainV=MakeDomainArray(zone,2);
		domainW=MakeDomainArray(zone,3);
		MakeSlabPtrArray(zone,slabU,1);
		MakeSlabPtrArray(zone,slabV,2);
		MakeSlabPtrArray(zone,slabW,3);

		domainUcrossB_X=MakeDomainArray(zone,22);
		domainUcrossB_Y=MakeDomainArray(zone,23);
		domainUcrossB_Z=MakeDomainArray(zone,24);
		MakeSlabPtrArray_com(zone,slabUcrossB_X,22);
		MakeSlabPtrArray_com(zone,slabUcrossB_Y,23);
		MakeSlabPtrArray_com(zone,slabUcrossB_Z,24);

             	{
			domainB_X=MakeDomainArray(zone,25);
			domainB_Y=MakeDomainArray(zone,26);
			domainB_Z=MakeDomainArray(zone,27);
		
		
			MakeSlabPtrArray_com(zone,slabB_X,25);
			MakeSlabPtrArray_com(zone,slabB_Y,26);
			MakeSlabPtrArray_com(zone,slabB_Z,27);
		}
		
					
		indexP=0;
		for(k=0;k<NK-1;k++)
		for(j=0;j<NJ-1;j++)
		for(i=0;i<NI-1;i++)
		{
			//indexP=(NJ-1)*(NI-1)*k+(NI-1)*j+i;
			
			*(domainUcrossB_X+indexP)= (*(domainV+indexP)* *(domainB_Z+indexP))- (*(domainW+indexP)* *(domainB_Y+indexP));
			*(domainUcrossB_Y+indexP)= (*(domainW+indexP)* *(domainB_X+indexP))- (*(domainU+indexP)* *(domainB_Z+indexP));
			*(domainUcrossB_Z+indexP)= (*(domainU+indexP)* *(domainB_Y+indexP))- (*(domainV+indexP)* *(domainB_X+indexP));
			
			if(i==(NI-2))
			{
				*(slabUcrossB_X[0]+(NJ-1)*k+j) = (*(slabV[0]+(NJ-1)*k+j)* *(slabB_Z[0]+(NJ-1)*k+j))-(*(slabW[0]+(NJ-1)*k+j)* *(slabB_Y[0]+(NJ-1)*k+j));
				*(slabUcrossB_Y[0]+(NJ-1)*k+j) = (*(slabW[0]+(NJ-1)*k+j)* *(slabB_X[0]+(NJ-1)*k+j))-(*(slabU[0]+(NJ-1)*k+j)* *(slabB_Z[0]+(NJ-1)*k+j));
				*(slabUcrossB_Z[0]+(NJ-1)*k+j) = (*(slabU[0]+(NJ-1)*k+j)* *(slabB_Y[0]+(NJ-1)*k+j))-(*(slabV[0]+(NJ-1)*k+j)* *(slabB_X[0]+(NJ-1)*k+j));
			}
			if(i==0)
			{
				*(slabUcrossB_X[1]+(NJ-1)*k+j) = (*(slabV[1]+(NJ-1)*k+j)* *(slabB_Z[1]+(NJ-1)*k+j))-(*(slabW[1]+(NJ-1)*k+j)* *(slabB_Y[1]+(NJ-1)*k+j));
				*(slabUcrossB_Y[1]+(NJ-1)*k+j) = (*(slabW[1]+(NJ-1)*k+j)* *(slabB_X[1]+(NJ-1)*k+j))-(*(slabU[1]+(NJ-1)*k+j)* *(slabB_Z[1]+(NJ-1)*k+j));
				*(slabUcrossB_Z[1]+(NJ-1)*k+j) = (*(slabU[1]+(NJ-1)*k+j)* *(slabB_Y[1]+(NJ-1)*k+j))-(*(slabV[1]+(NJ-1)*k+j)* *(slabB_X[1]+(NJ-1)*k+j));
			}
			if(j==(NJ-2))
			{
				*(slabUcrossB_X[2]+(NI-1)*k+i) = (*(slabV[2]+(NI-1)*k+i)* *(slabB_Z[2]+(NI-1)*k+i))-(*(slabW[2]+(NI-1)*k+i)* *(slabB_Y[2]+(NI-1)*k+i));
				*(slabUcrossB_Y[2]+(NI-1)*k+i) = (*(slabW[2]+(NI-1)*k+i)* *(slabB_X[2]+(NI-1)*k+i))-(*(slabU[2]+(NI-1)*k+i)* *(slabB_Z[2]+(NI-1)*k+i));
				*(slabUcrossB_Z[2]+(NI-1)*k+i) = (*(slabU[2]+(NI-1)*k+i)* *(slabB_Y[2]+(NI-1)*k+i))-(*(slabV[2]+(NI-1)*k+i)* *(slabB_X[2]+(NI-1)*k+i));
			}
			if(j==0)
			{
				*(slabUcrossB_X[3]+(NI-1)*k+i) = (*(slabV[3]+(NI-1)*k+i)* *(slabB_Z[3]+(NI-1)*k+i))-(*(slabW[3]+(NI-1)*k+i)* *(slabB_Y[3]+(NI-1)*k+i));
				*(slabUcrossB_Y[3]+(NI-1)*k+i) = (*(slabW[3]+(NI-1)*k+i)* *(slabB_X[3]+(NI-1)*k+i))-(*(slabU[3]+(NI-1)*k+i)* *(slabB_Z[3]+(NI-1)*k+i));
				*(slabUcrossB_Z[3]+(NI-1)*k+i) = (*(slabU[3]+(NI-1)*k+i)* *(slabB_Y[3]+(NI-1)*k+i))-(*(slabV[3]+(NI-1)*k+i)* *(slabB_X[3]+(NI-1)*k+i));
			}
			if(k==(NK-2))
			{
				*(slabUcrossB_X[4]+(NI-1)*j+i) = (*(slabV[4]+(NI-1)*j+i)* *(slabB_Z[4]+(NI-1)*j+i))-(*(slabW[4]+(NI-1)*j+i)* *(slabB_Y[4]+(NI-1)*j+i));
				*(slabUcrossB_Y[4]+(NI-1)*j+i) = (*(slabW[4]+(NI-1)*j+i)* *(slabB_X[4]+(NI-1)*j+i))-(*(slabU[4]+(NI-1)*j+i)* *(slabB_Z[4]+(NI-1)*j+i));
				*(slabUcrossB_Z[4]+(NI-1)*j+i) = (*(slabU[4]+(NI-1)*j+i)* *(slabB_Y[4]+(NI-1)*j+i))-(*(slabV[4]+(NI-1)*j+i)* *(slabB_X[4]+(NI-1)*j+i));
			}
			if(k==0)
			{
				*(slabUcrossB_X[5]+(NI-1)*j+i) = (*(slabV[5]+(NI-1)*j+i)* *(slabB_Z[5]+(NI-1)*j+i))-(*(slabW[5]+(NI-1)*j+i)* *(slabB_Y[5]+(NI-1)*j+i));
				*(slabUcrossB_Y[5]+(NI-1)*j+i) = (*(slabW[5]+(NI-1)*j+i)* *(slabB_X[5]+(NI-1)*j+i))-(*(slabU[5]+(NI-1)*j+i)* *(slabB_Z[5]+(NI-1)*j+i));
				*(slabUcrossB_Z[5]+(NI-1)*j+i) = (*(slabU[5]+(NI-1)*j+i)* *(slabB_Y[5]+(NI-1)*j+i))-(*(slabV[5]+(NI-1)*j+i)* *(slabB_X[5]+(NI-1)*j+i));
			}
			indexP++;
		}
	}
	return 0;
}



int CalcUcrossB2d(Base *base)
{
	int nZones,numZones;
	int k,j,i,indexP,NI,NJ,NK,*ptrSize;
	double *domainU,*domainV,*slabU[6],*slabV[6];
	double *domainB_X,*domainB_Y,*domainB_Z,*slabB_X[6],*slabB_Y[6],*slabB_Z[6];
	double *domainUcrossB_X,*domainUcrossB_Y,*domainUcrossB_Z,*slabUcrossB_X[6],*slabUcrossB_Y[6],*slabUcrossB_Z[6];
	Zone* zone;

	//	FILE *fp; 
	//	fp =fopen("intermediate.dat","w");

	nZones=(base->Number_of_zones);
	for(numZones=0;numZones<nZones;numZones++)
	{
		zone=(base->CgnsZone+numZones);
		
		ptrSize=zone->Size;

		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);
		
		domainU=MakeDomainArray(zone,1);
		domainV=MakeDomainArray(zone,2);
		//domainW=MakeDomainArray(zone,3);
		MakeSlabPtrArray(zone,slabU,1);
		MakeSlabPtrArray(zone,slabV,2);
		//MakeSlabPtrArray(zone,slabW,3);

		domainUcrossB_X=MakeDomainArray(zone,22);
		domainUcrossB_Y=MakeDomainArray(zone,23);
		domainUcrossB_Z=MakeDomainArray(zone,24);
		MakeSlabPtrArray_com(zone,slabUcrossB_X,22);
		MakeSlabPtrArray_com(zone,slabUcrossB_Y,23);
		MakeSlabPtrArray_com(zone,slabUcrossB_Z,24);
	
		{
			domainB_X=MakeDomainArray(zone,25);
			domainB_Y=MakeDomainArray(zone,26);
			domainB_Z=MakeDomainArray(zone,27);
		
		
			MakeSlabPtrArray_com(zone,slabB_X,25);
			MakeSlabPtrArray_com(zone,slabB_Y,26);
			MakeSlabPtrArray_com(zone,slabB_Z,27);
		}
					
		indexP=0;
		for(k=0;k<NK-1;k++)
		for(j=0;j<NJ-1;j++)
		for(i=0;i<NI-1;i++)
		{
			//indexP=(NJ-1)*(NI-1)*k+(NI-1)*j+i;

		  //	  fprintf(fp, "%lf \t %lf \t %lf \t %lf \t %lf \t %lf\n",*(slabUcrossB_Z[0]+(NJ-1)*k+j),*(slabUcrossB_Z[1]+(NJ-1)*k+j),	*(slabUcrossB_Z[2]+(NI-1)*k+i),	*(slabUcrossB_Z[3]+(NI-1)*k+i),	*(slabUcrossB_Z[4]+(NI-1)*j+i),	*(slabUcrossB_Z[5]+(NI-1)*j+i));
			
			*(domainUcrossB_X+indexP)= *(domainV+indexP) * *(domainB_Z+indexP);
			*(domainUcrossB_Y+indexP)= -1.0* *(domainU+indexP) * *(domainB_Z+indexP);
			*(domainUcrossB_Z+indexP)= (*(domainU+indexP) * *(domainB_Y+indexP))-(*(domainV+indexP) * *(domainB_X+indexP));
			
			if(i==(NI-2))
			{
				*(slabUcrossB_X[0]+(NJ-1)*k+j) = *(slabV[0]+(NJ-1)*k+j) * *(slabB_Z[0]+(NJ-1)*k+j);
				*(slabUcrossB_Y[0]+(NJ-1)*k+j) = -1.0* *(slabU[0]+(NJ-1)*k+j) * *(slabB_Z[0]+(NJ-1)*k+j);
				*(slabUcrossB_Z[0]+(NJ-1)*k+j) = ((*(slabU[0]+(NJ-1)*k+j) * *(slabB_Y[0]+(NJ-1)*k+j))-(*(slabV[0]+(NJ-1)*k+j) * *(slabB_X[0]+(NJ-1)*k+j)));
			}
			if(i==0)
			{
				*(slabUcrossB_X[1]+(NJ-1)*k+j) = *(slabV[1]+(NJ-1)*k+j) * *(slabB_Z[1]+(NJ-1)*k+j);
				*(slabUcrossB_Y[1]+(NJ-1)*k+j) = -1.0* *(slabU[1]+(NJ-1)*k+j) * *(slabB_Z[1]+(NJ-1)*k+j);
				*(slabUcrossB_Z[1]+(NJ-1)*k+j) = ((*(slabU[1]+(NJ-1)*k+j) * *(slabB_Y[1]+(NJ-1)*k+j))-(*(slabV[1]+(NJ-1)*k+j) * *(slabB_X[1]+(NJ-1)*k+j)));
			}
			if(j==(NJ-2))
			{
				*(slabUcrossB_X[2]+(NI-1)*k+i) = *(slabV[2]+(NI-1)*k+i) * *(slabB_Z[2]+(NI-1)*k+i);
				*(slabUcrossB_Y[2]+(NI-1)*k+i) = -1.0* *(slabU[2]+(NI-1)*k+i) * *(slabB_Z[2]+(NI-1)*k+i);
				*(slabUcrossB_Z[2]+(NI-1)*k+i) = ((*(slabU[2]+(NI-1)*k+i) * *(slabB_Y[2]+(NI-1)*k+i))-(*(slabV[2]+(NI-1)*k+i) * *(slabB_X[2]+(NI-1)*k+i)));
			}
			if(j==0)
			{
				*(slabUcrossB_X[3]+(NI-1)*k+i) = *(slabV[3]+(NI-1)*k+i) * *(slabB_Z[3]+(NI-1)*k+i);
				*(slabUcrossB_Y[3]+(NI-1)*k+i) = -1.0* *(slabU[3]+(NI-1)*k+i) * *(slabB_Z[3]+(NI-1)*k+i);
				*(slabUcrossB_Z[3]+(NI-1)*k+i) = ((*(slabU[3]+(NI-1)*k+i) * *(slabB_Y[3]+(NI-1)*k+i))-(*(slabV[3]+(NI-1)*k+i) * *(slabB_X[3]+(NI-1)*k+i)));
			}
			if(k==(NK-2))
			{
				*(slabUcrossB_X[4]+(NI-1)*j+i) = *(slabV[4]+(NI-1)*j+i) * *(slabB_Z[4]+(NI-1)*j+i);
				*(slabUcrossB_Y[4]+(NI-1)*j+i) = -1.0* *(slabU[4]+(NI-1)*j+i) * *(slabB_Z[4]+(NI-1)*j+i);
				*(slabUcrossB_Z[4]+(NI-1)*j+i) = ((*(slabU[4]+(NI-1)*j+i) * *(slabB_Y[4]+(NI-1)*j+i))-(*(slabV[4]+(NI-1)*j+i) * *(slabB_X[4]+(NI-1)*j+i)));
			}
			if(k==0)
			{
				*(slabUcrossB_X[5]+(NI-1)*j+i) = *(slabV[5]+(NI-1)*j+i) * *(slabB_Z[5]+(NI-1)*j+i);
				*(slabUcrossB_Y[5]+(NI-1)*j+i) = -1.0* *(slabU[5]+(NI-1)*j+i) * *(slabB_Z[5]+(NI-1)*j+i);
				*(slabUcrossB_Z[5]+(NI-1)*j+i) = ((*(slabU[5]+(NI-1)*j+i) * *(slabB_Y[5]+(NI-1)*j+i))-(*(slabV[5]+(NI-1)*j+i) * *(slabB_X[5]+(NI-1)*j+i)));
			}




			indexP++;
		}
	}
	//		fclose(fp);
	return 0;
}
