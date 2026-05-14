#include "APKCalGradScalerXY.h"
/**
this function will calculate gradient of scaler in x and y and z diredtion  with help of CalcGradVar() and store in GradScalX,GradScalY and GradScalZ
**/
int CalGradScalerXY(Base *base,int variableTag,int number)
{
	double *Scaler,*slabScaler[6];
	double  D[3];
	Zone *zone;
	int i,j,k,NJ,NK,NI,*ptrSize,numZones,iZone,index1,index[3];
        double *GradX,*GradY,*GradZ;
	numZones=base->Number_of_zones;
		
	for(iZone=0;iZone<numZones;iZone++)
	{      	
	        zone=base->CgnsZone+iZone;
	        GradX=zone->Flowsol->GradScalX[number];
	        GradY=zone->Flowsol->GradScalY[number];
	        GradZ=zone->Flowsol->GradScalZ[number];
		Scaler=zone->Flowsol->Scal[number];
		MakeSlabPtrArray(zone,slabScaler,variableTag);
        	
		ptrSize=zone->Size;
		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);
		
		for(k=0;k<NK-1;k++)
		for(i=0;i<NI-1;i++)
		for(j=0;j<NJ-1;j++)
		{
			index1=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			CalcGradVar(base,zone,Scaler,slabScaler,index,D);
			*(GradX+index1)=D[0];
			*(GradY+index1)=D[1];
			*(GradZ+index1)=D[2];	               
		}
	}
	
	for(iZone=0;iZone<numZones;iZone++)
	{  
		zone=base->CgnsZone+iZone;
				
		Extract_Var_Coeff_PorousMedia(60000+number,base,zone);  
		Extract_Var_Coeff_PorousMedia(60100+number,base,zone);  
		Extract_Var_Coeff_PorousMedia(60200+number,base,zone);  
	}
	return 0;
}
