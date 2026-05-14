 
#include "APKCalcFaceValue.h"
/** This function is for finding out the values at face centroids of cells.  **/

double CalcRhoFaceValue(Base *base,Zone *zone,double *variable)
{
 	int i,j,k,NI,NJ,NK,*ptrSize;
	int K_Increment,J_Increment,I_Increment;
	int iFace;
	int indexE,indexW,indexS,indexN,indexP;
	double vare,varw,varn,vars;
	double var_P,var_E,var_W,var_N,var_S;
	double vol_P,vol_E,vol_W,vol_N,vol_S;

	double *variable_face,*slabVolume[6],*domainVolume;
	
	ptrSize= zone->Size;
	
	MakeSlabPtrArray(zone,slabVolume,5);
	domainVolume=MakeDomainArray(zone,5);
	
	variable_face = zone->PrevFlowsol->RhoG_face;
	
	NI=*(ptrSize+0);
	NJ=*(ptrSize+1);
	NK=*(ptrSize+2);
	
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
		 
	  	var_P=*(variable+indexP);
	       	vol_P=*(domainVolume+indexP);
			
		if(i==0)
		{		
			vol_W =*(slabVolume[1]+(NJ-1)*k+j);
			var_W =*(variable+indexP);
		}
		else
		{
			var_W=*(variable+indexW);				
			vol_W=*(domainVolume+indexW);
		}
		
			if(i==(NI-2))
			{
				var_E=*(variable+indexP);					
				vol_E=*(slabVolume[0]+(NJ-1)*k+j);
			}
			else
			{
				var_E=*(variable+indexE);
				vol_E=*(domainVolume+indexE);
			}
			if(j==0)
			{
				var_S=*(variable+indexP);				
				vol_S=*(slabVolume[3]+(NI-1)*k+i);
			}
			else
			{
				var_S=*(variable+indexS);			
				vol_S=*(domainVolume+indexS);				
			}
			if(j==(NJ-2))
			{
				var_N= *(variable+indexP);
				vol_N=*(slabVolume[2]+(NI-1)*k+i);
			}
			else
			{
				var_N=*(variable+indexN);
				vol_N=*(domainVolume+indexN);
			}
		
			

			vare=((var_E*vol_P+var_P*vol_E)/(vol_P+vol_E));
			varw=((var_W*vol_P+var_P*vol_W)/(vol_P+vol_W));
			vars=((var_S*vol_P+var_P*vol_S)/(vol_P+vol_S));
			varn=((var_N*vol_P+var_P*vol_N)/(vol_P+vol_N));

			*(variable_face+(indexP*6)+0)=varw;
			*(variable_face+(indexP*6)+1)=vare;
			*(variable_face+(indexP*6)+2)=vars;
			*(variable_face+(indexP*6)+3)=varn;
			*(variable_face+(indexP*6)+4)=0;
			*(variable_face+(indexP*6)+5)=0;		

	}//end of for loop

return (0.0);
	
}

