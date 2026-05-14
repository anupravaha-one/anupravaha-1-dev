#include "APKConvfluxExt2D.h"
/**
This function computes the convective flux for QUICK scheme for 2 dimensional geometries for those cells
which have more than one boundary or interface among its four faces.
This function returns both the Upwind flux(calculated using current iteration values) and the Deferred correction flux
(which is calculated using the previous iteration values). See Eqn(2.21) in Project Report No 1 Dec 2005, for further details.
*/ 
double ConvfluxExt2D(Zone *zone,double *F,double *phicellNew,double *phicell,double *volume,double *FDn,double *FUn,int index[3],double NCoeff[6],double *slabvolume[6],double *slabPhi[6],int face,double Alpha)
{ 
 	int NI,NJ,NK,*ptrSize;
	int i,j,k;
	double volN,volC,phiU,phiN,phiCU,phiCe,phiC;//,PhiSou,volUU,phiUU;
 	double dummy;
	
	int INDEX,I_Increment,J_Increment,K_Increment;
	
	ptrSize=zone->Size;
	NI=*(ptrSize+0);
	NJ=*(ptrSize+1);
	NK=*(ptrSize+2);
	
	I_Increment=1;
	J_Increment=NI-1;
	K_Increment=(NI-1)*(NJ-1);
	
	k=index[0];
	j=index[1];
	i=index[2];
	
	INDEX = (NJ-1)*(NI-1)*k+(NI-1)*j+i;
 	
 	switch(face)
	{
		case 0:
		{
			volC=*(volume+INDEX);
			phiC=*(phicell+INDEX);

			if(*(F+INDEX*6+face)<=0.0)
			{
				if(i!=0)
				{
					phiU = *(phicell+INDEX-1);
					phiCU=*(phicellNew+INDEX-1);
				}
				else
				{
					phiU=*(slabPhi[1]+(NJ-1)*k+j);
					phiCU=*(slabPhi[1]+(NJ-1)*k+j);
				}
				dummy=0.0;				
				NCoeff[0]=0.0;
			}
			else
			{
				phiU = *(phicell+INDEX);
				phiCU= *(phicellNew+INDEX);
				dummy=*(F+INDEX*6+face);
				NCoeff[0]=0.0;
			}
			
			if(i!=0)
			{
				volN=*(volume +INDEX-1);
				phiN=*(phicell+INDEX-1);
			}
			else
			{
				volN=*(slabvolume[1]+(NJ-1)*k+j);
				phiN=*(slabPhi[1]+(NJ-1)*k+j);
			}
			
			phiCe=(phiN*volC+phiC*volN)/(volN+volC);

			*(FDn+INDEX*6+face)=*(F+INDEX*6+face)*(phiCe-phiU)*Alpha;
			*(FUn+INDEX*6+face)=*(F+INDEX*6+face)*(phiCU);
		}
		break;
		
		case 1:
		{
			volC=*(volume+INDEX);
			phiC=*(phicell+INDEX);

			if(*(F+INDEX*6+face)<=0.0)
			{
				if(i!=NI-2)
				{
					phiU = *(phicell+INDEX+1);
					phiCU=*(phicellNew+INDEX+1);
				}
				else
				{
					phiU=*(slabPhi[0]+(NJ-1)*k+j);
					phiCU=*(slabPhi[0]+(NJ-1)*k+j);
				}
				dummy=0.0;				
				NCoeff[1]=0.0;
			}
			else
			{
				phiU = *(phicell+INDEX);
				phiCU= *(phicellNew+INDEX);
				dummy=*(F+INDEX*6+face);
				NCoeff[1]=0.0;
			}
			
			if(i!=NI-2)
			{
				volN=*(volume +INDEX+1);
				phiN=*(phicell+INDEX+1);
			}
			else
			{
				volN=*(slabvolume[0]+(NJ-1)*k+j);
				phiN=*(slabPhi[0]+(NJ-1)*k+j);
			}
			
			phiCe=(phiN*volC+phiC*volN)/(volN+volC);

			*(FDn+INDEX*6+face)=*(F+INDEX*6+face)*(phiCe-phiU)*Alpha;
			*(FUn+INDEX*6+face)=*(F+INDEX*6+face)*(phiCU);
		}
		break;

		case 2:
		{
			volC=*(volume+INDEX);
			phiC=*(phicell+INDEX);

			if(*(F+INDEX*6+face)<=0.0)
			{
				if(j!=0)
				{
					phiU = *(phicell+INDEX-J_Increment);
					phiCU=*(phicellNew+INDEX-J_Increment);
				}
				else
				{
					phiU=*(slabPhi[3]+(NI-1)*k+i);
					phiCU=*(slabPhi[3]+(NI-1)*k+i);
				}
				dummy=0.0;				
				NCoeff[2]=0.0;
			}
			else
			{
				phiU = *(phicell+INDEX);
				phiCU= *(phicellNew+INDEX);
				dummy=*(F+INDEX*6+face);
				NCoeff[2]=0.0;
			}
			
			if(j!=0)
			{
				volN=*(volume +INDEX-J_Increment);
				phiN=*(phicell+INDEX-J_Increment);
			}
			else
			{
				volN=*(slabvolume[3]+(NI-1)*k+i);
				phiN=*(slabPhi[3]+(NI-1)*k+i);
			}
			
			phiCe=(phiN*volC+phiC*volN)/(volN+volC);

			*(FDn+INDEX*6+face)=*(F+INDEX*6+face)*(phiCe-phiU)*Alpha;
			*(FUn+INDEX*6+face)=*(F+INDEX*6+face)*(phiCU);
		}
		break;
 
		case 3:
		{
			volC=*(volume+INDEX);
			phiC=*(phicell+INDEX);

			if(*(F+INDEX*6+face)<=0.0)
			{
				if(j!=NJ-2)
				{
					phiU = *(phicell+INDEX+J_Increment);
					phiCU=*(phicellNew+INDEX+J_Increment);
				}
				else
				{
					phiU=*(slabPhi[2]+(NI-1)*k+i);
					phiCU=*(slabPhi[2]+(NI-1)*k+i);
				}
				dummy=0.0;				
				NCoeff[3]=0.0;
			}
			else
			{
				phiU = *(phicell+INDEX);
				phiCU= *(phicellNew+INDEX);
				dummy=*(F+INDEX*6+face);
				NCoeff[3]=0.0;
			}
			
			if(j!=NJ-2)
			{
				volN=*(volume +INDEX+J_Increment);
				phiN=*(phicell+INDEX+J_Increment);
			}
			else
			{
				volN=*(slabvolume[2]+(NI-1)*k+i);
				phiN=*(slabPhi[2]+(NI-1)*k+i);
			}
			
			phiCe=(phiN*volC+phiC*volN)/(volN+volC);

			*(FDn+INDEX*6+face)=*(F+INDEX*6+face)*(phiCe-phiU)*Alpha;
			*(FUn+INDEX*6+face)=*(F+INDEX*6+face)*(phiCU);
		}	
		break;
		
		case 4:
				*(FDn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=0.0;
				*(FUn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=0.0;

				dummy=0.0;				
				NCoeff[4]=0.0;
		break;
		
		case 5:
			*(FDn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=0.0;
			*(FUn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=0.0;

			dummy=0.0;				
			NCoeff[5]=0.0;
		break;
 	}

	return(dummy);
}
/**************************************************************************************************/
/*double SOUIntpo(double volu,double voluu,double phiu,double phiuu)
{
	double phiSou,coeffU,coeffUU;
	
	
	coeffU=(2*volu+voluu)/(volu+voluu);
	coeffUU=volu/(volu+voluu);
	
	phiSou=coeffU*phiu-coeffUU*phiuu;

return(phiSou);

}*/

