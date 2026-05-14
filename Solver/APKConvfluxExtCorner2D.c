#include "APKConvfluxExtCorner2D.h"
/**
This function computes the convective flux for central difference scheme for 2 dimensional geometries for those cells
which have more than one boundary or interface among its four faces.
This function returns both the Upwind flux(calculated using current iteration values) and the Deferred correction flux
(which is calculated using the previous iteration values). See Eqn(2.21) in Project Report No 1 Dec 2005, for further details.
*/ 
double ConvfluxExtCorner2D(Zone *zone,double *F,double *phicellNew,double *phicell,double *volume,double *FDn,double *FUn,int index[3],double NCoeff[6],double *slabvolume[6],double *slabPhi[6],int face,double Alpha)
{ 
 	int NI,NJ,NK,*ptrSize;
	int i,j,k;
 	double volN,volC,phiU,phiN,phiCU,phiCe,phiC;
 	double dummy;
	
	ptrSize=zone->Size;
 
 	NI=*(ptrSize+0);
 	NJ=*(ptrSize+1);
 	NK=*(ptrSize+2);
 	
	k=index[0];
	j=index[1];
	i=index[2];
	
 	//printf("face=%d\t,index1=%d\t,index2=%d\tindex3=%d\n",face,index[0],index[1],index[2]);
 	switch(face)
	{
		case 1:
			
			if(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)<=0.0)
			{
				volN=*(slabvolume[0]+(NJ-1)*k+j);
				volC=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i);

				phiN=*(slabPhi[0]+(NJ-1)*k+j);
				phiC=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
				phiU=*(slabPhi[0]+(NJ-1)*k+j);
				phiCe=(phiN*volC+phiC*volN)/(volN+volC);
				
				phiCU=*(slabPhi[0]+(NJ-1)*k+j);
				
				
				*(FDn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCe-phiU))*Alpha;
				*(FUn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCU);

				dummy=0.0;				
				NCoeff[1]=0.0;
				
			}
			 else if(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)>0.0)				
			{
				volN=*(slabvolume[0]+(NJ-1)*k+j);
				volC=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
				
				phiN=*(slabPhi[0]+(NJ-1)*k+j);
				phiC=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
				phiU=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
				phiCe=(phiN*volC+phiC*volN)/(volN+volC);
				phiCU=*(phicellNew+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
								
				*(FDn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCe-phiU))*Alpha;
				*(FUn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCU);

				dummy=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face);				
				NCoeff[1]=0.0;
			
			}
		break;
		
		case 0:
			if(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)>0.0)
			{
				volC=*(slabvolume[1]+(NJ-1)*k+j);
				volN=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i);

				phiC=*(slabPhi[1]+(NJ-1)*k+j);
				phiN=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
				phiU=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
				phiCe=(phiN*volC+phiC*volN)/(volN+volC);
				phiCU=*(phicellNew+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
				
				
				*(FDn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCe-phiU))*Alpha;
				*(FUn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCU);

				dummy=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face);				
				NCoeff[0]=0.0;
			}
			else if(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)<=0.0)				
			{
				volN=*(slabvolume[1]+(NJ-1)*k+j);
				volC=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
								
				phiN=*(slabPhi[1]+(NJ-1)*k+j);
				phiC=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
				phiU=*(slabPhi[1]+(NJ-1)*k+j);
				phiCe=(phiN*volC+phiC*volN)/(volN+volC);
				phiCU=*(slabPhi[1]+(NJ-1)*k+j);
				
				
				*(FDn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCe-phiU))*Alpha;
				*(FUn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCU);

				dummy=0.0;
				NCoeff[0]=0.0;
			}	
		break;
	
		case 3:
			if(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)<=0.0)
			{
				volN=*(slabvolume[2]+(NI-1)*k+i);
				volC=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i);

				phiN=*(slabPhi[2]+(NI-1)*k+i);
				phiC=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
				phiU=*(slabPhi[2]+(NI-1)*k+i);
				phiCe=(phiN*volC+phiC*volN)/(volN+volC);
				phiCU=*(slabPhi[2]+(NI-1)*k+i);
				
				
				*(FDn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCe-phiU))*Alpha;
				*(FUn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCU);

				dummy=0.0;				
				NCoeff[3]=0.0;
			}
			else if(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)>0.0)				
			{
				volN=*(slabvolume[2]+(NI-1)*k+i);
				volC=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
							
				phiN=*(slabPhi[2]+(NI-1)*k+i);
				phiC=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
				phiU=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
				phiCe=(phiN*volC+phiC*volN)/(volN+volC);
				phiCU=*(phicellNew+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
				
				
				*(FDn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCe-phiU))*Alpha;
				*(FUn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCU);

				dummy=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face);				
				NCoeff[3]=0.0;
			
			}	
		break;
 
		case 2:
			if(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)>=0.0)
			{
				volN=*(slabvolume[3]+(NI-1)*k+i);
				volC=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i);

				phiN=*(slabPhi[3]+(NI-1)*k+i);
				phiC=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
				phiU=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
				phiCe=(phiN*volC+phiC*volN)/(volN+volC);
				phiCU=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
								
				*(FDn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCe-phiU))*Alpha;
				*(FUn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCU);

				dummy=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face);				
				NCoeff[2]=0.0;
			}
			else if(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)<0.0)				
			{
				volN=*(slabvolume[3]+(NI-1)*k+i);
				volC=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
				
				phiN=*(slabPhi[3]+(NI-1)*k+i);
				phiC=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
				phiU=*(slabPhi[3]+(NI-1)*k+i);
				phiCe=(phiN*volC+phiC*volN)/(volN+volC);
				phiCU=*(slabPhi[3]+(NI-1)*k+i);
				
				
				*(FDn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCe-phiU))*Alpha;
				*(FUn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCU);

				dummy=0.0;
				NCoeff[2]=0.0;
			
			}	
		break;
 
 		case 5:
 			*(FDn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=0.0;
			*(FUn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=0.0;

			dummy=0.0;				
			NCoeff[5]=0.0;
			
		break;
 		
		case 4:
			*(FDn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=0.0;
			*(FUn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=0.0;

			dummy=0.0;				
			NCoeff[4]=0.0;
			
		break;
	 
 	}

	return(dummy);
}
/**************************************************************************************************/
