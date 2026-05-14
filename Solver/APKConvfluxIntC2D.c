#include "APKConvfluxIntC2D.h"
/**
This function computes the convective flux for central difference scheme for 2 dimensional geometries for internal cells.
This function returns both the Upwind flux(calculated using current iteration values) and the Deferred correction flux
(which is calculated using the previous iteration values). See Eqn(2.21) in Project Report No 1 Dec 2005, for further details.
*/ 
int ConvfluxIntC2D(Zone *zone,double *F,double *phicellNew,double *phicell,double *volume,double *FDn,double *FUn,int index[3],double NCoeff[6],double Diagonal[6],double *slabvolume[6],double *slabPhi[6],double Alpha)
{ 
 	int face,NI,NJ,NK,*ptrSize;
	int i,j,k;
 	double volU,volD,phiU,phiD,phiCU,phiC;
	ptrSize=zone->Size;
 
 	NI=*(ptrSize+0);
 	NJ=*(ptrSize+1);
 	NK=*(ptrSize+2);
 	
	k=index[0];
	j=index[1];
	i=index[2];
	
 for(face=0;face<6;face++)
 {
 	switch(face)
	{
	
		case 1:/*east*/
			if(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)>0.0)
			{
				volU=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
				volD=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i+1);

				phiU=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
				phiD=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i+1);
				
				
				
				phiC=(volU*phiD+volD*phiU)/(volU+volD);
				phiCU=*(phicellNew+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
				
							
				*(FDn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiC-phiU))*Alpha;
				*(FUn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCU);	
				
				Diagonal[1]=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face);
				NCoeff[1]=0.0;
			}
			else if(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)<=0.0)
			{
					volU=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i+1);
					volD=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i);

					phiU=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i+1);
					phiD=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
					phiC=(volU*phiD+volD*phiU)/(volU+volD);
					
					phiCU=*(phicellNew+(NJ-1)*(NI-1)*k+(NI-1)*j+i+1);
				
					
				*(FDn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiC-phiU))*Alpha;
				*(FUn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCU);	
				
					Diagonal[1]=0.0;				
					NCoeff[1]=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face);				
			}
		break;
		
		case 0:/*west*/
			if(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)>=0.0)
			{
				volU=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
				volD=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i-1);

				phiU=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
				phiD=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i-1);
				phiC=(volU*phiD+volD*phiU)/(volU+volD);
				
				phiCU=*(phicellNew+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
				
				
				*(FDn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiC-phiU))*Alpha;
				*(FUn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCU);	

				Diagonal[0]=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face);				
				NCoeff[0]=0.0;

			}
			else if(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)<0.0)
			{
				volU=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i-1);
				volD=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i);

				phiU=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i-1);
				phiD=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
				phiC=(volU*phiD+volD*phiU)/(volU+volD);
				
				phiCU=*(phicellNew+(NJ-1)*(NI-1)*k+(NI-1)*j+i-1);
				
				
				*(FDn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=Alpha*(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiC-phiU));
				*(FUn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCU);	
				
				Diagonal[0]=0.0;				
				NCoeff[0]=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face);		
			}
		break;
	
		case 3:
			if(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)>0.0)
			{
				volU=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
				volD=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+i);

				phiU=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
				phiD=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+i);
				phiC=(volU*phiD+volD*phiU)/(volU+volD);
				phiCU=*(phicellNew+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
				
				
				*(FDn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiC-phiU))*Alpha;
				*(FUn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCU);	

				Diagonal[3]=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face);				
				NCoeff[3]=0.0;

			}
			else if(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)<=0.0)
			{
				volU=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+i);
				volD=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i);

				phiU=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+i);
				phiD=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
				phiC=(volU*phiD+volD*phiU)/(volU+volD);
				phiCU=*(phicellNew+(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+i);
				
				
				*(FDn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiC-phiU))*Alpha;
				*(FUn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCU);	

				Diagonal[3]=0.0;				
				NCoeff[3]=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face);	
			}
		break;
 
		case 2:
			if(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)>=0.0)
			{
				volU=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*(j)+i);
				volD=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i);

				phiU=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*(j)+i);
				phiD=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i);
				phiC=(volU*phiD+volD*phiU)/(volU+volD);
				
				phiCU=*(phicellNew+(NJ-1)*(NI-1)*k+(NI-1)*(j)+i);
				
				
				*(FDn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiC-phiU))*Alpha;
				*(FUn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCU);	

				Diagonal[2]=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face);			
				NCoeff[2]=0.0;

			}
			else if(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)<0.0)
			{
				volU=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i);
				volD=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*(j)+i);

				phiU=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i);
				phiD=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*(j)+i);
				phiC=(volU*phiD+volD*phiU)/(volU+volD);
				phiCU=*(phicellNew+(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i);
				
				
				*(FDn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiC-phiU))*Alpha;
				*(FUn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCU);	

				Diagonal[2]=0.0;			
				NCoeff[2]=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face);
			}
		break;
 
 		case 5:
 			*(FDn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=0.0;
			*(FUn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=0.0;
				
			Diagonal[5]=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face);				
			NCoeff[5]=0.0;
				
			
		break;
 		
		case 4:
			
			*(FDn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=0.0;
			*(FUn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=0.0;

			Diagonal[4]=0.0;			
			NCoeff[4]=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face);
				
			
		break;
	 
 	}
}
	return(0);
}
/**************************************************************************************************/

