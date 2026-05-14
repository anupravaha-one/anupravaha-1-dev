#include "APKConvfluxExtCother3D.h"
/**
This function computes the convective flux for central difference scheme for 3 dimensional geometries for those cells
which have only one boundary or interface among its six faces.
This function returns both the Upwind flux(calculated using current iteration values) and the Deferred correction flux
(which is calculated using the previous iteration values). See Eqn(2.21) in Project Report No 1 Dec 2005, for further details.
*/ 
double ConvfluxExtCother3D(Zone *zone,double *F,double *phicellNew,double *phicell,double *volume,double *FDn,double *FUn,int index[3],double NCoeff[6],double *slabvolume[6],double *slabPhi[6],int face,double Alpha)
{ 
	int NI,NJ,NK,*ptrSize;
	int i,j,k;
	double volU,volD,phiU,phiD,phiCU,phiCe;
	double dummy;
	ptrSize=zone->Size;
 
	NI=*(ptrSize+0);
	NJ=*(ptrSize+1);
	NK=*(ptrSize+2);
	
	k=index[0];
	j=index[1];
	i=index[2];
	//printf("%d %d %d %d\n",k,j,i,face);
	switch(face)
	{
		case 1:
			if(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)>0.0)
			{
					volU=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
					volD=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i+1);

					phiU=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
					phiD=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i+1);

					phiCU=*(phicellNew+(NJ-1)*(NI-1)*k+(NI-1)*j+i);

					phiCe=(volU*phiD+volD*phiU)/(volU+volD);

					*(FDn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCe-phiU)*Alpha;
					*(FUn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCU);	

					dummy=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face);				
					NCoeff[1]=0.0;
				
			}
			else if(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)<=0.0)
			{
					volU=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i+1);
					volD=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i);

					phiU=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i+1);
					phiD=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
					
					phiCU=*(phicellNew+(NJ-1)*(NI-1)*k+(NI-1)*j+i+1);
					
					phiCe=(volU*phiD+volD*phiU)/(volU+volD);
					
					*(FDn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCe-phiU)*Alpha;
					*(FUn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCU);	
				
					dummy=0.0;				
					NCoeff[1]=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face);				
				
			}
		break;
		
		case 0:
			if(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)<=0.0)
			{
					volU=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i-1);
					volD=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i);

					phiU=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i-1);
					phiD=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i);

					phiCU=*(phicellNew+(NJ-1)*(NI-1)*k+(NI-1)*j+i-1);

					phiCe=(volU*phiD+volD*phiU)/(volU+volD);

					*(FDn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCe-phiU)*Alpha;
					*(FUn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCU);	

					dummy=0.0;				
					NCoeff[0]=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face);
				
			}
			else if(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)>0.0)
			{
				volU=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
				volD=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i-1);

				phiU=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
				phiD=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i-1);

				phiCU=*(phicellNew+(NJ-1)*(NI-1)*k+(NI-1)*j+i);

				phiCe=(volU*phiD+volD*phiU)/(volU+volD);

				*(FDn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCe-phiU)*Alpha;
				*(FUn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCU);	

				dummy=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face);
				NCoeff[0]=0.0;	
			
			}
		break;
	
		case 3:
			if(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)>0.0)
			{
				volU=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
				volD=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+i);

				phiU=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
				phiD=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+i);

				phiCU=*(phicellNew+(NJ-1)*(NI-1)*k+(NI-1)*j+i);

				phiCe=(volU*phiD+volD*phiU)/(volU+volD);

				*(FDn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCe-phiU)*Alpha;
				*(FUn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCU);	

				dummy=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face);				
				NCoeff[3]=0.0;
				
			}
			else if(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)<=0.0)
			{
				volU=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+i);
				volD=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i);

				phiU=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+i);
				phiD=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i);

				phiCU=*(phicellNew+(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+i);

				phiCe=(volU*phiD+volD*phiU)/(volU+volD);

				*(FDn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCe-phiU)*Alpha;
				*(FUn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCU);	

				dummy=0.0;				
				NCoeff[3]=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face);	
			}
		break;
 
		case 2:
			if(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)<=0.0)
			 {
				volU=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i);
				volD=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i);

				phiU=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i);
				phiD=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i);

				phiCU=*(phicellNew+(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i);

				phiCe=(volU*phiD+volD*phiU)/(volU+volD);

				*(FDn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCe-phiU)*Alpha;
				*(FUn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCU);	

				dummy=0.0;			
				NCoeff[2]=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face);
				
			}
			else if(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)>0.0)
			{
				volU=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*(j)+i);
				volD=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i);

				phiU=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*(j)+i);
				phiD=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i);

				phiCU=*(phicellNew+(NJ-1)*(NI-1)*k+(NI-1)*(j)+i);

				phiCe=(volU*phiD+volD*phiU)/(volU+volD);

				*(FDn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCe-phiU)*Alpha;
				*(FUn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCU);	

				dummy=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face);				
				NCoeff[2]=0.0;	
				
			}
		break;
 
 		case 5:
 			if(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)>0.0)
			{
				volU=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
				volD=*(volume+(NJ-1)*(NI-1)*(k+1)+(NI-1)*(j)+i);

				phiU=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
				phiD=*(phicell+(NJ-1)*(NI-1)*(k+1)+(NI-1)*(j)+i);

				phiCU=*(phicellNew+(NJ-1)*(NI-1)*k+(NI-1)*j+i);

				phiCe=(volU*phiD+volD*phiU)/(volU+volD);

				*(FDn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCe-phiU)*Alpha;
				*(FUn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCU);	

				dummy=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face);				
				NCoeff[5]=0.0;
			
			}
			else if(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)<=0.0)
			{
				
				volU=*(volume+(NJ-1)*(NI-1)*(k+1)+(NI-1)*(j)+i);
				volD=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i);

				phiU=*(phicell+(NJ-1)*(NI-1)*(k+1)+(NI-1)*(j)+i);
				phiD=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i);

				phiCU=*(phicellNew+(NJ-1)*(NI-1)*(k+1)+(NI-1)*(j)+i);

				phiCe=(volU*phiD+volD*phiU)/(volU+volD);

				*(FDn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCe-phiU)*Alpha;
				*(FUn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCU);	

				dummy=0.0;				
				NCoeff[5]=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face);	
			}
		break;
 		
		case 4:
			if(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)<=0.0)
			 {
				volU=*(volume+(NJ-1)*(NI-1)*(k-1)+(NI-1)*j+i);
				volD=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i);

				phiU=*(phicell+(NJ-1)*(NI-1)*(k-1)+(NI-1)*(j)+i);
				phiD=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i);

				phiCU=*(phicellNew+(NJ-1)*(NI-1)*(k-1)+(NI-1)*(j)+i);

				phiCe=(volU*phiD+volD*phiU)/(volU+volD);

				*(FDn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCe-phiU)*Alpha;
				*(FUn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCU);	

				dummy=0.0;			
				NCoeff[4]=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face);
				
			}
			else if(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)>0.0)
			{
				volU=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*(j)+i);
				volD=*(volume+(NJ-1)*(NI-1)*(k-1)+(NI-1)*(j)+i);

				phiU=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*(j)+i);
				phiD=*(phicell+(NJ-1)*(NI-1)*(k-1)+(NI-1)*(j)+i);

				phiCU=*(phicellNew+(NJ-1)*(NI-1)*k+(NI-1)*(j)+i);

				phiCe=(volU*phiD+volD*phiU)/(volU+volD);

				*(FDn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCe-phiU)*Alpha;
				*(FUn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCU);	

				dummy=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face);				
				NCoeff[4]=0.0;	
			
			}
		break;
	 
 	}

	return(dummy);
}
/**************************************************************************************************/


