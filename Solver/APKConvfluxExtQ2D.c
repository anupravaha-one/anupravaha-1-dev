#include "APKConvfluxExtQ2D.h"
/**
This function computes the convective flux for QUICK scheme for 2 dimensional geometries for those cells
which have only one boundary or interface among its four faces.
This function returns both the Upwind flux(calculated using current iteration values) and the Deferred correction flux
(which is calculated using the previous iteration values). See Eqn(2.21) in Project Report No 1 Dec 2005, for further details.
*/ 
 
double ConvfluxExtQ2D(Zone *zone,double *F,double *phicellNew,double *phicell,double *volume,double *FDn,double *FUn,int index[3],double NCoeff[6],double *slabvolume[6],double *slabPhi[6],int face,double Alpha)
{ 
	int NI,NJ,NK,*ptrSize;
	int i,j,k;
	double volU,volUU,volD,phiU,phiUU,phiD,phiCU,phiQ;
	double dummy;
	
	ptrSize=zone->Size; 
	NI=*(ptrSize+0);
	NJ=*(ptrSize+1);
	NK=*(ptrSize+2);
	
	k=index[0];
	j=index[1];
	i=index[2];
	
	switch(face)
	{
		case 1:
			if(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)>0.0)
			{
				if(i!=0)
				{
					volU=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
					volUU=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i-1);
					volD=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i+1);

					phiU=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
					phiUU=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i-1);
					phiD=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i+1);

					phiCU=*(phicellNew+(NJ-1)*(NI-1)*k+(NI-1)*j+i);

					phiQ=QVIntpo(volU,volUU,volD,phiU,phiUU,phiD);

					*(FDn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiQ-phiU)*Alpha;
					*(FUn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCU);	

					dummy=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face);
					NCoeff[1]=0.0;
				}
				else
				{
					volU=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
					volUU=*(slabvolume[1]+(NJ-1)*k+j);
					volD=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i+1);

					phiU=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
					phiUU=*(slabPhi[1]+(NJ-1)*k+j);
					phiD=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i+1);

					phiCU=*(phicellNew+(NJ-1)*(NI-1)*k+(NI-1)*j+i);

					phiQ=QVIntpo(volU,volUU,volD,phiU,phiUU,phiD);

					*(FDn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiQ-phiU)*Alpha;
					*(FUn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCU);	

					dummy=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face);
					NCoeff[1]=0.0;
				}
			}
			else if(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)<=0.0)
			{
				if(i!=NI-3)
				{
					volU=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i+1);
					volUU=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i+2);
					volD=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i);

					phiU=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i+1);
					phiUU=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i+2);
					phiD=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
					
					phiCU=*(phicellNew+(NJ-1)*(NI-1)*k+(NI-1)*j+i+1);
					
					phiQ=QVIntpo(volU,volUU,volD,phiU,phiUU,phiD);
					
					*(FDn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiQ-phiU)*Alpha;
					*(FUn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCU);	
				
					dummy=0.0;				
					NCoeff[1]=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face);				
				}
				else
				{
					volU=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i+1);
					volUU=*(slabvolume[0]+(NJ-1)*k+j);
					volD=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i);

					phiU=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i+1);
					phiUU=*(slabPhi[0]+(NJ-1)*k+j);
					phiD=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
					
					phiCU=*(phicellNew+(NJ-1)*(NI-1)*k+(NI-1)*j+i+1);
					
					phiQ=QVIntpo(volU,volUU,volD,phiU,phiUU,phiD);
					
					*(FDn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiQ-phiU)*Alpha;
					*(FUn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCU);
					
					dummy=0.0;				
					NCoeff[1]=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face);					
				}
			}
		break;
		
		case 0:
			if(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)<=0.0)
			{
				if(i!=1)
				{
					volU=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i-1);
					volUU=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i-2);
					volD=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i);

					phiU=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i-1);
					phiUU=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i-2);
					phiD=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i);

					phiCU=*(phicellNew+(NJ-1)*(NI-1)*k+(NI-1)*j+i-1);

					phiQ=QVIntpo(volU,volUU,volD,phiU,phiUU,phiD);

					*(FDn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiQ-phiU)*Alpha;
					*(FUn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCU);	

					dummy=0.0;				
					NCoeff[0]=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face);
				}
				else
				{
					volU=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i-1);
					volUU=*(slabvolume[1]+(NJ-1)*k+j);
					volD=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i);

					phiU=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i-1);
					phiUU=*(slabPhi[1]+(NJ-1)*k+j);
					phiD=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i);

					phiCU=*(phicellNew+(NJ-1)*(NI-1)*k+(NI-1)*j+i-1);

					phiQ=QVIntpo(volU,volUU,volD,phiU,phiUU,phiD);

					*(FDn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiQ-phiU)*Alpha;
					*(FUn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCU);	

					dummy=0.0;				
					NCoeff[0]=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face);
				}
			}
			else if(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)>0.0)
			{
				if(i!=NI-2)
				{
					volU=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
					volUU=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i+1);
					volD=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i-1);

					phiU=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
					phiUU=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i+1);
					phiD=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i-1);
					
					phiCU=*(phicellNew+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
					
					phiQ=QVIntpo(volU,volUU,volD,phiU,phiUU,phiD);
					
					*(FDn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiQ-phiU)*Alpha;
					*(FUn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCU);	
				
					dummy=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face);				
					NCoeff[0]=0.0;	
				}
				else
				{
					volU=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
					volUU=*(slabvolume[0]+(NJ-1)*k+j);
					volD=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i-1);

					phiU=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
					phiUU=*(slabPhi[0]+(NJ-1)*k+j);
					phiD=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i-1);
					
					phiCU=*(phicellNew+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
					
					phiQ=QVIntpo(volU,volUU,volD,phiU,phiUU,phiD);
					
					*(FDn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiQ-phiU)*Alpha;
					*(FUn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCU);	
				
					dummy=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face);				
					NCoeff[0]=0.0;
				}		
				
			}
		break;
	
		case 3:
			if(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)>0.0)
			{
				if(j!=0)
				{
					volU=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
					volUU=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i);
					volD=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+i);

					phiU=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
					phiUU=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i);
					phiD=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+i);

					phiCU=*(phicellNew+(NJ-1)*(NI-1)*k+(NI-1)*j+i);

					phiQ=QVIntpo(volU,volUU,volD,phiU,phiUU,phiD);

					*(FDn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiQ-phiU)*Alpha;
					*(FUn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCU);	

					dummy=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face);				
					NCoeff[3]=0.0;
				}
				else
				{
					volU=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
					volUU=*(slabvolume[3]+(NI-1)*k+i);
					volD=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+i);

					phiU=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
					phiUU=*(slabPhi[3]+(NI-1)*k+i);
					phiD=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+i);

					phiCU=*(phicellNew+(NJ-1)*(NI-1)*k+(NI-1)*j+i);

					phiQ=QVIntpo(volU,volUU,volD,phiU,phiUU,phiD);

					*(FDn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiQ-phiU)*Alpha;
					*(FUn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCU);	

					dummy=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face);				
					NCoeff[3]=0.0;
				}

			}
			else if(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)<=0.0)
			{
				if(j!=NJ-3)
				{	
					volU=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+i);
					volUU=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*(j+2)+i);
					volD=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i);

					phiU=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+i);
					phiUU=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*(j+2)+i);
					phiD=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
					
					phiCU=*(phicellNew+(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+i);
					
					phiQ=QVIntpo(volU,volUU,volD,phiU,phiUU,phiD);
					
					*(FDn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiQ-phiU)*Alpha;
					*(FUn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCU);	
				
					dummy=0.0;				
					NCoeff[3]=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face);	
				}
				else
				{	
					volU=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+i);
					volUU=*(slabvolume[2]+(NI-1)*k+i);
					volD=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i);

					phiU=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+i);
					phiUU=*(slabPhi[2]+(NI-1)*k+i);
					phiD=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
					
					phiCU=*(phicellNew+(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+i);
					
					phiQ=QVIntpo(volU,volUU,volD,phiU,phiUU,phiD);
					
					*(FDn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiQ-phiU)*Alpha;
					*(FUn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCU);	
				
					dummy=0.0;				
					NCoeff[3]=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face);
				}			
				
			}
		break;
 
		case 2:
			if(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)<=0.0)
			{
				if(j!=1)
				{
					volU=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i);
					volUU=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*(j-2)+i);
					volD=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i);

					phiU=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i);
					phiUU=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*(j-2)+i);
					phiD=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i);

					phiCU=*(phicellNew+(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i);

					phiQ=QVIntpo(volU,volUU,volD,phiU,phiUU,phiD);

					*(FDn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiQ-phiU)*Alpha;
					*(FUn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCU);	

					dummy=0.0;			
					NCoeff[2]=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face);
				}
				else
				{
					volU=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i);
					volUU=*(slabvolume[3]+(NI-1)*k+i);
					volD=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*j+i);

					phiU=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i);
					phiUU=*(slabPhi[3]+(NI-1)*k+i);
					phiD=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*j+i);

					phiCU=*(phicellNew+(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i);

					phiQ=QVIntpo(volU,volUU,volD,phiU,phiUU,phiD);

					*(FDn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiQ-phiU)*Alpha;
					*(FUn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCU);	

					dummy=0.0;			
					NCoeff[2]=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face);
				}
			}
			else if(*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)>0.0)
			{
				if(j!=NJ-2)	
				{
					volU=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*(j)+i);
					volUU=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+i);
					volD=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i);

					phiU=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*(j)+i);
					phiUU=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+i);
					phiD=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i);

					phiCU=*(phicellNew+(NJ-1)*(NI-1)*k+(NI-1)*(j)+i);

					phiQ=QVIntpo(volU,volUU,volD,phiU,phiUU,phiD);

					*(FDn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiQ-phiU)*Alpha;
					*(FUn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCU);	

					dummy=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face);				
					NCoeff[2]=0.0;	
				}
				else
				{
					volU=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*(j)+i);
					volUU=*(slabvolume[2]+(NI-1)*k+i);
					volD=*(volume+(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i);

					phiU=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*(j)+i);
					phiUU=*(slabPhi[2]+(NI-1)*k+i);
					phiD=*(phicell+(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i);

					phiCU=*(phicellNew+(NJ-1)*(NI-1)*k+(NI-1)*(j)+i);

					phiQ=QVIntpo(volU,volUU,volD,phiU,phiUU,phiD);

					*(FDn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiQ-phiU)*Alpha;
					*(FUn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)*(phiCU);	

					dummy=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face);				
					NCoeff[2]=0.0;
				}	
				
			}
		break;
 
 		case 5:
 			
			*(FDn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=0.0;
			*(FUn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=0.0;

			dummy=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face);				
			NCoeff[5]=0.0;
				
			
		break;
 		
		case 4:
			

			*(FDn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=0.0;
			*(FUn+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=0.0;

			dummy=0.0;			
			NCoeff[4]=*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face);

		break;
	 
 	}

	return(dummy);
}
/**************************************************************************************************/


