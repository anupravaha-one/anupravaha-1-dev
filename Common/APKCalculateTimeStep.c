#include"APKCalculateTimeStep.h"

/**
This function is used to calculate timestep using the CFL criteria, if user does not want to define time step
and wants the code to fix timestep approproately.

This function comes in picture, only when the user has given the option to the code to 
find the time step from CFL criteria. For the calculation of time step, Courant-Fricedrichs-lewy formulation has 
been used, present function is wholly based on the Book by K. Murlidhar and T. Sundararajan (Page. 465).
*/
double CalculateTimeStep(Base *base)
{
	int nZone,numZone,*ptrSize;
	GridCoordinate *Zgrid;
	Zone *zone;
	double CFL=1.0;
	double time;
	double maxNU=0,maxNV=0,maxNW=0,maxTXC,maxTYC,maxTZC,maxTXD,maxTYD,maxTZD;
	double maxU=0,maxV=0,maxW=0,timeSc,ConvecTime,timeSd,diffTimeScale;
	double *maxTimeX,*maxTimeY,*maxTimeZ;
     
	nZone=base->Number_of_zones;
      
	maxTimeX=(double *)malloc(nZone*2*sizeof(double));
	maxTimeY=(double *)malloc(nZone*2*sizeof(double));
	maxTimeZ=(double *)malloc(nZone*2*sizeof(double));
	
	for(numZone=0;numZone<nZone;numZone++)
	{
		zone=(base->CgnsZone+numZone);
		Zgrid=zone->CgnsGrid;
		ptrSize=zone->Size;
		
		//Nu=zone->Flowsol->Viscoeff;
		FindmaXConv(zone,ptrSize,&maxU,&maxV,&maxW,base->BaseParameters->flag2d);
		
		FindmaXDiff(zone,ptrSize,&maxNU,&maxNV,&maxNW,base->BaseParameters->flag2d);
		
		*(maxTimeX+numZone*2+0)=maxU;
		*(maxTimeY+numZone*2+0)=maxV;
		*(maxTimeZ+numZone*2+0)=maxW;
		
		*(maxTimeX+numZone*2+1)=maxNU;
		*(maxTimeY+numZone*2+1)=maxNV;
		*(maxTimeZ+numZone*2+1)=maxNW; 	
	}
	
	maxTXC=*(maxTimeX+0);
	maxTYC=*(maxTimeY+0);
	maxTZC=*(maxTimeZ+0);
	
	maxTXD=*(maxTimeX+1);
	maxTYD=*(maxTimeY+1);
	maxTZD=*(maxTimeZ+1);
	
	for(numZone=0;numZone<nZone;numZone++)
	{
		
		if(maxTXC<*(maxTimeX+numZone*2+0))
		maxTXC=*(maxTimeX+numZone*2+0);
		
		if(maxTYC<*(maxTimeY+numZone*2+0))
		maxTYC=*(maxTimeY+numZone*2+0);
	
		if(maxTZC<*(maxTimeZ+numZone*2+0))
		maxTZC=*(maxTimeZ+numZone*2+0);
		
		if(maxTXD<*(maxTimeX+numZone*2+1))
		maxTXD=*(maxTimeX+numZone*2+1);
		
		if(maxTYD<*(maxTimeY+numZone*2+1))
		maxTYD=*(maxTimeY+numZone*2+1);
	
		if(maxTZD<*(maxTimeZ+numZone*2+1))
		maxTZD=*(maxTimeZ+numZone*2+1);
	
	}
		timeSc=(maxTXC+maxTYC+maxTZC);
		ConvecTime=CFL/timeSc;
		//printf("Convectime=%lf\n",ConvecTime);
		
		timeSd=(maxTXD+maxTYD+maxTZD);
		//printf("%lf %lf %lf\n",maxTXD,maxTYD,maxTZD);
		diffTimeScale=0.5/(timeSd);
		//printf("difftime=%lf\n",timeSd);
		
		time=ConvecTime;
		if(time>diffTimeScale)
		time=diffTimeScale;
		
	free(maxTimeX);
	free(maxTimeY);
	free(maxTimeZ);
	return(time);      
}
      

