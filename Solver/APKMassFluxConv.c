#include "APKMassFluxConv.h"

/**
This function checks for the massflux conservation in the domain, i.e. total amount coming in and
total amount leaving the domain. This function sums up the massflux at all the faces in all zones
and returns the total difference between the incoming and outgoing fluxes.
*/

double MassFluxConv(Base *base)
{
	int iface,NI,NJ,NK,*ptrSize,nZone;
	int i,j,k,count,dir,index;
	char tag,name[33];
	void *ptr;
	double TotalFlux,*F;
	Zone *zone;
	Interface *facePtr;
	ZoneBoundary *bcPtr;
	TotalFlux=0.0;
	count=base->Number_of_zones;
	
	for(nZone=0;nZone<base->Number_of_zones;nZone++)
	{
		zone=(base->CgnsZone+nZone);
		ptrSize=zone->Size;
		F = zone->Flowsol->Massflux; 
		NI = *(ptrSize+0);
		NJ = *(ptrSize+1);
		NK = *(ptrSize+2);
		for(iface=0;iface < 6;iface++)
		{
			
			ptr = (void *)FindInterfaceBCPtrMod(zone,iface+1,name,&tag);
			if(tag=='i')
			{
				facePtr=(Interface *)ptr;
			}
			else
			if(tag=='b')
			{
				bcPtr=(ZoneBoundary *)ptr;
				dir=bcPtr->Direction;
				switch(dir)
				{
					case 1:
						for(k=0;k<NK-1;k++)
						for(j=0;j<NJ-1;j++)
						{
							i=NI-2;
							index = ((NJ-1)*(NI-1)*k+(NI-1)*j+i);
							TotalFlux=TotalFlux + *(F+(index)*6+1);
							//printf("case1 =%lf\n",*(F+(index)*6+1));
						}
					break;	
					case 2:
						for(k=0;k<NK-1;k++)
						for(j=0;j<NJ-1;j++)
						{
							i=0;
							index = ((NJ-1)*(NI-1)*k+(NI-1)*j+i);
							TotalFlux=TotalFlux + *(F+(index)*6+0);
							//printf("case2 =%lf\n",*(F+(index)*6+0));
						}
					break;
					case 3:
						for(k=0;k<NK-1;k++)
						for(i=0;i<NI-1;i++)
						{
							j=NJ-2;
							index = ((NJ-1)*(NI-1)*k+(NI-1)*j+i);
							TotalFlux=TotalFlux + *(F+(index)*6+3);
							//printf("case3 =%lf\n",*(F+(index)*6+3));
						}
					break;
					case 4:
						for(k=0;k<NK-1;k++)
						for(i=0;i<NI-1;i++)
						{
							j=0;
							index = ((NJ-1)*(NI-1)*k+(NI-1)*j+i);
							TotalFlux=TotalFlux + *(F+(index)*6+2);
							//printf("case4 =%lf\n",*(F+(index)*6+2));
						}
					break;
					case 5:
						for(j=0;j<NJ-1;j++)
						for(i=0;i<NI-1;i++)
						{
							k=NK-2;
							index = ((NJ-1)*(NI-1)*k+(NI-1)*j+i);
							TotalFlux=TotalFlux + *(F+(index)*6+5);
							//printf("case5 =%lf\n",*(F+(index)*6+5));
						}
					break;
					case 6:
						for(j=0;j<NJ-1;j++)
						for(i=0;i<NI-1;i++)
						{
							k=0;
							index = ((NJ-1)*(NI-1)*k+(NI-1)*j+i);
							TotalFlux=TotalFlux + *(F+(index)*6+4);
							//printf("case6 =%lf\n",*(F+(index)*6+4));
						}
					break;
				}
				printf("Flux at face No. = %d\tZone number%d\t is =%e\n",iface,nZone,TotalFlux);
			}
		}
	}
	//printf("Value before Return=%12.8e\n",TotalFlux);
	//outlet(base);
	return (TotalFlux);
}
	


