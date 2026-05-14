#include "APKCalcPsiSource.h"
int CalcPsiSource(Base *base)
{
	Zone *zone;
	int i,j,k,NI,NJ,NK,*ptrSize,indexP,iFace;//rajendra
	double *MagneticFlux,sumMagneticFlux;
	int iZone,numZones;
	numZones=base->Number_of_zones;
		if(base->BaseParameters->flag2d==1)
		{
			for(iZone=0;iZone<numZones;iZone++)//rajendra
			{
				zone=base->CgnsZone+iZone;
				ptrSize = (zone->grid+0)->GSize;
				NI=*(ptrSize+0);NJ=*(ptrSize+1);NK=*(ptrSize+2);
				
				for(k=0;k<NK-1;k++)
				for(j=0;j<NJ-1;j++)
				for(i=0;i<NI-1;i++)
				{
					indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				
					MagneticFlux=zone->Flowsol->Magneticflux+indexP*6;
					
					sumMagneticFlux=0.0;
					for(iFace=0;iFace<4;iFace++)
					{
						sumMagneticFlux=sumMagneticFlux +(*(MagneticFlux+iFace));
						
					}
					*((zone->grid)->GRHS+indexP) = sumMagneticFlux;
					
					
					
					
				}
			}
		}
		else
		{
			for(iZone=0;iZone<numZones;iZone++)//rajendra
			{
				zone=base->CgnsZone+iZone;
				ptrSize = (zone->grid+0)->GSize;
				NI=*(ptrSize+0);NJ=*(ptrSize+1);NK=*(ptrSize+2);
				indexP=0;
				for(k=0;k<NK-1;k++)
				for(j=0;j<NJ-1;j++)
				for(i=0;i<NI-1;i++)
				{
					//indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				
					MagneticFlux=zone->Flowsol->Magneticflux+indexP*6;
					sumMagneticFlux=0.0;
					for(iFace=0;iFace<6;iFace++)
					{
						sumMagneticFlux=sumMagneticFlux +(*(MagneticFlux+iFace));
					}
					*((zone->grid)->GRHS+indexP) = sumMagneticFlux;
					indexP++;
				}
			}
 		}
	return 0;			
}

