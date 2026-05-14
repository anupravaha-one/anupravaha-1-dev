#include"APKInitializeU.h"
/**
This function is used for calculating the fully developed velocity for the case of Hartmann flow using the expression (4.5) from the book titled "Magnetofluiddynamics in channels and containers" by U.Muller & L.Buhler. This is used only for checking the profiles given in the same book.
*/
int InitializeU(Base *base)
{
	int NI,NJ,NK,nZones,numZones;
	int i,j,k,*ptrSize;
	double *coordX,*coordY,*coordZ;
	double *u,*prevu,*umass,*prevumass;
	double *ptrPoint1X,*ptrPoint1Y,*ptrPoint1Z,centroid[3];
	double Ha=5.0,ucap;
	Zone* zone;
	//FILE* Uvelocity;
	
	int FindCellCentre();
	
	//Uvelocity = fopen("SOLUTIONS/Uvelocity.dat","w");
	nZones=((base)->Number_of_zones);
	ucap=(1/(Ha)*(Ha));//for conducting walls case
	for(numZones=0;numZones<nZones;numZones++)
	{
		zone=(base->CgnsZone+numZones);

		ptrSize=zone->Size;

		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2); 
		
		coordX=((base->CgnsZone+numZones)->CgnsGrid)->X;
		coordY=((base->CgnsZone+numZones)->CgnsGrid)->Y; 
		coordZ=((base->CgnsZone+numZones)->CgnsGrid)->Z; 
		u=((base->CgnsZone+numZones)->Flowsol)->u;
		prevu=((base->CgnsZone+numZones)->PrevFlowsol)->u;
		umass=((base->CgnsZone+numZones)->Flowsol)->umass;
		prevumass=((base->CgnsZone+numZones)->PrevFlowsol)->umass;
		
		for(k=0;k<NK-1;k++)
		for(j=0;j<NJ-1;j++)
		for(i=0;i<NI-1;i++)
		{
			ptrPoint1X=coordX+NI*NJ*k+NI*j+i;
			ptrPoint1Y=coordY+NI*NJ*k+NI*j+i;
			ptrPoint1Z=coordZ+NI*NJ*k+NI*j+i;
			
			FindCellCentre(ptrSize,ptrPoint1X,ptrPoint1Y,ptrPoint1Z,centroid);
			*(u+(NI-1)*(NJ-1)*k+(NI-1)*j+i)=ucap *(1- (cosh(Ha*centroid[1])/cosh(Ha)));
			*(prevu+(NI-1)*(NJ-1)*k+(NI-1)*j+i)=ucap *(1- (cosh(Ha*centroid[1])/cosh(Ha)));
			*(umass+(NI-1)*(NJ-1)*k+(NI-1)*j+i)=ucap *(1- (cosh(Ha*centroid[1])/cosh(Ha)));
			*(prevumass+(NI-1)*(NJ-1)*k+(NI-1)*j+i)=ucap *(1- (cosh(Ha*centroid[1])/cosh(Ha)));
			
			//fprintf(Uvelocity,"%lf\t%lf\n",centroid[1],*(u+(NI-1)*(NJ-1)*k+(NI-1)*j+i));
		}
		
	}
	
	return 0;
}

