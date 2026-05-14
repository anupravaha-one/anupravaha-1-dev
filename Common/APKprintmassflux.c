#include "APKprintmassflux.h"
/**
This is also a checking temporary function, which prints centroidal coordinates of a cell and  volume flux values of its six faces
into a data file.
*/
int printmassflux(Base *CgnsBase)
{
	FILE *flux;
	flux=fopen("SOLUTIONS/massfluxes.dat","w");
	int NI,NJ,NK,nZones,numZones;
	int i,j,k,*ptrSize;
	double *coordX,*coordY,*coordZ,*F;
	double *ptrPoint1X,*ptrPoint1Y,*ptrPoint1Z,centroid[3];

	nZones=((CgnsBase)->Number_of_zones);
	for(numZones=0;numZones<nZones;numZones++)
	{
		fprintf(flux,"~~~~~~~ZONENAME :%s ~~~~~~~~\n",((CgnsBase)->CgnsZone+numZones)->zonename);
		
		ptrSize=((CgnsBase)->CgnsZone+numZones)->Size;
      
		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);
		
		coordX=((CgnsBase->CgnsZone+numZones)->CgnsGrid)->X;
		coordY=((CgnsBase->CgnsZone+numZones)->CgnsGrid)->Y; 
		coordZ=((CgnsBase->CgnsZone+numZones)->CgnsGrid)->Z; 
		F = ((CgnsBase->CgnsZone+numZones)->Flowsol)->Massflux;
		
		for(k=0;k<NK-1;k++)
		for(j=0;j<NJ-1;j++)
		for(i=0;i<NI-1;i++)	
		{
			ptrPoint1X=coordX+NI*NJ*k+NI*j+i;
			ptrPoint1Y=coordY+NI*NJ*k+NI*j+i;
			ptrPoint1Z=coordZ+NI*NJ*k+NI*j+i;
			
			FindCellCentre(ptrSize,ptrPoint1X,ptrPoint1Y,ptrPoint1Z,centroid);
			if(((CgnsBase->CgnsZone+numZones)->Flowsol)->u!=NULL)
			{
				fprintf(flux,"\n%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\n",centroid[0],centroid[1],centroid[2],*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+0),*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+1),*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+2),*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+3),*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+4),*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+5));
			}
			
		}/*End of NI,NJ,NK loop*/
	}/*End of looping over Zones*/
	fclose(flux);
	return 0;
}

