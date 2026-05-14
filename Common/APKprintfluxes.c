#include "APKprintfluxes.h"
/**
This is a temporary checking function, which prints volume flux values at six faces of a cell.
*/
int printfluxes(Base *CgnsBase)
{
	FILE *pflux;
	pflux = fopen("SOLUTIONS/printfluxes.dat","w");
	int NI,NJ,NK,nZones,numZones,face;
	int i,j,k,*ptrSize;
	double *volume,*surfarea,*massflux/*,surfacearea,diff*/;
	
	nZones=((CgnsBase)->Number_of_zones);
	for(numZones=0;numZones<nZones;numZones++)
	{
		ptrSize=((CgnsBase)->CgnsZone+numZones)->Size;
      
		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);
		
		for(k=0;k<NK-1;k++)
		for(j=0;j<NJ-1;j++)
		for(i=0;i<NI-1;i++)	
		{
			volume = (((CgnsBase)->CgnsZone+numZones)->ZGeom->Volume+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			surfarea = (((CgnsBase)->CgnsZone+numZones)->ZGeom->SurfaceArea+((NI-1)*(NJ-1)*k+(NI-1)*j+i)*6*3);
			massflux = (((CgnsBase)->CgnsZone+numZones)->Flowsol->Massflux+((NI-1)*(NJ-1)*k+(NI-1)*j+i)*6);
			fprintf(pflux,"k=%d j=%d i=%d volume=%12.8lf\n",k,j,i,*volume);
			for(face=0;face<6;face++)
			{
				fprintf(pflux,"k=%d j=%d i=%d face=%d umass=%12.8lf\n",k,j,i,face,*(massflux+face));
				fprintf(pflux,"k=%d j=%d i=%d face=%d surfarea[3*%d+0]=%lf surfarea[3*%d+0]=%lf surfarea[3*%d+0]=%lf\n",k,j,i,face,face,*(surfarea+3*face+0),face,*(surfarea+3*face+1),face,*(surfarea+3*face+2));
			}
		}
	}
	fclose(pflux);
	return 0;
}

			
