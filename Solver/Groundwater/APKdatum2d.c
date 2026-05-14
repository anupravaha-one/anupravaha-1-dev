#include "APKdatum2d.h"

/**
This will calculate the flux due to datum height in case of saturated medium. It will multily hydraulic conductivity
by the surface area component in Y direction. Similar to this datum3d function is there for 3D case.

*/

int  datum2d(Zone *zone,double *F ,double K)
{
	
	int k,j,i;
	int Index;
	double *SurfArea;
	
	int NI,NJ,NK,*ptrSize,face;
	ptrSize=zone->Size;
	
	SurfArea=zone->ZGeom->SurfaceArea;
				
	NI= *(ptrSize+0);
	NJ= *(ptrSize+1);
	NK= *(ptrSize+2);
	
	for(k=0;k<NK-1;k++)
	{
	for(j=0;j<NJ-1;j++)
	{
	for(i=0;i<NI-1;i++)
		{  
	Index=((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6*3;
	
	//printf("\n %d",Index);
	face = 0;  // west face
		*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = K*(*(SurfArea+Index+3*0+1));
		//fprintf(fp1,"west=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
	face = 1;  // east face
		*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = K*(*(SurfArea+Index+3*1+1));
		//fprintf(fp1,"east=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
	face = 2;  // south face
		*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = K*(*(SurfArea+Index+3*2+1));
		//fprintf(fp1,"south=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
	face = 3;  // north face
		*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) =K*(*(SurfArea+Index+3*3+1));
		//fprintf(fp1,"north=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
	face = 4;  // bottom face
		*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = 0.0;
		//fprintf(fp1,"bottom=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
	face = 5;  // top face
		*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = 0.0;
		//fprintf(fp1,"top=%lf\n",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));;
		}
		}
		
		}
	
	return 0; 
	
}
		
			
	
		
int  datum3d(Zone *zone,double *F ,double K)
{
	int k,j,i;
	int Index;
	double *SurfArea;
	
	int NI,NJ,NK,*ptrSize,face;
	ptrSize=zone->Size;
	
	SurfArea=zone->ZGeom->SurfaceArea;
				
	NI= *(ptrSize+0);
	NJ= *(ptrSize+1);
	NK= *(ptrSize+2);
	
	for(k=0;k<NK-1;k++)
	{
	for(j=0;j<NJ-1;j++)
	{
	for(i=0;i<NI-1;i++)
		{  
	Index=((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6*3;
	
	//printf("\n %d",Index);
	face = 0;  // west face
		*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = K*(*(SurfArea+Index+3*0+1));
		//fprintf(fp1,"west=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
	face = 1;  // east face
		*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = K*(*(SurfArea+Index+3*1+1));
		//fprintf(fp1,"east=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
	face = 2;  // south face
		*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = K*(*(SurfArea+Index+3*2+1));
		//fprintf(fp1,"south=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
	face = 3;  // north face
		*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) =K*(*(SurfArea+Index+3*3+1));
		//fprintf(fp1,"north=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
	face = 4;  // bottom face
		*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = K*(*(SurfArea+Index+3*4+1));
		//fprintf(fp1,"bottom=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
	face = 5;  // top face
		*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = K*(*(SurfArea+Index+3*5+1));
		//fprintf(fp1,"top=%lf\n",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));;
		}
		}
		
		}
	
	return 0; 

		
}
