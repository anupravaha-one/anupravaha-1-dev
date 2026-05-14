#include "APKdatumUnsat2d.h"
/**
 This function is used to calculate the  datum flux in case for unsatuarted medium. It interpolate the unsaturated permeabilty at face and multiply it, by the surface area in Y direction for 2d. 
It has other function also datumUnsat3d() for calculating datum flux for 3D.
*/
       	   
int datumUnsat2d(Base *base,Zone *zone,double *volumes,double *slabVolume[6],double *diffusionCoeffdomain,double *slabDiffusionCoeff[6],double *F)
{        
       double UnSatPer_p[7],UnSatPer_face[6],vol_p,vol_w,vol_e,vol_s,vol_n,vol_b,vol_t;   
       double    *SurfArea;                                     
       int NI,NJ,NK,i,j,k,*ptrSize,face,Index;                                                          
	
	SurfArea=zone->ZGeom->SurfaceArea;                                                                            
	ptrSize = zone->Size;                                                                     
	                                                                                       
	
	
	NI=*(ptrSize+0);
	NJ=*(ptrSize+1);
	NK=*(ptrSize+2);
	 
	
         for(k=0;k<NK-1;k++)
       	 for(j=0;j<NJ-1;j++)
	 for(i=0;i<NI-1;i++)
	 {  
	 Index=((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6*3;
       
	 vol_p=*(volumes+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
	 UnSatPer_p[0]=*(diffusionCoeffdomain+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
	
	if(i!=0)
	{	vol_w=*(volumes+(NJ-1)*(NI-1)*k+(NI-1)*j+(i-1));
		UnSatPer_p[1]=*(diffusionCoeffdomain+(NJ-1)*(NI-1)*k+(NI-1)*j+(i-1));
	}
	else
	{	vol_w =*(slabVolume[1]+(NJ-1)*k+j);
		UnSatPer_p[1] =*(slabDiffusionCoeff[1]+(NJ-1)*k+j);
		
	}
	if(i!=(NI-2))
	{       vol_e=*(volumes+(NJ-1)*(NI-1)*k+(NI-1)*j+(i+1));
		UnSatPer_p[2]=*(diffusionCoeffdomain+(NJ-1)*(NI-1)*k+(NI-1)*j+(i+1));
	}
	else
	{
		vol_e=*(slabVolume[0]+(NJ-1)*k+j);
		UnSatPer_p[2]=*(slabDiffusionCoeff[0]+(NJ-1)*k+j);
	}
	if(j!=0)
	{	vol_s=*(volumes+(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i);
		UnSatPer_p[3]=*(diffusionCoeffdomain+(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i);
	}
	else
	{
		vol_s=*(slabVolume[3]+(NI-1)*k+i);
		UnSatPer_p[3]=*(slabDiffusionCoeff[3]+(NI-1)*k+i);
	}
	if(j!=(NJ-2))
	{
		vol_n=*(volumes+(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+i);
		UnSatPer_p[4]=*(diffusionCoeffdomain+(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+i);
	}
	else 
	{
		vol_n=*(slabVolume[2]+(NI-1)*k+i);
		UnSatPer_p[4]=*(slabDiffusionCoeff[2]+(NI-1)*k+i);
		
	}
	if(k!=0)
	{
		vol_b=*(volumes+(NJ-1)*(NI-1)*(k-1)+(NI-1)*j+i);
		UnSatPer_p[5]=*(diffusionCoeffdomain+(NJ-1)*(NI-1)*(k-1)+(NI-1)*j+i);
	}
	else
	{
	        vol_b=*(slabVolume[5]+(NI-1)*j+i);
		UnSatPer_p[5]=*(slabDiffusionCoeff[5]+(NI-1)*j+i);
	}
	if(k!=(NK-2))
	{
		vol_t=*(volumes+(NJ-1)*(NI-1)*(k+1)+(NI-1)*j+i);
		UnSatPer_p[6]=*(diffusionCoeffdomain+(NJ-1)*(NI-1)*(k+1)+(NI-1)*j+i);
	}
	else
	{
		vol_t=*(slabVolume[4]+(NI-1)*j+i);
		UnSatPer_p[6]=*(slabDiffusionCoeff[4]+(NI-1)*j+i);
	}
	
	
	
	UnSatPer_face[1] = (UnSatPer_p[2]*vol_p + UnSatPer_p[0]*vol_e)/(vol_p+vol_e);
	UnSatPer_face[0] = (UnSatPer_p[1]*vol_p + UnSatPer_p[0]*vol_w)/(vol_p+vol_w);
	UnSatPer_face[3] = (UnSatPer_p[4]*vol_p + UnSatPer_p[0]*vol_n)/(vol_p+vol_n);
	UnSatPer_face[2] = (UnSatPer_p[3]*vol_p + UnSatPer_p[0]*vol_s)/(vol_p+vol_s);
	UnSatPer_face[5] = (UnSatPer_p[6]*vol_p + UnSatPer_p[0]*vol_t)/(vol_p+vol_t);
	UnSatPer_face[4] = (UnSatPer_p[5]*vol_p + UnSatPer_p[0]*vol_b)/(vol_p+vol_b);
	
	face = 0;  // west face
	*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) =UnSatPer_face[0]*(*(SurfArea+Index+3*0+1));	
	
	face = 1;  // east face
	*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) =UnSatPer_face[1]*(*(SurfArea+Index+3*1+1));
	
	face = 2;  // south face
	*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) =UnSatPer_face[2]*(*(SurfArea+Index+3*2+1));
 				
 	face = 3;  // north face
 	*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) =UnSatPer_face[3]*(*(SurfArea+Index+3*3+1));
 	
 	face = 4;  // bottom face
 	*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = 0.0;
 			
 	face = 5;  // top face
 	*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = 0.0;
	}
	
	return(0);
	}
   
     
int datumUnsat3d(Base *base,Zone *zone,double *volumes,double *slabVolume[6],double *diffusionCoeffdomain,double *slabDiffusionCoeff[6],double *F)
  {        
       double UnSatPer_p[7],UnSatPer_face[6],vol_p,vol_w,vol_e,vol_s,vol_n,vol_b,vol_t;   
       double    *SurfArea;                                     
       int NI,NJ,NK,i,j,k,*ptrSize,face,Index;                                                          
	
         SurfArea=zone->ZGeom->SurfaceArea;                                                                            
	 ptrSize = zone->Size;                                                                     
	                                                                                       
	
	
	NI=*(ptrSize+0);
	NJ=*(ptrSize+1);
	NK=*(ptrSize+2);
	 
	
         for(k=0;k<NK-1;k++)
       	 for(j=0;j<NJ-1;j++)
	 for(i=0;i<NI-1;i++)
	 {  
	 Index=((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6*3;
       
	 vol_p=*(volumes+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
	 UnSatPer_p[0]=*(diffusionCoeffdomain+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
	
	if(i!=0)
	{	vol_w=*(volumes+(NJ-1)*(NI-1)*k+(NI-1)*j+(i-1));
		UnSatPer_p[1]=*(diffusionCoeffdomain+(NJ-1)*(NI-1)*k+(NI-1)*j+(i-1));
	}
	else
	{	vol_w =*(slabVolume[1]+(NJ-1)*k+j);
		UnSatPer_p[1] =*(slabDiffusionCoeff[1]+(NJ-1)*k+j);
		
	}
	if(i!=(NI-2))
	{       vol_e=*(volumes+(NJ-1)*(NI-1)*k+(NI-1)*j+(i+1));
		UnSatPer_p[2]=*(diffusionCoeffdomain+(NJ-1)*(NI-1)*k+(NI-1)*j+(i+1));
	}
	else
	{
		vol_e=*(slabVolume[0]+(NJ-1)*k+j);
		UnSatPer_p[2]=*(slabDiffusionCoeff[0]+(NJ-1)*k+j);
	}
	if(j!=0)
	{	vol_s=*(volumes+(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i);
		UnSatPer_p[3]=*(diffusionCoeffdomain+(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i);
	}
	else
	{
		vol_s=*(slabVolume[3]+(NI-1)*k+i);
		UnSatPer_p[3]=*(slabDiffusionCoeff[3]+(NI-1)*k+i);
	}
	if(j!=(NJ-2))
	{
		vol_n=*(volumes+(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+i);
		UnSatPer_p[4]=*(diffusionCoeffdomain+(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+i);
	}
	else 
	{
		vol_n=*(slabVolume[2]+(NI-1)*k+i);
		UnSatPer_p[4]=*(slabDiffusionCoeff[2]+(NI-1)*k+i);
		
	}
	if(k!=0)
	{
		vol_b=*(volumes+(NJ-1)*(NI-1)*(k-1)+(NI-1)*j+i);
		UnSatPer_p[5]=*(diffusionCoeffdomain+(NJ-1)*(NI-1)*(k-1)+(NI-1)*j+i);
	}
	else
	{
	        vol_b=*(slabVolume[5]+(NI-1)*j+i);
		UnSatPer_p[5]=*(slabDiffusionCoeff[5]+(NI-1)*j+i);
	}
	if(k!=(NK-2))
	{
		vol_t=*(volumes+(NJ-1)*(NI-1)*(k+1)+(NI-1)*j+i);
		UnSatPer_p[6]=*(diffusionCoeffdomain+(NJ-1)*(NI-1)*(k+1)+(NI-1)*j+i);
	}
	else
	{
		vol_t=*(slabVolume[4]+(NI-1)*j+i);
		UnSatPer_p[6]=*(slabDiffusionCoeff[4]+(NI-1)*j+i);
	}
	
	
	
	UnSatPer_face[1] = (UnSatPer_p[2]*vol_p + UnSatPer_p[0]*vol_e)/(vol_p+vol_e);// East face
	UnSatPer_face[0] = (UnSatPer_p[1]*vol_p + UnSatPer_p[0]*vol_w)/(vol_p+vol_w);// west face
	UnSatPer_face[3] = (UnSatPer_p[4]*vol_p + UnSatPer_p[0]*vol_n)/(vol_p+vol_n);// north face
	UnSatPer_face[2] = (UnSatPer_p[3]*vol_p + UnSatPer_p[0]*vol_s)/(vol_p+vol_s);// south face
	UnSatPer_face[5] = (UnSatPer_p[6]*vol_p + UnSatPer_p[0]*vol_t)/(vol_p+vol_t);// top face
	UnSatPer_face[4] = (UnSatPer_p[5]*vol_p + UnSatPer_p[0]*vol_b)/(vol_p+vol_b);// bottom face
		
	
	
	face = 0;  // west face
	*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) =UnSatPer_face[0]*(*(SurfArea+Index+3*0+1));	
	
	face = 1;  // east face
	*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) =UnSatPer_face[1]*(*(SurfArea+Index+3*1+1));
	
	face = 2;  // south face
	*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) =UnSatPer_face[2]*(*(SurfArea+Index+3*2+1));
 				
 	face = 3;  // north face
 	*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) =UnSatPer_face[3]*(*(SurfArea+Index+3*3+1));
 	
 	face = 4;  // bottom face
 	*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = UnSatPer_face[4]*(*(SurfArea+Index+3*4+1));
 			
 	face = 5;  // top face
 	*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = UnSatPer_face[5]*(*(SurfArea+Index+3*5+1));
	}
	
	return(0);
}


