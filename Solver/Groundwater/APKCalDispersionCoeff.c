#include "APKCalDispersionCoeff.h"

/**
This routine which is called in the initialize function performs the task of creating and alloting the sufficient
memory for that particular variable in Zone 'slab' structure for both Boundary and Interface.
*/
int CalDispersionCoeff2d(Base *base,int cntvariable)
{        
	int k,i,j,index,numZones,iZone,*ptrSize,NI,NJ,NK;
	double *u1,*v1;
	Zone *zone;
	double U,*Disp;
	double u,v,L,T;
	
	numZones=(base->Number_of_zones);
	
	for(iZone=0;iZone<numZones;iZone++)
	{      
	      zone=base->CgnsZone+iZone;
	      Disp=((zone)->Flowsol)->Scalcoeff[cntvariable];
	      
	      u1 = zone->Flowsol->u; 
	      v1 = zone->Flowsol->v;
	      if(zone->CalDispersion==0)
	      { 
	                continue; 
	       
	      }
	     else
	      {
	            
				L=zone->alphaL;
				T=zone->alphaT;
					
				ptrSize=zone->Size;
				
				NI=*(ptrSize+0);
				NJ=*(ptrSize+1);
				NK=*(ptrSize+2);
					
					
				for(k=0;k<NK-1;k++)
				for(i=0;i<NI-1;i++)
				for(j=0;j<NJ-1;j++)
				{ 
				index=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				
				u=*(u1+index);      
				v=*(v1+index);      
				
				U=sqrt((u*u+v*v));
					
				*( Disp+9*index+0)=  T*U+(L-T)*u*u/U;  //Dxx
				*( Disp+9*index+1)=  (L-T)*u*v/U;      //Dxy
				*( Disp+9*index+2)=  0.0;              //Dxz
				
				*( Disp+9*index+3)=  *( Disp+9*index+1);     //Dxy
				*( Disp+9*index+4)=  T*U+(L-T)*v*v/U;   //Dya
				*( Disp+9*index+5)=  0.0;      
				
				*( Disp+9*index+6)= 0.0;     //Dzx
				*( Disp+9*index+7)= 0.0;   //Dzy
				*( Disp+9*index+8)=0.0;      //Dzz
				
						
				}
			
		}
	
	}
	
	
	return(0);
}



int CalDispersionCoeff3d(Base *base,int cntvariable)
{
	int k,i,j,index,numZones,iZone,*ptrSize,NI,NJ,NK;
	double *u1,*v1,*w1;
	Zone *zone;
	double U,*Disp;
	double u,v,w,L,T;
	
	numZones=(base->Number_of_zones);
	
	
	for(iZone=0;iZone<numZones;iZone++)
	{      
		zone=base->CgnsZone+iZone;
		Disp=((zone)->Flowsol)->Scalcoeff[cntvariable];
			
		u1 = zone->Flowsol->u; 
		v1 = zone->Flowsol->v;
		w1 = zone->Flowsol->w;
			
		if(zone->CalDispersion==0)
		{    
			continue; 
		}
		else
		{
			L=zone->alphaL;
			T=zone->alphaT;
			
			ptrSize=zone->Size;
			NI=*(ptrSize+0);
			NJ=*(ptrSize+1);
			NK=*(ptrSize+2);
					
			for(k=0;k<NK-1;k++)
			for(i=0;i<NI-1;i++)
			for(j=0;j<NJ-1;j++)
			{ 
				index=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				
				u=*(u1+index);      
				v=*(v1+index);     
				
				U=sqrt((u*u+v*v+w*w));//avaerage velcoity
				
				*( Disp+9*index+0)=   T*U+(L-T)*u*u/U;  //Dxx
				*( Disp+9*index+1)=  (L-T)*u*v/U;      //Dxy
				*( Disp+9*index+2)=  (L-T)*u*w/U;              //Dxz
				
				*( Disp+9*index+3)=  *( Disp+9*index+1);     //Dxy
				*( Disp+9*index+4)=   T*U+(L-T)*v*v/U;   //Dya
				*( Disp+9*index+5)=  (L-T)*v*w/U;        //Dyz
				
				*( Disp+9*index+6)=   *( Disp+9*index+2);     //Dzx
				*( Disp+9*index+7)=   *( Disp+9*index+5);   //Dzy
				*( Disp+9*index+8)=   T*U+(L-T)*w*w/U;     //Dzz
				
			}
		}	
		
	}
	
	return(0);
}
