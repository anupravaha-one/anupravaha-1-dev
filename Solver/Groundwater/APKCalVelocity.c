#include "APKCalVelocity.h"
/**
this function will calculate velocity based on the pressure for groundwater flow.(saturated , unsaturated and variably saturated region).
**/
int CalVelocity(Base *base)
{
	double *H1,*slabPressureHead[6],*u,*v,*KunSat;
	double  k1,D[3];
	Zone *zone;
	int i,j,k,NJ,NK,NI,*ptrSize,numZones,iZone,index1,index[3];
        double *Gradu,*Gradv;
	numZones=base->Number_of_zones;
	if(base->flagVelocity==0)	
	 {
		for(iZone=0;iZone<numZones;iZone++)
		{      	
			zone=base->CgnsZone+iZone;
			MakeSlabPtrArray(zone,slabPressureHead,19);
			ptrSize=zone->Size;
			
			NI=*(ptrSize+0);
			NJ=*(ptrSize+1);
			NK=*(ptrSize+2);
				
			H1 = zone->Flowsol->H;
			
			u = zone->Flowsol->u; 
			v = zone->Flowsol->v;
			
			Gradu=zone->Flowsol->GradU;
			Gradv=zone->Flowsol->GradV;
			
			
			k1 =  *(zone->Flowsol->SatHydConductivity);
			
			if(base->flagSaturation==0)
			{
			if(base->BaseParameters->flagdat==0)// mean along datum plane 
			{	
				for(k=0;k<NK-1;k++)
				for(i=0;i<NI-1;i++)
				for(j=0;j<NJ-1;j++)
				{
				index1=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				index[0]=k,index[1]=j,index[2]=i;
				CalcGradVar(base,zone,H1,slabPressureHead,index,D);
					
				if(fabs(D[0])<1e-07)
				D[0]=0.0;
				else if(fabs(D[1])<1e-07)	               
				D[1]=0.0;
				
	                
				*(u+index1)=-k1*D[0];
				*(v+index1)=-k1*D[1];
				
				*(Gradu+index1)=D[0];
				*(Gradv+index1)=D[1];
				
				}
			}
	      	 	 else
			{		
				for(k=0;k<NK-1;k++)
				for(i=0;i<NI-1;i++)
				for(j=0;j<NJ-1;j++)
				{
				index1=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				index[0]=k,index[1]=j,index[2]=i;
				CalcGradVar(base,zone,H1,slabPressureHead,index,D);
				
				if(fabs(D[0])<1e-07)
				D[0]=0.0;
				else if(fabs(D[1])<1e-07)	               
				D[1]=0.0;
				
				*(u+index1)=-k1*D[0];
				*(v+index1)=-k1*(D[1]+1.0);
				
				*(Gradu+index1)=D[0];
				*(Gradv+index1)=D[1];
					
				}
			
		        } 
	        
		        }
	   	 else 
	   	 {  
	     	
		 KunSat=zone->Flowsol->UnSatPermeability; 
	      
	     	  if(base->BaseParameters->flagdat==0)// mean along datum plane 
		   {	
			for(k=0;k<NK-1;k++)
			for(i=0;i<NI-1;i++)
			for(j=0;j<NJ-1;j++)
			{
			index1=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			CalcGradVar(base,zone,H1,slabPressureHead,index,D);
			
			if(fabs(D[0])<1e-07)
			D[0]=0.0;
			else if(fabs(D[1])<1e-07)	               
			D[1]=0.0;
			
			
			*(u+index1)=-(*(KunSat+index1))*k1*D[0];
			*(v+index1)=-(*(KunSat+index1))*k1*D[1];
			
			*(Gradu+index1)=D[0];
			*(Gradv+index1)=D[1];
				
	                }
	
	        }
	        else
	        {
			for(k=0;k<NK-1;k++)
			for(i=0;i<NI-1;i++)
			for(j=0;j<NJ-1;j++)
			{
			index1=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			CalcGradVar(base,zone,H1,slabPressureHead,index,D);
			
			if(fabs(D[0])<1e-07)
			D[0]=0.0;
			else if(fabs(D[1])<1e-07)	               
			D[1]=0.0;
			
			
			*(Gradu+index1)=D[0];
			*(Gradv+index1)=D[1];
				     
			
			*(u+index1)=-(*(KunSat+index1))*k1*D[0];
			*(v+index1)=-(*(KunSat+index1))*k1*(D[1]+1.0);
	                
			}
		
		}
	
	    
	       }
	
	} 
	}
	
	RefreshVariableCoefficients(1,base);
	RefreshVariableCoefficients(2,base);  
	
	return 0;
	
}

int CalVelocity3d(Base *base)
{
	double *H1,*slabPressureHead[6],*u,*v,*w,*KunSat;
	double  k1,D[3];
	Zone *zone;
	int i,j,k,NJ,NK,NI,*ptrSize,numZones,iZone,index1,index[3];
        double *Gradu,*Gradv,*Gradw;
	numZones=base->Number_of_zones;
	if(base->flagVelocity==0)	
	{	
	for(iZone=0;iZone<numZones;iZone++)
	{      	
	        zone=base->CgnsZone+iZone;
		MakeSlabPtrArray(zone,slabPressureHead,19);
        	ptrSize=zone->Size;
		
		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);
			
		H1 = zone->Flowsol->H;
		
		u = zone->Flowsol->u; 
		v = zone->Flowsol->v;
		w= zone->Flowsol->w;
		
		Gradu=zone->Flowsol->GradU;
	        Gradv=zone->Flowsol->GradV;
		Gradw=zone->Flowsol->GradW;
		
		k1 =  *(zone->Flowsol->SatHydConductivity);
	        
		if(base->flagSaturation==0)
	        {
		     	for(k=0;k<NK-1;k++)
			for(i=0;i<NI-1;i++)
			for(j=0;j<NJ-1;j++)
			{
			index1=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			CalcGradVar(base,zone,H1,slabPressureHead,index,D);
			     
			if(fabs(D[0])<1e-07)
			D[0]=0.0;
			else if(fabs(D[1])<1e-07)	               
			D[1]=0.0;
			else if(fabs(D[2])<1e-07)	               
			D[2]=0.0;
			
			*(Gradu+index1)=D[0];
			*(Gradv+index1)=D[1];
			*(Gradw+index1)=D[2];	
			
			
			*(u+index1)=-k1*D[0];
			*(v+index1)=-k1*(D[1]+1.0);
			*(u+index1)=-k1*D[2];
				
			}
		
		  } 
	    else 
	    {  
	     	
		 KunSat=zone->Flowsol->UnSatPermeability; 
	      
	     		
			for(k=0;k<NK-1;k++)
			for(i=0;i<NI-1;i++)
			for(j=0;j<NJ-1;j++)
			{
			index1=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			CalcGradVar(base,zone,H1,slabPressureHead,index,D);
			
			
			if(fabs(D[0])<1e-07)
			D[0]=0.0;
			else if(fabs(D[1])<1e-07)	               
			D[1]=0.0;
			else if(fabs(D[2])<1e-07)	               
			D[2]=0.0;
			
			
			*(Gradu+index1)=D[0];
			*(Gradv+index1)=D[1];
			*(Gradw+index1)=D[2];	
			
			     
			*(u+index1)=-(*(KunSat+index1))*k1*D[0];
			*(v+index1)=-(*(KunSat+index1))*k1*(D[1]+1.0);
	                *(w+index1)=-(*(KunSat+index1))*k1*D[2];
			}
		      }
	} 
	}
	RefreshVariableCoefficients(1,base);
	RefreshVariableCoefficients(2,base);  
	RefreshVariableCoefficients(3,base);  
	
	return 0;

}
