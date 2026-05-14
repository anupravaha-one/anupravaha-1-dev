#include "APKupdateBCslabsVar_Coeff_PorousMedia.h"

/**
This function is used to calculate the boundary face values of gradient of scalars and dispersion coefficients
and also to update the boundary values of u,v,w. It has function UpdateBcDarcyVelocity and  updateBcDispersionCoeff, to update the disperison coeffcients based on boundary velocity(Darcy velocity) and calculate the velocity based on gradient at boundary face.
*/


int updateBCslabsVar_Coeff_PorousMedia(Base *base,int variableTag,Zone *zone,ZoneBoundary *bcPtr,int dir)
{
        int k,j,i,NI,NJ,NK;
	double *variable,*domain;
	int *ptrRange,index;
	double *coordX,*coordY,*coordZ;
	int c,nv,number;
	double a;
	
	ptrRange = bcPtr->Point_Range;
	NI = *(zone->Size+0);
	NJ = *(zone->Size+1);
	NK = *(zone->Size+2);
	
	coordX=(zone->CgnsGrid)->X;
	coordY=(zone->CgnsGrid)->Y; 
	coordZ=(zone->CgnsGrid)->Z;
	
	variable = FindSlabPtrBC(bcPtr,variableTag);
		
	
	if(variableTag>=60000 && variableTag<=60009)
	{
        	number=variableTag-60000;
		domain=zone->Flowsol->GradScalX[number];
				
		nv = find_nv(base,600+number);
	
		if(((bcPtr)->Data_Set+nv)->Bc_Type==0)
		a=0.0;
		else
		a=1.0;
		
		
		
		switch(dir)
			{
				case 1: 
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(NJ-1)+j;
						i=(NI-2);
						
					*(variable+index) =a*(*(domain+(NJ-1)*(NI-1)*k+(NI-1)*j+i));
			
					}
					//exit(1);
				break;

				case 2: 
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(NJ-1)+j;
						i=0;
						*(variable+index) =a*(*(domain+(NJ-1)*(NI-1)*k+(NI-1)*j+i));
					}
				break;

				case 3:
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = k*(NI-1)+i;
						j=(NJ-2);
						*(variable+index) =a*(*(domain+(NJ-1)*(NI-1)*k+(NI-1)*j+i));
						
					}
				break;

				case 4:
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = k*(NI-1)+i;
						j=0;
						*(variable+index) =a*(*(domain+(NJ-1)*(NI-1)*k+(NI-1)*j+i));		
					}
				break;

				case 5:
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = j*(NI-1)+i;
						k=(NK-2);
						*(variable+index) =a*(*(domain+(NJ-1)*(NI-1)*k+(NI-1)*j+i));
					}
				break;

				case 6:
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = j*(NI-1)+i;
						k=0;
						*(variable+index) =a*(*(domain+(NJ-1)*(NI-1)*k+(NI-1)*j+i));
					}
				break;
				
				default:
					printf("Case not found,Terminating the Program from  updateBCslabsVar_Coeff_PorousMedia\n");
					exit(1);
				break;
			}
	          }	
		else if(variableTag>=60100 && variableTag<=60109)	
		 {	
		     number=variableTag-60100;
		     domain=zone->Flowsol->GradScalY[number];	
		     
		        nv = find_nv(base,600+number);
	
			if(((bcPtr)->Data_Set+nv)->Bc_Type==0)
			a=0.0;
			else
			a=1.0;
		
		     
		     switch(dir)
			{
				case 1: 
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(NJ-1)+j;
						i=(NI-2);
						
					*(variable+index) =a*(*(domain+(NJ-1)*(NI-1)*k+(NI-1)*j+i));
					}
					
				break;

				case 2: 
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(NJ-1)+j;
						i=0;
						*(variable+index) =a*(*(domain+(NJ-1)*(NI-1)*k+(NI-1)*j+i));
					}
				break;

				case 3:
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = k*(NI-1)+i;
						j=(NJ-2);
						*(variable+index) =a*(*(domain+(NJ-1)*(NI-1)*k+(NI-1)*j+i));
					}
				break;

				case 4:
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = k*(NI-1)+i;
						j=0;
						*(variable+index) =a*(*(domain+(NJ-1)*(NI-1)*k+(NI-1)*j+i));	
					}
				break;

				case 5:
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = j*(NI-1)+i;
						k=(NK-2);
						*(variable+index) =a*(*(domain+(NJ-1)*(NI-1)*k+(NI-1)*j+i));
					}
				break;

				case 6:
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = j*(NI-1)+i;
						k=0;
						
						*(variable+index) =a*(*(domain+(NJ-1)*(NI-1)*k+(NI-1)*j+i));
					}
				break;
				
				default:
					printf("Case not found,Terminating the Program from  updateBCslabsVar_Coeff_PorousMedia\n");
					exit(1);
				break;
			}
		}	
		else if(variableTag>=60200 && variableTag<=60209)	
		{	
			number=variableTag-60200;
		  	domain=zone->Flowsol->GradScalZ[number];	
		
		        nv = find_nv(base,600+number);
	  
			if(((bcPtr)->Data_Set+nv)->Bc_Type==0)
			a=0.0;
			else
			a=1.0;
		     
		     		
		switch(dir)
			{
				case 1: 
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(NJ-1)+j;
						i=(NI-2);
						
					*(variable+index) =a*(*(domain+(NJ-1)*(NI-1)*k+(NI-1)*j+i));
					}
					
				break;

				case 2: 
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(NJ-1)+j;
						i=0;
						*(variable+index) =a*(*(domain+(NJ-1)*(NI-1)*k+(NI-1)*j+i));
			//printf("\n var=%lf   \tdomain=%lf",*(variable+index),*(domain+(NJ-1)*(NI-1)*k+(NI-1)*j+i));
					}
				break;

				case 3:
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = k*(NI-1)+i;
						j=(NJ-2);
						*(variable+index) =a*(*(domain+(NJ-1)*(NI-1)*k+(NI-1)*j+i));
						
					}
				break;

				case 4:
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = k*(NI-1)+i;
						j=0;
						*(variable+index) =a*(*(domain+(NJ-1)*(NI-1)*k+(NI-1)*j+i));
		//printf("\n var=%lf  \tdomain=%lf",*(variable+index),*(domain+(NJ-1)*(NI-1)*k+(NI-1)*j+i));		
					}
				break;

				case 5:
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = j*(NI-1)+i;
						k=(NK-2);
						*(variable+index) =a*(*(domain+(NJ-1)*(NI-1)*k+(NI-1)*j+i));
					}
				break;

				case 6:
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = j*(NI-1)+i;
						k=0;
						
					*(variable+index) =a*(*(domain+(NJ-1)*(NI-1)*k+(NI-1)*j+i));
					}
				break;
				
				default:
					printf("Case not found,Terminating the Program from  updateBCslabsVar_Coeff_PorousMedia\n");
					exit(1);
				break;
			}
		}	
		else if(variableTag>=700 && variableTag<=709)	
		{  	   
		   number=variableTag-700;
		   domain=zone->Flowsol->Scalcoeff[number];	
		   
		   if(zone->CalDispersion==0)
	           { 
	          switch(dir)
			{
				case 1: 
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(NJ-1)+j;
						i=(NI-2);
					
						for(c=0;c<9;c++)
					{*(variable+9*index+c) =*(domain+9*((NJ-1)*(NI-1)*k+(NI-1)*j+i)+c);
			//printf("\n var=%lf   \tdomain=%lf",*(variable+index),*(domain+(NJ-1)*(NI-1)*k+(NI-1)*j+i));
					}
					}
					
				break;

				case 2: 
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(NJ-1)+j;
						i=0;
					for(c=0;c<9;c++)
					{*(variable+9*index+c) =*(domain+9*((NJ-1)*(NI-1)*k+(NI-1)*j+i)+c);
			//printf("\n var=%lf   \tdomain=%lf",*(variable+index),*(domain+(NJ-1)*(NI-1)*k+(NI-1)*j+i));
					}
					}
				break;

				case 3:
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = k*(NI-1)+i;
						j=(NJ-2);
					for(c=0;c<9;c++)
					{*(variable+9*index+c) =*(domain+9*((NJ-1)*(NI-1)*k+(NI-1)*j+i)+c);
					//printf("  %lf",*(variable+9*index+c));
					}
					//printf("\n ");
					}
				break;

				case 4:
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = k*(NI-1)+i;
						j=0;
					for(c=0;c<9;c++)
					*(variable+9*index+c) =*(domain+9*((NJ-1)*(NI-1)*k+(NI-1)*j+i)+c);
					
					
						
					}
				break;

				case 5:
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = j*(NI-1)+i;
						k=(NK-2);
					for(c=0;c<9;c++)
					*(variable+9*index+c) =*(domain+9*((NJ-1)*(NI-1)*k+(NI-1)*j+i)+c);
					}
				break;

				case 6:
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = j*(NI-1)+i;
						k=0;
						for(c=0;c<9;c++)
					*(variable+9*index+c) =*(domain+9*((NJ-1)*(NI-1)*k+(NI-1)*j+i)+c);
					}
				break;
				
				default:
					printf("Case not found,terminating from updateBCslabsVar_Coeff_PorousMedia\n");
					exit(1);
				break;
			}
		}
		  else
		   {
		   updateBcDispersionCoeff(base,variableTag,zone,bcPtr,dir);
		   }
		   }	
		  else if(variableTag>=1&&variableTag<=3)
		   {	
		   	   
		        if(base->flagVelocity==0)	
		        {
			
			UpdateBcDarcyVelocity(base,variableTag,zone,bcPtr, dir);
		     	
			
			}
			else 
		        {
			if(variableTag==1)
			domain=zone->Flowsol->u;
			else if(variableTag==2)
			domain=zone->Flowsol->v;
			else if(variableTag==3)
			domain=zone->Flowsol->w;
			
			switch(dir)
			{
				case 1: 
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(NJ-1)+j;
						i=(NI-2);
						
					*(variable+index) =*(domain+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
			
					}
					//exit(1);
				break;

				case 2: 
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(NJ-1)+j;
						i=0;
						*(variable+index) =*(domain+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
					}
				break;

				case 3:
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = k*(NI-1)+i;
						j=(NJ-2);
						*(variable+index) =*(domain+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
						
					}
				break;

				case 4:
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = k*(NI-1)+i;
						j=0;
						*(variable+index) =*(domain+(NJ-1)*(NI-1)*k+(NI-1)*j+i);		
					}
				break;

				case 5:
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = j*(NI-1)+i;
						k=(NK-2);
						*(variable+index) =*(domain+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
					}
				break;

				case 6:
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = j*(NI-1)+i;
						k=0;
						*(variable+index) =*(domain+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
					}
				break;
				
				default:
					printf("Case not found,Terminating the Program from  updateBCslabsVar_Coeff_PorousMedia\n");
					exit(1);
				break;
		
		            }
		
		     }
		 }
		
		return(0);
		}
				
	
	
	int UpdateBcDarcyVelocity(Base *base,int variableTag,Zone *zone,ZoneBoundary *bcPtr,int dir)
	{
		
	int k,j,i,NI,NJ,NK;
	
	int *ptrRange,index;
	
	int nv;
	double a,K,c,P;
	double *U,*V,*W;
	double *Gradu,*Gradv,*Gradw;// these store gardient of pressure heaad for x, y, z direction i.e in u, v and w 
	double *Kunsat;
		
	
	ptrRange = bcPtr->Point_Range;
	
	NI = *(zone->Size+0);
	NJ = *(zone->Size+1);
	NK = *(zone->Size+2);
	
	        
		if(base->BaseParameters->flag2d == 1) 
		{
		if(base->flagSaturation==0)
	        {
		        if(base->BaseParameters->flagdat==0)// mean along datum plane 
		      	c=0.0;
			else
			c=1.0;
		      
		        Gradu=zone->Flowsol->GradU;
	       	        Gradv=zone->Flowsol->GradV;
	                		
			nv = find_nv(base,19);
	  		
			if(((bcPtr)->Data_Set+nv)->Bc_Type==0)
			a=0.0;
			else
			a=1.0;
		       
		      
		         U =  bcPtr->ZSlab->u;
		         V=   bcPtr->ZSlab->v;
		         K =  *(zone->Flowsol->SatHydConductivity);  
		       
		    switch(dir)      
		     {  
		              case 1: 
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(NJ-1)+j;
						i=(NI-2);
					
					
					*(U+index) =-a*K*(*(Gradu+(NJ-1)*(NI-1)*k+(NI-1)*j+i));
					*(V+index) =-a*K*(*(Gradv+(NJ-1)*(NI-1)*k+(NI-1)*j+i)+c);			
					}
		       
		                         break;

				case 2: 
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(NJ-1)+j;
						i=0;
					*(U+index) =-a*K*(*(Gradu+(NJ-1)*(NI-1)*k+(NI-1)*j+i));
					*(V+index) =-a*K*(*(Gradv+(NJ-1)*(NI-1)*k+(NI-1)*j+i)+c);	
					}
					
				break;

				case 3:
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = k*(NI-1)+i;
						j=(NJ-2);
					*(U+index) =-a*K*(*(Gradu+(NJ-1)*(NI-1)*k+(NI-1)*j+i));
					*(V+index) =-a*K*(*(Gradv+(NJ-1)*(NI-1)*k+(NI-1)*j+i)+c);	
					//printf("\n ");
					}
				break;

				case 4:
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = k*(NI-1)+i;
						j=0;
					*(U+index) =-a*K*(*(Gradu+(NJ-1)*(NI-1)*k+(NI-1)*j+i));
					*(V+index) =-a*K*(*(Gradv+(NJ-1)*(NI-1)*k+(NI-1)*j+i)+c);	
					}
				break;

				case 5:
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = j*(NI-1)+i;
						k=(NK-2);
					*(U+index) =-a*K*(*(Gradu+(NJ-1)*(NI-1)*k+(NI-1)*j+i));
					*(V+index) =-a*K*(*(Gradv+(NJ-1)*(NI-1)*k+(NI-1)*j+i)+c);	
					}
				break;

				case 6:
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = j*(NI-1)+i;
						k=0;
					*(U+index) =-a*K*(*(Gradu+(NJ-1)*(NI-1)*k+(NI-1)*j+i));
					*(V+index) =-a*K*(*(Gradv+(NJ-1)*(NI-1)*k+(NI-1)*j+i)+c);	
					}
				break;
				
				default:
					printf("Case not found,Terminating the Program from UpdateBcDarcyVelocity\n");
					exit(1);
				break;
	               
		      }    
		       
		       
		        
	               
		       
		}
		else
		{
		
			 if(base->BaseParameters->flagdat==0)// mean along datum plane 
		      	c=0.0;
			else
			c=1.0;
		      
		        Gradu=zone->Flowsol->GradU;
	       	        Gradv=zone->Flowsol->GradV;
	                		
			nv = find_nv(base,19);
	  		
			if(((bcPtr)->Data_Set+nv)->Bc_Type==0)
			a=0.0;
			else
			a=1.0;
		       
		   
			
			  Kunsat=bcPtr->ZSlab->UnSatPermeability; 
		         U =  bcPtr->ZSlab->u;
		         V=   bcPtr->ZSlab->v;
		         K =  *(zone->Flowsol->SatHydConductivity);  
		       
			
		    
		    switch(dir)      
		     {  
		              case 1: 
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(NJ-1)+j;
						i=(NI-2);
					
					   P  = *(Kunsat+index);
					*(U+index) =-a*K*(*(Gradu+(NJ-1)*(NI-1)*k+(NI-1)*j+i))*P;
					*(V+index) =-a*K*(*(Gradv+(NJ-1)*(NI-1)*k+(NI-1)*j+i)+c)*P;	
					}
		       
		                         break;

				case 2: 
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(NJ-1)+j;
						i=0;
					 P  = *(Kunsat+index);
					*(U+index) =-a*K*(*(Gradu+(NJ-1)*(NI-1)*k+(NI-1)*j+i))*P;
					*(V+index) =-a*K*(*(Gradv+(NJ-1)*(NI-1)*k+(NI-1)*j+i)+c)*P;	
					}
					
				break;

				case 3:
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = k*(NI-1)+i;
						j=(NJ-2);
					 P  = *(Kunsat+index);
					*(U+index) =-a*K*(*(Gradu+(NJ-1)*(NI-1)*k+(NI-1)*j+i))*P;
					*(V+index) =-a*K*(*(Gradv+(NJ-1)*(NI-1)*k+(NI-1)*j+i)+c)*P;		
					
					}
				break;

				case 4:
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = k*(NI-1)+i;
						j=0;
					 P  = *(Kunsat+index);
					*(U+index) =-a*K*(*(Gradu+(NJ-1)*(NI-1)*k+(NI-1)*j+i))*P;
					*(V+index) =-a*K*(*(Gradv+(NJ-1)*(NI-1)*k+(NI-1)*j+i)+c)*P;		
					}
				break;

				case 5:
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = j*(NI-1)+i;
						k=(NK-2);
					 P  = *(Kunsat+index);
					*(U+index) =-a*K*(*(Gradu+(NJ-1)*(NI-1)*k+(NI-1)*j+i))*P;
					*(V+index) =-a*K*(*(Gradv+(NJ-1)*(NI-1)*k+(NI-1)*j+i)+c)*P;		
					}
				break;

				case 6:
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = j*(NI-1)+i;
						k=0;
					 P  = *(Kunsat+index);
					*(U+index) =-a*K*(*(Gradu+(NJ-1)*(NI-1)*k+(NI-1)*j+i))*P;
					*(V+index) =-a*K*(*(Gradv+(NJ-1)*(NI-1)*k+(NI-1)*j+i)+c)*P;		
					}
				break;
				
				default:
					printf("Case not found,Terminating the Program from UpdateBcDarcyVelocity\n");
					exit(1);
				break;
	               
		      }    
		       
		
		
		
		}
		}
		else
		{
		if(base->flagSaturation==0)
	        {
		        Gradu=zone->Flowsol->GradU;
	       	        Gradv=zone->Flowsol->GradV;
	                Gradw=zone->Flowsol->GradW;		
			
			nv = find_nv(base,19);
	  		
			if(((bcPtr)->Data_Set+nv)->Bc_Type==0)
			a=0.0;
			else
			a=1.0;
		       
		      
		         U =  bcPtr->ZSlab->u;
		         V=   bcPtr->ZSlab->v;
			 W=   bcPtr->ZSlab->w;
			 
		         K =  *(zone->Flowsol->SatHydConductivity);  
		       
		    switch(dir)      
		     {  
		              case 1: 
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(NJ-1)+j;
						i=(NI-2);
					
					
					*(U+index) =-a*K*(*(Gradu+(NJ-1)*(NI-1)*k+(NI-1)*j+i));
					*(V+index) =-a*K*(*(Gradv+(NJ-1)*(NI-1)*k+(NI-1)*j+i)+1.0);		*(W+index) =-a*K*(*(Gradw+(NJ-1)*(NI-1)*k+(NI-1)*j+i));	
					}
		       
		                         break;

				case 2: 
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(NJ-1)+j;
						i=0;
					*(U+index) =-a*K*(*(Gradu+(NJ-1)*(NI-1)*k+(NI-1)*j+i));
					*(V+index) =-a*K*(*(Gradv+(NJ-1)*(NI-1)*k+(NI-1)*j+i)+1.0);		*(W+index) =-a*K*(*(Gradw+(NJ-1)*(NI-1)*k+(NI-1)*j+i));	
					}
					
				break;

				case 3:
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = k*(NI-1)+i;
						j=(NJ-2);
					*(U+index) =-a*K*(*(Gradu+(NJ-1)*(NI-1)*k+(NI-1)*j+i));
					*(V+index) =-a*K*(*(Gradv+(NJ-1)*(NI-1)*k+(NI-1)*j+i)+1.0);		*(W+index) =-a*K*(*(Gradw+(NJ-1)*(NI-1)*k+(NI-1)*j+i));	
					}
				break;

				case 4:
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = k*(NI-1)+i;
						j=0;
					*(U+index) =-a*K*(*(Gradu+(NJ-1)*(NI-1)*k+(NI-1)*j+i));
					*(V+index) =-a*K*(*(Gradv+(NJ-1)*(NI-1)*k+(NI-1)*j+i)+1.0);		*(W+index) =-a*K*(*(Gradw+(NJ-1)*(NI-1)*k+(NI-1)*j+i));	
					}
				break;

				case 5:
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = j*(NI-1)+i;
						k=(NK-2);
					*(U+index) =-a*K*(*(Gradu+(NJ-1)*(NI-1)*k+(NI-1)*j+i));
					*(V+index) =-a*K*(*(Gradv+(NJ-1)*(NI-1)*k+(NI-1)*j+i)+1.0);		*(W+index) =-a*K*(*(Gradw+(NJ-1)*(NI-1)*k+(NI-1)*j+i));	
					}
				break;

				case 6:
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = j*(NI-1)+i;
						k=0;
					*(U+index) =-a*K*(*(Gradu+(NJ-1)*(NI-1)*k+(NI-1)*j+i));
					*(V+index) =-a*K*(*(Gradv+(NJ-1)*(NI-1)*k+(NI-1)*j+i)+1.0);		*(W+index) =-a*K*(*(Gradw+(NJ-1)*(NI-1)*k+(NI-1)*j+i));		
					}
				break;
				
				default:
					printf("Case not found,Terminating the Program from UpdateBcDarcyVelocity\n");
					exit(1);
				break;
	               
		      }    
		       
		       
		        
	               
		       
		}
		else
		{
		
		 
		        Gradu=zone->Flowsol->GradU;
	       	        Gradv=zone->Flowsol->GradV;
	                Gradw=zone->Flowsol->GradW; 		
			
			
			nv = find_nv(base,19);
	  		
			if(((bcPtr)->Data_Set+nv)->Bc_Type==0)
			a=0.0;
			else
			a=1.0;
		       
		      
		         U =  bcPtr->ZSlab->u;
		         V=   bcPtr->ZSlab->v;
			 W=   bcPtr->ZSlab->w;
			 
		         K =  *(zone->Flowsol->SatHydConductivity);  
		      
			 
			  Kunsat=bcPtr->ZSlab->UnSatPermeability; 
		    switch(dir)      
		     {  
		              case 1: 
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(NJ-1)+j;
						i=(NI-2);
					
					   P  = *(Kunsat+index);
					*(U+index) =-a*K*(*(Gradu+(NJ-1)*(NI-1)*k+(NI-1)*j+i))*P;
					*(V+index) =-a*K*(*(Gradv+(NJ-1)*(NI-1)*k+(NI-1)*j+i)+1.0)*P;		*(W+index) =-a*K*(*(Gradw+(NJ-1)*(NI-1)*k+(NI-1)*j+i))*P;		
					}
		       
		                         break;

				case 2: 
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(NJ-1)+j;
						i=0;
						  P  = *(Kunsat+index);
					*(U+index) =-a*K*(*(Gradu+(NJ-1)*(NI-1)*k+(NI-1)*j+i))*P;
					*(V+index) =-a*K*(*(Gradv+(NJ-1)*(NI-1)*k+(NI-1)*j+i)+1.0)*P;		*(W+index) =-a*K*(*(Gradw+(NJ-1)*(NI-1)*k+(NI-1)*j+i))*P;	
					}
					
				break;

				case 3:
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = k*(NI-1)+i;
						j=(NJ-2);
						  P  = *(Kunsat+index);
					*(U+index) =-a*K*(*(Gradu+(NJ-1)*(NI-1)*k+(NI-1)*j+i))*P;
					*(V+index) =-a*K*(*(Gradv+(NJ-1)*(NI-1)*k+(NI-1)*j+i)+1.0)*P;		*(W+index) =-a*K*(*(Gradw+(NJ-1)*(NI-1)*k+(NI-1)*j+i))*P;			
					//printf("\n ");
					}
				break;

				case 4:
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = k*(NI-1)+i;
						j=0;
						  P  = *(Kunsat+index);
					*(U+index) =-a*K*(*(Gradu+(NJ-1)*(NI-1)*k+(NI-1)*j+i))*P;
					*(V+index) =-a*K*(*(Gradv+(NJ-1)*(NI-1)*k+(NI-1)*j+i)+1.0)*P;		*(W+index) =-a*K*(*(Gradw+(NJ-1)*(NI-1)*k+(NI-1)*j+i))*P;		
					}
				break;

				case 5:
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = j*(NI-1)+i;
						k=(NK-2);
					 P  = *(Kunsat+index);
					*(U+index) =-a*K*(*(Gradu+(NJ-1)*(NI-1)*k+(NI-1)*j+i))*P;
					*(V+index) =-a*K*(*(Gradv+(NJ-1)*(NI-1)*k+(NI-1)*j+i)+c)*P;
					*(W+index) =-a*K*(*(Gradw+(NJ-1)*(NI-1)*k+(NI-1)*j+i))*P;
							
					}
				break;

				case 6:
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = j*(NI-1)+i;
						k=0;
					 P  = *(Kunsat+index);
					*(U+index) =-a*K*(*(Gradu+(NJ-1)*(NI-1)*k+(NI-1)*j+i))*P;
					*(V+index) =-a*K*(*(Gradv+(NJ-1)*(NI-1)*k+(NI-1)*j+i)+c)*P;	
					*(W+index) =-a*K*(*(Gradw+(NJ-1)*(NI-1)*k+(NI-1)*j+i))*P;
						
					}
				break;
				
				default:
					printf("Case not found,Terminating the Program from UpdateBcDarcyVelocity\n");
					exit(1);
				break;
	               
		      }    
		       
		
		
		
		}
		
		
		
		}
	return(0);
	}
		
	int  updateBcDispersionCoeff(Base *base,int variableTag,Zone *zone,ZoneBoundary *bcPtr,int dir)
	{
	 
	int k,j,i,NI,NJ,NK;
	
	int *ptrRange,index;
	//double surfPoints[4][3],centNeighbour[3];
	//int nv,number;
	double L,T;
	double *u1,*v1,*w1;
	double U,u,v,w;
	double *Disp;
	 
	 Disp = FindSlabPtrBC(bcPtr,variableTag);
	  
	 ptrRange = bcPtr->Point_Range;
	
	NI = *(zone->Size+0);
	NJ = *(zone->Size+1);
	NK = *(zone->Size+2);
	      
	      
	              
	L=zone->alphaL;
	T=zone->alphaT;
	      	      
	      
	
		   
		if(base->BaseParameters->flag2d == 1) 
		{     
			
		u1 =  bcPtr->ZSlab->u;
		v1=   bcPtr->ZSlab->v;
	           
		
		   switch(dir)
			{
				case 1: 
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(NJ-1)+j;
						i=(NI-2);
					
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
	      				  *( Disp+9*index+8)= 0.0;      //Dzz
	       
		 
								
					}
					
				break;

				case 2: 
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(NJ-1)+j;
						i=0;
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
	      				  *( Disp+9*index+8)= 0.0;      //Dzz
					}
				break;

				case 3:
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = k*(NI-1)+i;
						j=(NJ-2);
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
	      				  *( Disp+9*index+8)= 0.0;      //Dzz
					}
				break;

				case 4:
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = k*(NI-1)+i;
						j=0;
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
	      				  *( Disp+9*index+8)= 0.0;      //Dzz
					
					
						
					}
				break;

				case 5:
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = j*(NI-1)+i;
						k=(NK-2);
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
	      				  *( Disp+9*index+8)= 0.0;      //Dzz
					}
				break;

				case 6:
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = j*(NI-1)+i;
						k=0;
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
	      				  *( Disp+9*index+8)= 0.0;      //Dzz
					}
				break;
				
				default:
					printf("Case not found,Terminating the Program from updateBcDispersionCoeff\n");
					exit(1);
				break;
			}
			}
			else
			{
			
			u1 = zone->Flowsol->u; 
	                v1 = zone->Flowsol->v;
	                w1 = zone->Flowsol->w;
			
			 switch(dir)
			{
				case 1: 
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(NJ-1)+j;
						i=(NI-2);
					
					 u=*(u1+index);      
	   				 v=*(v1+index);      
	                                 w=*(w1+index); 
	     				 U=sqrt((u*u+v*v+w*w));
	      	      
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
					
				break;

				case 2: 
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(NJ-1)+j;
						i=0;
					 u=*(u1+index);      
	   				 v=*(v1+index);      
	                                 w=*(w1+index); 
	     				 U=sqrt((u*u+v*v+w*w));
	      	      
	    				*( Disp+9*index+0)=   T*U+(L-T)*u*u/U;  //Dxx
	        			*( Disp+9*index+1)=  (L-T)*u*v/U;      //Dxy
	        			*( Disp+9*index+2)=  (L-T)*u*w/U;              //Dxz
	      
	        			*( Disp+9*index+3)=  *( Disp+9*index+1);     //Dxy
	        			*( Disp+9*index+4)=   T*U+(L-T)*v*v/U;   //Dya
	       				 *( Disp+9*index+5)=  (L-T)*v*w/U;        //Dyz
	      
	     				*( Disp+9*index+6)=   *( Disp+9*index+2);     //Dzx
	      				*( Disp+9*index+7)=   *( Disp+9*index+5);   //Dzy
	        			*( Disp+9*index+8)=   T*U+(L-T)*w*w/U;    
					}
				break;

				case 3:
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = k*(NI-1)+i;
						j=(NJ-2);
					 u=*(u1+index);      
	   				 v=*(v1+index);      
	                                 w=*(w1+index); 
	     				 U=sqrt((u*u+v*v+w*w));
	      	      
	    				*( Disp+9*index+0)=   T*U+(L-T)*u*u/U;  //Dxx
	        			*( Disp+9*index+1)=  (L-T)*u*v/U;      //Dxy
	        			*( Disp+9*index+2)=  (L-T)*u*w/U;              //Dxz
	      
	        			*( Disp+9*index+3)=  *( Disp+9*index+1);     //Dxy
	        			*( Disp+9*index+4)=   T*U+(L-T)*v*v/U;   //Dya
	       				 *( Disp+9*index+5)=  (L-T)*v*w/U;        //Dyz
	      
	     				*( Disp+9*index+6)=   *( Disp+9*index+2);     //Dzx
	      				*( Disp+9*index+7)=   *( Disp+9*index+5);   //Dzy
	        			*( Disp+9*index+8)=   T*U+(L-T)*w*w/U;    
					}
				break;

				case 4:
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = k*(NI-1)+i;
						j=0;
					  u=*(u1+index);      
	   				 v=*(v1+index);      
	                                 w=*(w1+index); 
	     				 U=sqrt((u*u+v*v+w*w));
	      	      
	    				*( Disp+9*index+0)=   T*U+(L-T)*u*u/U;  //Dxx
	        			*( Disp+9*index+1)=  (L-T)*u*v/U;      //Dxy
	        			*( Disp+9*index+2)=  (L-T)*u*w/U;              //Dxz
	      
	        			*( Disp+9*index+3)=  *( Disp+9*index+1);     //Dxy
	        			*( Disp+9*index+4)=   T*U+(L-T)*v*v/U;   //Dya
	       				*( Disp+9*index+5)=  (L-T)*v*w/U;        //Dyz
	      
	     				*( Disp+9*index+6)=   *( Disp+9*index+2);     //Dzx
	      				*( Disp+9*index+7)=   *( Disp+9*index+5);   //Dzy
	        			*( Disp+9*index+8)=   T*U+(L-T)*w*w/U;    
					
					
						
					}
				break;

				case 5:
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = j*(NI-1)+i;
						k=(NK-2);
					 u=*(u1+index);      
	   				 v=*(v1+index);      
	                                 w=*(w1+index); 
	     				 U=sqrt((u*u+v*v+w*w));
	      	      
	    				*( Disp+9*index+0)=   T*U+(L-T)*u*u/U;  //Dxx
	        			*( Disp+9*index+1)=  (L-T)*u*v/U;      //Dxy
	        			*( Disp+9*index+2)=  (L-T)*u*w/U;              //Dxz
	      
	        			*( Disp+9*index+3)=  *( Disp+9*index+1);     //Dxy
	        			*( Disp+9*index+4)=   T*U+(L-T)*v*v/U;   //Dya
	       				 *( Disp+9*index+5)=  (L-T)*v*w/U;        //Dyz
	      
	     				*( Disp+9*index+6)=   *( Disp+9*index+2);     //Dzx
	      				*( Disp+9*index+7)=   *( Disp+9*index+5);   //Dzy
	        			*( Disp+9*index+8)=   T*U+(L-T)*w*w/U;    
					}
				break;

				case 6:
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = j*(NI-1)+i;
						k=0;
					 u=*(u1+index);      
	   				 v=*(v1+index);      
	                                 w=*(w1+index); 
	     				 U=sqrt((u*u+v*v+w*w));
	      	      
	    				*( Disp+9*index+0)=   T*U+(L-T)*u*u/U;  //Dxx
	        			*( Disp+9*index+1)=  (L-T)*u*v/U;      //Dxy
	        			*( Disp+9*index+2)=  (L-T)*u*w/U;              //Dxz
	      
	        			*( Disp+9*index+3)=  *( Disp+9*index+1);     //Dxy
	        			*( Disp+9*index+4)=   T*U+(L-T)*v*v/U;   //Dya
	       				 *( Disp+9*index+5)=  (L-T)*v*w/U;        //Dyz
	      
	     				*( Disp+9*index+6)=   *( Disp+9*index+2);     //Dzx
	      				*( Disp+9*index+7)=   *( Disp+9*index+5);   //Dzy
	        			*( Disp+9*index+8)=   T*U+(L-T)*w*w/U;    
					}
				break;
				
				default:
					printf("Case not found,Terminating the Program from updateBcDispersionCoeff\n");
					exit(1);
				break;
			}
			
		}
	return(0);
}


