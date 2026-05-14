#include "APKInputRadiationBoundarycondition.h"

/**
The function InputRadiationBoundarycondition() interacts with user to store the desired boundary condition on a 
particular BC patch. This function is called inside a zone loop for each zone present in the computational 
domain.

The input parameters to this function are, count of total no of variables, pointer of the zone, flag for 2D 
solver. 

It loops over all the BC patches present in the zone. To give an idea of the BC patch loacation BC tag, 
point range of the Bc patch and corrosponding vertex coordinates are printed for this purpose point_coord() 
is called. 
 
Checks one by one the entire range of possible varibles, offered by this solver 
if they are initialized by the user (ie. user wants to solve them)
 
If a particular variable is initialize, for that varible BC options (Symmetric, Adiabatic or Isothermal or Isoflux) are printed, the users responce along with variable name is then stored in the data structue under 
the node Data_set.
 
If user wants to enter a function Call_user_func() is called. 
*/

int InputRadiationBoundarycondition(Base *base, int NS, int RadVar, Zone *zone, int flag)
	{
	char comp1[]="Kmin",comp2[]="Kmax",*varname;
	int BC_index, number, maxn, ni, nj, nk, c, c1, c2, count, s, choice;
	double value,emissivity;
	
	ni=*(zone->Size+0);
	nj=*(zone->Size+1);
	nk=*(zone->Size+2);
	
	printf("\n Number of Radiation Variables =%d", RadVar);
	printf("\n This is zone %s", zone->zonename);
	
	for(BC_index=0;BC_index<(zone->Number_of_BC_Patches);BC_index++)
		{
		if(base->BaseParameters->flagRadiation==1)
			{
		 	count=0;
		 	(zone->ZoneBc+BC_index)->Data_Set=(Dataset *)malloc(RadVar*sizeof(Dataset));
			}
		else
		 	count=NS;
		 c1=strcmp((zone->ZoneBc+BC_index)->BCName,comp1);
		 c2=strcmp((zone->ZoneBc+BC_index)->BCName,comp2);
			//Determine the span of BC Patches; 
		 point_coordRadiation(zone,(zone->ZoneBc+BC_index)->Point_Range);
			 
		 for(number=0;number<RadVar;number++)
		 {
			 count=count+1;
			 varname=find_name(base,count);
			
			if((c1==0||c2==0)&&(flag==1))
				{
			 	c=0;
			  	((zone->ZoneBc+BC_index)->Data_Set+(count-1))->funct=NULL;
				((zone->ZoneBc+BC_index)->Data_Set+(count-1))->funct_beta=NULL;
				((zone->ZoneBc+BC_index)->Data_Set+(count-1))->funct_gama=NULL;
				*(((zone->ZoneBc+BC_index)->Data_Set+(count-1))->identity)=0;
				*(((zone->ZoneBc+BC_index)->Data_Set+(count-1))->identity+1)=0;
				*(((zone->ZoneBc+BC_index)->Data_Set+(count-1))->identity+2)=0;
				*(((zone->ZoneBc+BC_index)->Data_Set+(count-1))->identity+3)=0;
			  	((zone->ZoneBc+BC_index)->Data_Set+(count-1))->Bc_Type=c;
			 	}
				 
			 else
			 {
			 printf("\nI am solving for %s",varname);
			
			 ((zone->ZoneBc+BC_index)->Data_Set+(count-1))->Name_of_variable[0]=varname[0];
			 
			  printf("\nEnter the Boundary Condition for Radiative Heat Transfer\n");
			  printf("\nFor Inlet and Outlet Please Select Iso Boudary Condition With Wall Emissivity 1\n");
			  printf("Enter 0->Symmetric Boundary Condition 1 -> Adiabatic Boundary Condition, 2->Isothermal or IsoFlux Boundary Condition\n");
			  scanf("%d",&c);
			  ((zone->ZoneBc+BC_index)->Data_Set+(count-1))->Bc_Type=c;
			  if(c==0)
			  	{
			  	((zone->ZoneBc+BC_index)->Data_Set+(count-1))->funct=NULL;
				((zone->ZoneBc+BC_index)->Data_Set+(count-1))->funct_beta=NULL;
				((zone->ZoneBc+BC_index)->Data_Set+(count-1))->funct_gama=NULL;
				*(((zone->ZoneBc+BC_index)->Data_Set+(count-1))->identity)=0;
				*(((zone->ZoneBc+BC_index)->Data_Set+(count-1))->identity+1)=0;
				*(((zone->ZoneBc+BC_index)->Data_Set+(count-1))->identity+2)=0;
				*(((zone->ZoneBc+BC_index)->Data_Set+(count-1))->identity+3)=0;
				 }
			if(c==1)
			  	{
				printf("\nEnter The Emissivity of this Surface\n");
				scanf("%lf",&emissivity);
			        (zone->ZoneBc+BC_index)->ZSlab->epsilon=emissivity;
			  	printf("\nEnter the constant value of Refractive Index\n");
				printf("As I am Assuming Refrective Index is constant On this Surface\n");
				scanf("%lf",((zone->ZoneBc+BC_index)->ZSlab->RefIndex));
				 
				((zone->ZoneBc+BC_index)->Data_Set+(count-1))->funct=NULL;
				((zone->ZoneBc+BC_index)->Data_Set+(count-1))->funct_beta=NULL;
				((zone->ZoneBc+BC_index)->Data_Set+(count-1))->funct_gama=NULL;
				*(((zone->ZoneBc+BC_index)->Data_Set+(count-1))->identity)=0;
				*(((zone->ZoneBc+BC_index)->Data_Set+(count-1))->identity+1)=0;
				*(((zone->ZoneBc+BC_index)->Data_Set+(count-1))->identity+2)=0;
				*(((zone->ZoneBc+BC_index)->Data_Set+(count-1))->identity+3)=0;
				 }
			  if(c==2)
			    {
			    if(base->BaseParameters->flagflowRadiation==1)
			    {
			     printf("\nEnter The Emissivity of this Surface\n");
			     scanf("%lf",&emissivity);
			     (zone->ZoneBc+BC_index)->ZSlab->epsilon=emissivity;
			      printf("Do you want to enter a function [1/yes or 0/no]\n");
			  	scanf("%d",&s);
				if(s==0)
				{
			     ((zone->ZoneBc+BC_index)->Data_Set+(count-1))->Data_Array=(double *)malloc(sizeof(double));
			     ((zone->ZoneBc+BC_index)->Data_Set+(count-1))->Data_Array=((zone->ZoneBc+BC_index)->Data_Set+find_nv(base,6))->Data_Array;
			     printf("%e\n",*(((zone->ZoneBc+BC_index)->Data_Set+(count-1))->Data_Array));
			     ((zone->ZoneBc+BC_index)->Data_Set+(count-1))->funct=NULL;
				((zone->ZoneBc+BC_index)->Data_Set+(count-1))->funct_beta=NULL;
				((zone->ZoneBc+BC_index)->Data_Set+(count-1))->funct_gama=NULL;
				*(((zone->ZoneBc+BC_index)->Data_Set+(count-1))->identity)=0;
				*(((zone->ZoneBc+BC_index)->Data_Set+(count-1))->identity+1)=0;
				*(((zone->ZoneBc+BC_index)->Data_Set+(count-1))->identity+2)=0;
				*(((zone->ZoneBc+BC_index)->Data_Set+(count-1))->identity+3)=0;
			     }
			     if(s==1)
				{
					if(*(((zone->ZoneBc+BC_index)->Data_Set+find_nv(base,6))->identity)==0)
				{
				*(((zone->ZoneBc+BC_index)->Data_Set+(count-1))->identity)=1;
				*(((zone->ZoneBc+BC_index)->Data_Set+(count-1))->identity+2)=0;
				*(((zone->ZoneBc+BC_index)->Data_Set+(count-1))->identity+3)=0;
				((zone->ZoneBc+BC_index)->Data_Set+(count-1))->funct=NULL;
				((zone->ZoneBc+BC_index)->Data_Set+(count-1))->funct_beta=NULL;
				((zone->ZoneBc+BC_index)->Data_Set+(count-1))->funct_gama=NULL;
				
				switch((zone->ZoneBc+BC_index)->Direction)
				{
				case 1:
				 maxn=nj*nk;
				 break;
				case 2:
				 maxn=nj*nk;
				 break;
				case 3:
				 maxn=ni*nk;
				 break;
				case 4:
				 maxn=ni*nk;
				 break;
				case 5:
				 maxn=ni*nj;
				 break;
				case 6:
				 maxn=ni*nj;
				 break;
				}
				((zone->ZoneBc+BC_index)->Data_Set+(count-1))->Data_Array=(double*)malloc((maxn)*sizeof(double));
				if(((zone->ZoneBc+BC_index)->Data_Set+find_nv(base,6))->Bc_Type==2||((zone->ZoneBc+BC_index)->Data_Set+find_nv(base,6))->Bc_Type==3)
					((zone->ZoneBc+BC_index)->Data_Set+(count-1))->Data_Array=((zone->ZoneBc+BC_index)->Data_Set+find_nv(base,6))->Data_Array;
				else
				((zone->ZoneBc+BC_index)->Data_Set+(count-1))->Data_Array=(zone->ZoneBc+BC_index)->ZSlab->T;
				printf("S=%e\n%d\n",*(((zone->ZoneBc+BC_index)->Data_Set+(count-1))->Data_Array),count);
				}
				if(*(((zone->ZoneBc+BC_index)->Data_Set+find_nv(base,6))->identity)==1)
				{
				*(((zone->ZoneBc+BC_index)->Data_Set+(count-1))->identity)=1;
				*(((zone->ZoneBc+BC_index)->Data_Set+(count-1))->identity+2)=0;
				*(((zone->ZoneBc+BC_index)->Data_Set+(count-1))->identity+3)=0;
				((zone->ZoneBc+BC_index)->Data_Set+(count-1))->funct=NULL;
				((zone->ZoneBc+BC_index)->Data_Set+(count-1))->funct_beta=NULL;
				((zone->ZoneBc+BC_index)->Data_Set+(count-1))->funct_gama=NULL;
				
				switch((zone->ZoneBc+BC_index)->Direction)
				{
				case 1:
				 maxn=nj*nk;
				 break;
				case 2:
				 maxn=nj*nk;
				 break;
				case 3:
				 maxn=ni*nk;
				 break;
				case 4:
				 maxn=ni*nk;
				 break;
				case 5:
				 maxn=ni*nj;
				 break;
				case 6:
				 maxn=ni*nj;
				 break;
				}
				((zone->ZoneBc+BC_index)->Data_Set+(count-1))->Data_Array=(double*)malloc((maxn)*sizeof(double));
				((zone->ZoneBc+BC_index)->Data_Set+(count-1))->Data_Array=((zone->ZoneBc+BC_index)->Data_Set+find_nv(base,6))->funct;
				}
				}
				printf("Enter the constant value of Refractive Index\n");
				printf("\nAs I am Assuming Refrective Index is constant On this Surface\n");
				scanf("%lf",((zone->ZoneBc+BC_index)->ZSlab->RefIndex));
			    }
			    if(base->BaseParameters->flagRadiation==1)
			    {
			     printf("\nEnter The Emissivity of this Surface\n");
			     scanf("%lf",&emissivity);
			     (zone->ZoneBc+BC_index)->ZSlab->epsilon=emissivity;
			     printf("Do you want to enter a function [1/yes or 0/no]\n");
			  	scanf("%d",&s);
					
			  while(!(s==0||s==1))
			  	{
				printf("\nPlease Enter 0 for no or 1 for yes\n");
				scanf("%d",&s);
				}
					 
				//printf("\n Do You Want to give Boundary condition in terms of Temperature or Radiative flux\n");
				//printf("Please Enter 1->Temperature Boundary Condition, 0-> Radiative Flux Boundary Condition\n");
				//scanf("%d",&(base->flagRadiativeTemperature));
				if(s==0)
				{
					((zone->ZoneBc+BC_index)->Data_Set+(count-1))->Data_Array=(double *)malloc(sizeof(double));
					
				
				if(base->BaseParameters->flagRadiativeTemperature==1)
					{
					 printf("\nEnter the Constant value of Temperature\n");
					 scanf("%lf",&value);
					 }
				 if(base->BaseParameters->flagRadiativeTemperature==0)
					 {
					 printf("\nEnter the Constant Value of Radiative Flux\n");
					 scanf("%lf",&value);
					 }
				 *((zone->ZoneBc+BC_index)->Data_Set+(count-1))->Data_Array=value;
					
				
				 
				 
					 
				 printf("\nEnter the constant value of Refractive Index\n");
				 printf("As I am Assuming Refrective Index is constant On this Surface\n");
				 scanf("%lf",((zone->ZoneBc+BC_index)->ZSlab->RefIndex));
					
					
				
				
				((zone->ZoneBc+BC_index)->Data_Set+(count-1))->funct=NULL;
				((zone->ZoneBc+BC_index)->Data_Set+(count-1))->funct_beta=NULL;
				((zone->ZoneBc+BC_index)->Data_Set+(count-1))->funct_gama=NULL;
				*(((zone->ZoneBc+BC_index)->Data_Set+(count-1))->identity)=0;
				*(((zone->ZoneBc+BC_index)->Data_Set+(count-1))->identity+1)=0;
				*(((zone->ZoneBc+BC_index)->Data_Set+(count-1))->identity+2)=0;
				*(((zone->ZoneBc+BC_index)->Data_Set+(count-1))->identity+3)=0;
				}
				if(s==1)
				{
				*(((zone->ZoneBc+BC_index)->Data_Set+(count-1))->identity)=1;
				*(((zone->ZoneBc+BC_index)->Data_Set+(count-1))->identity+2)=0;
				*(((zone->ZoneBc+BC_index)->Data_Set+(count-1))->identity+3)=0;
				((zone->ZoneBc+BC_index)->Data_Set+(count-1))->funct_beta=NULL;
				((zone->ZoneBc+BC_index)->Data_Set+(count-1))->funct_gama=NULL;
				printf("\nENTER THE FUNCTION OF YOUR CHOICE FOR TEMPERATURE/FLUX\n");
				scanf("%d",&choice);
					
				printf("Enter the constant value of Refractive Index\n");
				printf("\nAs I am Assuming Refrective Index is constant On this Surface\n");
					
				scanf("%lf",((zone->ZoneBc+BC_index)->ZSlab->RefIndex));
				switch((zone->ZoneBc+BC_index)->Direction)
				{
				case 1:
				 maxn=nj*nk;
				 break;
				case 2:
				 maxn=nj*nk;
				 break;
				case 3:
				 maxn=ni*nk;
				 break;
				case 4:
				 maxn=ni*nk;
				 break;
				case 5:
				 maxn=ni*nj;
				 break;
				case 6:
				 maxn=ni*nj;
				 break;
				}
				((zone->ZoneBc+BC_index)->Data_Set+(count-1))->funct=(double*)malloc((maxn)*sizeof(double));
				Call_user_func(choice,((zone->ZoneBc+BC_index)->Data_Set+(count-1))->funct,zone,(zone->ZoneBc+BC_index),0.0);
				}
			}
			}
			}
		
			
		}
		}
	return 0;
	}

