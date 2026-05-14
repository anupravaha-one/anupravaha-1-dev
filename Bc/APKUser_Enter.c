#include "APKUser_Enter.h"

/**
Function user enter is used to store a single valued boundary condition ie, constant values.
The function input parameters are Bc type (Drichelet, Nueman etc.) and the name of variable. 
This is also a user interactive function. Switch case is used for various Bc types.  
Inside any case for a particular BC relevent values are enterd by the user in a local variable. 
The pointer to this variable is returned to the calling function. 
*/

double *User_Enter(Base *base,int Bctype,char name)
{
        FILE *var2,*var3;
	double *variable;
        
        var2=fopen("Settings/coefficients.dat","a");	
        var3=fopen("Settings/ProgramReadInputValues.dat","a");
	switch(Bctype)
	{
	
	case 1:
	printf("\nENTER THE INWARD GRADIENT OF %c\n",name);
	variable=(double*)malloc(sizeof(double));
	
	if(variable==NULL)
	{
	printf("\nMemory Allocation failed. Not enough memory available exiting from APKUser_Enter\n");
	exit(1);
	}
	scanf("%lf",variable);        
        fprintf(var2,"%lf\n",*variable);        
        fprintf(var3,"%lf \t \n\n",*variable);
	break;
	
	case 10:
	printf("\nENTER THE INWARD FLUX OF %c\n",name);
	variable=(double*)malloc(sizeof(double));
	
	
	if(variable==NULL)
	{
	printf("\nMemory Allocation failed. Not enough memory available exiting from APKUser_Enter\n");
	exit(1);
	}
	scanf("%lf",variable);        
        fprintf(var2,"%lf\n",*variable);        
        fprintf(var3,"%lf \t \n\n",*variable);
	break;
	
	case 0:
	variable=(double*)malloc(sizeof(double));
	*variable=0;
	break;
	
	case 2:
	printf("\nENTER THE VALUE OF %c AT THE BOUNDARY NODE\n",name);
	variable=(double*)malloc(sizeof(double));
	
	if(variable==NULL)
	{
	printf("\nMemory Allocation failed. Not enough memory available exiting from APKUser_Enter\n");
	exit(1);
	}
	
	scanf("%lf",variable);        
        fprintf(var2,"%lf\n",*variable);        
        fprintf(var3,"%lf \t \n\n",*variable);
	break;
	
	case 3:
	variable=(double*)malloc(sizeof(double));
	*variable=0;
	break;
	
	case 41:
	if(base->BaseParameters->flagflowRadiation==0)
	{
	printf("\nENTER THE VARIABLE COFFICIENTS, COEFFICIENT OF GRADIENT AND BOUNDARY NODE VALUE OR GAMMA IN THE SAME ORDER FOR  %c\n",name);
	variable=(double*)malloc(3*sizeof(double));
	
	if(variable==NULL)
	{
	printf("\nMemory Allocation failed. Not enough memory available exiting from APKUser_Enter\n");
	exit(1);
	}
	
	scanf("%lf%lf%lf",(variable),(variable+1),(variable+2));        
        fprintf(var2,"%lf %lf %lf\n",*(variable),*(variable+1),*(variable+2));        
        fprintf(var3,"%lf %lf %lf \t \n\n",*(variable),*(variable+1),*(variable+2));             
	}
	if(base->BaseParameters->flagflowRadiation==1)
	{
	printf("\nENTER THE Flux  %c\n",name);
	variable=(double*)malloc(sizeof(double));
	
	if(variable==NULL)
	{
	printf("\nMemory Allocation failed. Not enough memory available exiting from APKUser_Enter\n");
	exit(1);
	}
	
	scanf("%lf",variable);        
        fprintf(var2,"%lf\n",*variable);        
        fprintf(var3,"%lf \t \n\n",*variable);
	}
	break;
	case 40:
	if(base->BaseParameters->flagflowRadiation==0)
	{
	printf("\nENTER THE VARIABLE COFFICIENTS, COEFFICIENT OF GRADIENT AND BOUNDARY NODE VALUE OR GAMMA IN THE SAME ORDER FOR  %c\n",name);
	variable=(double*)malloc(3*sizeof(double));
	
	if(variable==NULL)
	{
	printf("\nMemory Allocation failed. Not enough memory available exiting from APKUser_Enter\n");
	exit(1);
	}
	
	scanf("%lf%lf%lf",(variable),(variable+1),(variable+2));        
        fprintf(var2,"%lf %lf %lf\n",*(variable),*(variable+1),*(variable+2));        
        fprintf(var3,"%lf %lf %lf \t \n\n",*(variable),*(variable+1),*(variable+2));             
	}
	if(base->BaseParameters->flagflowRadiation==1)
	{
	printf("\nENTER THE Flux  %c\n",name);
	variable=(double*)malloc(sizeof(double));
	
	if(variable==NULL)
	{
	printf("\nMemory Allocation failed. Not enough memory available exiting from APKUser_Enter\n");
	exit(1);
	}
	
	scanf("%lf",variable);        
        fprintf(var2,"%lf\n",*variable);        
        fprintf(var3,"%lf \t \n\n",*variable);
	}
	break;
	case 5:
	printf("\nEnter the Average Velocity at the outlet \n");
        variable=(double*)malloc(sizeof(double));
	
	
	if(variable==NULL)
	{
	printf("\nMemory Allocation failed. Not enough memory available exiting from APKUser_Enter\n");
	exit(1);
	}
	scanf("%lf",variable);        
        fprintf(var2,"%lf\n",*variable);        
        fprintf(var3,"%lf \t \n\n",*variable);
	break;
	case 70:
	printf("\nENTER THE COFFICIENTS OF GRADIENT AND OF THE BOUNDARY NODE AND ALSO VALUE OF GAMMA IN THE SAME ORDER FOR  %c\n",name);
	variable=(double*)malloc(3*sizeof(double));
	
	if(variable==NULL)
	{
	printf("\nMemory Allocation failed. Not enough memory available exiting from APKUser_Enter\n");
	exit(1);
	}
	scanf("%lf%lf%lf",(variable),(variable+1),(variable+2));        
        fprintf(var2,"%lf %lf %lf\n",*(variable),*(variable+1),*(variable+2));        
        fprintf(var3,"%lf %lf %lf \t \n\n",*(variable),*(variable+1),*(variable+2));             
	
	break;
	
	case 71:
	printf("\nENTER THE COFFICIENTS OF GRADIENT AND OF THE BOUNDARY NODE AND ALSO VALUE OF GAMMA IN THE SAME ORDER FOR  %c\n",name);
	variable=(double*)malloc(3*sizeof(double));
	
	if(variable==NULL)
	{
	printf("\nMemory Allocation failed. Not enough memory available exiting from APKUser_Enter\n");
	exit(1);
	}
	
	scanf("%lf%lf%lf",(variable),(variable+1),(variable+2));        
        fprintf(var2,"%lf %lf %lf\n",*(variable),*(variable+1),*(variable+2));           
        fprintf(var3,"%lf %lf %lf \t \n\n",*(variable),*(variable+1),*(variable+2));             
	
	break;
	case 9:
	printf("zero flux Boundary Condition");
	variable=(double*)malloc(sizeof(double));
	break;
	case 90:
	printf("zero flux Boundary Condition");
	variable=(double*)malloc(sizeof(double));
	break;
	default:
	printf("\nINCORRECT INPUT\n");
	printf("\nPLEASE REFER THE ABOVE OPTIONS AND ENTER THE RIGHT SELECTION\n");
	variable=NULL;
	break;
	}
        fclose(var2);
        fclose(var3);
	return(variable);
}


