#include "APKInputBoundarycondition.h"

/**
The function InputBoundarycondition() interacts with user to store the desired boundary condition on a 
particular BC patch. This function is called inside a zone loop for each zone present in the computational 
domain.

The input parameters to this function are, count of total no of variables, pointer of the zone, flag for 2D 
solver. 

It loops over all the BC patches present in the zone. To give an idea of the BC patch loacation BC tag, 
point range of the Bc patch and corrosponding vertex coordinates are printed for this purpose point_coord() 
is called. 
 
Checks one by one the entire range of possible varibles, offered by this solver 
if they are initialized by the user (ie. user wants to solve them)
 
If a particular variable is initialize, for that varible BC options (Neumann Dirichlet Homogenous Dirichlet, 
Robin etc.) are printed, the users responce along with variable name is then stored in the data structue under 
the node Data_set.
 
If user wants to enter a consatant value User_Enter() is called.

If user wants to enter a function Call_user_func() is called. 
  
If 2D flag is positive (ie. solution is 2D) then default BC for w component of velocity at all the patches 
is assigned. Also default boundary conditions at BC patch Kmax and Kmin for all the variables are assigned.
 */

int InputBoundarycondition(Base *base,int NS,Zone *CgnsZone,int flag2d)
{
        FILE *var2,*var3;
	double dum;
	char comp1[]="Kmin";
	char comp2[]="Kmax";	
	char *varname;
	ZoneBoundary *zbc;
	zbc=CgnsZone->ZoneBc;
	int RadVar,che;
	
	/*char  variable[]={'u','v','w','P','T','\0'};
	char scal[10][4]={    {'0','s','\0'},
	{'1','s','\0'},
	{'2','s','\0'},
	{'3','s','\0'},
	{'4','s','\0'},
	{'5','s','\0'},
	{'6','s','\0'},
	{'7','s','\0'},
	{'8','s','\0'},
	{'9','s','\0'}
};*/
	int c,BC_index,s,choice,maxn,nki,njk,nij,ni,nj,nk,rc,ts=0,c1,c2,pa,conx,fof;
	int count=0,test,number,nv,nop;
	int count1,size;
	double *domain,*slab;
              
	nop=base->nphase;
	ni=*(CgnsZone->Size);
	nj=*(CgnsZone->Size+1);
	nk=*(CgnsZone->Size+2); 
	nij=ni*nj;
	nki=nk*ni;
	njk=nj*nk;
	if(nij>njk)
	{
		if(nij>nki)
			maxn=nij;

		else
			maxn=nki;
	}
	else
	{	if(njk>nki)
			maxn=njk;

		else
			maxn=nki;
	}
	printf("\n\nNumber of variables are %d\n\n",NS);
	printf("\n~~~~~~~~THIS IS ZONE: %s~~~~~~~~",CgnsZone->zonename);
	
	/*The looping over all boundary patches starts here */
	RadVar = base->Radvar;
	for(BC_index=0;BC_index<(CgnsZone->Number_of_BC_Patches);BC_index++)
	{
		count=0;
		/*The memory space for Dataset is alloted here*/
		(zbc+BC_index)->Data_Set=(Dataset*)malloc((NS+RadVar)*sizeof(Dataset));
		c1=strcmp((zbc+BC_index)->BCName,comp1);
		c2=strcmp((zbc+BC_index)->BCName,comp2);
		printf("\n~~~~~~~~THIS IS BOUNDARY PATCH NAMED: %s~~~~~~~~\n",(zbc+BC_index)->BCName);
		/*This function gives the Tag given to the boundary patch by the user while grid generation */
		print_tagname((zbc+BC_index)->Boundary_Tag);
		/*This function converts point range into X,Y,Z coordinates*/
		point_coord(CgnsZone,(zbc+BC_index)->Point_Range);	
		/*The check statement for checking the memory allocation from initialize function*/
		for (number=0;number<NS;number++)
		{
			count=count+1;
			nv=base->VARTAG[count-1];//BY SUSHIL
			//Added by sushil for multiphase flow...
			
			if ((base->VARTAG[count-1]!=31)&&(base->VARTAG[count-1]!=32)&&(base->VARTAG[count-1]!=33))//BY SUSHIL
			{	
				//End 
				
				varname = find_name(base,count);
				strcpy(((CgnsZone->ZoneBc+BC_index)->Data_Set+(count-1))->Name_of_variable,varname);
				if(c1==0||c2==0 ||(flag2d==1 && ((base->VARTAG[count-1]==3)||(base->VARTAG[count-1]==33)||(base->VARTAG[count-1]==43)||(base->VARTAG[count-1]==46))))
				{
				    c=0;
				}
				else
				{
					printf("\n~~~~~~~~THE BOUNDARY CONDITION for %s\n",find_name(base,count));
					
					if ((base->VARTAG[count-1]>=41)&&(base->VARTAG[count-1]<=55))
					{
						printf("\n\n Enter velocity boundary condition for each phase respectively\n");
					}	
					if ((base->VARTAG[count-1]==17))
					  {	printf("\n Do you want to apply the zero efflux condition for current [1/0] \n");
                                                scanf("%d",&c);    
                                                
                                                var2=fopen("Settings/coefficients.dat","a");   	
                                                var3=fopen("Settings/ProgramReadInputValues.dat","a");   
                                                fprintf(var2,"%d\n",c);        
                                                fprintf(var3,"%d \tDo you want to apply the zero efflux condition for current [Yes 1/No 0] \n\n",c);
                                                fclose(var2);        
                                                fclose(var3);
                                                
                                             	if(c==1)
						  {
                                                   printf("\n Specify whether applying orthogonal [0] / Non-orthogonal [1] patch \n");
                                                   scanf("%d",&pa); 
                                                   
                                                   var2=fopen("Settings/coefficients.dat","a");                                                 
                                                   var3=fopen("Settings/ProgramReadInputValues.dat","a");          
                                                   fprintf(var2,"%d\n",pa);        
                                                   fprintf(var3,"%d \tSpecify whether applying orthogonal [0] / Non-orthogonal [1] patch \n\n",pa);
                                                   fclose(var2);
                                                   fclose(var3);
                                                   
                                                   c = 81*pa + 9;

                                                  }
                                                else
                                                  {
	                                           //printf("\nENTER  1 for  Neumann..\n 2 for Dirichlet..\n 3 for Homogenous Dirichlet  \n");
						   printf("\nFor Neumann Boundary Conditions.................[Enter 1]");
						   printf("\nFor Dirichlet Boundary Conditions...............[Enter 2]");
						   printf("\nFor Homogenous Dirichlet Boundary Conditions....[Enter 3]\n");
					               	do 
				                	{
					                	//scanf("%d",&c);
					                	test=check_bc_type(base,(CgnsZone->ZoneBc+BC_index),&c,base->VARTAG[count-1]);
				                	}while(test!=1);
						  }
					  }
                                        else{
					//printf("\nENTER  1 for  Neumann..\n 2 for Dirichlet..\n 3 for Homogenous Dirichlet ..\n 4 for Robin/Including Radiative Flux..\n5 for Orlanski\n7 for TimeDependent \n");
					printf("\nFor Neumann Boundary Conditions.................[Enter 1]");
					printf("\nFor Dirichlet Boundary Conditions...............[Enter 2]");
					printf("\nFor Homogenous Dirichlet Boundary Conditions....[Enter 3]");
					printf("\nFor Robin Boundary Conditions...................[Enter 4]");
					printf("\nFor Orlanski Boundary Conditions................[Enter 5]");
					printf("\nFor TimeDependent Boundary Conditions.......... [Enter 7]\n");
	
					do 
					{
						//scanf("%d",&c);
						test=check_bc_type(base,(CgnsZone->ZoneBc+BC_index),&c,base->VARTAG[count-1]);
					}while(test!=1);
					}		
				}
	
                                         
				if(c==1)
				{
					printf("ENTER THE FOLLOWING INFORMATION FOR NEUMANN BC\n");
					printf("if you have ORTHOGONAL grid on the boundary patch AND zero gradient:................ [Enter 0]\n" );
					printf("if you have ORTHOGONAL grid on the boundary patch AND non-zero gradient:............ [Enter 1]\n" );
					printf("if you have NON ORTHOGONAL grid on the boundary patch AND non-zero OR zero gradient: [Enter 10]\n" );
					scanf("%d",&c); 
                                        
                                        var2=fopen("Settings/coefficients.dat","a");
                                        var3=fopen("Settings/ProgramReadInputValues.dat","a");                                        
                                        fprintf(var2,"%d\n",c);
                                        fprintf(var3,"%d \tNEUMANN BC ORTHOGONAL grid zero gradient 0, non-zero gradient 1, NON ORTHOGONAL grid 10\n\n",c);
                                        fclose(var2);
                                        fclose(var3);
	
				}
				if(c==4)
				{
					printf("ENTER THE FOLLOWING INFORMATION FOR ROBIN BOUNDARY CONDITION\n");
					printf("If you have an ORTHOGONAL grid on the Boundary Patch then:........................ [Enter 40]\n");
					printf("If you have a NON ORTHOGONAL grid on the Boundary Patch then:..................... [Enter 41]\n");
					scanf("%d",&c);                                        
                                        
                                        var2=fopen("Settings/coefficients.dat","a");
                                        var3=fopen("Settings/ProgramReadInputValues.dat","a");                                        
                                        fprintf(var2,"%d\n",c);
                                        fprintf(var3,"%d \tROBIN BC ORTHOGONAL grid 40, NON ORTHOGONAL grid 41\n\n",c);
                                        fclose(var2);
                                        fclose(var3);
				}
				if(c==7)
				{
					((CgnsZone->ZoneBc+BC_index)->Data_Set+(count-1))->flagTime=1;
					printf("THE USER HAS ENTERED A TIME DEPENDENT BOUNDARY CONDITION\n");		
					printf("BY DEFAULT TIME DEPENDENT BOUNDARY CONDITION IS ROBIN BOUNDARY SO PLEASE ENTER THE FOLLOWING INFORMATION\n");
					printf("If you have an ORTHOGONAL grid on the Boundary Patch then:........................ [Enter 70]\n");
					printf("If you have a NON ORTHOGONAL grid on the Boundary Patch then:..................... [Enter 71]\n");
					scanf("%d",&c);                                        
                                       
                                        var2=fopen("Settings/coefficients.dat","a");
                                        var3=fopen("Settings/ProgramReadInputValues.dat","a");                                        
                                        fprintf(var2,"%d\n",c);
                                        fprintf(var3,"%d \tTIME DEPENDENT BC ORTHOGONAL grid 70, NON ORTHOGONAL grid 71\n\n",c);
                                        fclose(var2);
                                        fclose(var3);
				}
				else
				{
					((CgnsZone->ZoneBc+BC_index)->Data_Set+(count-1))->flagTime=0;
				}
				if(c==1||c==2||c==41||c==40||c==10)	
				{
					printf("Do you want to enter a Function [1/yes or 0/no] ");
	                           	printf("\nOpt 'Function' for BC input from file\n");
					do 
					{
						scanf("%d",&s);
						test=check_function(s);
					}while(test!=1);	                                        
                                        
                                        var2=fopen("Settings/coefficients.dat","a");
                                        var3=fopen("Settings/ProgramReadInputValues.dat","a");                                        
                                        fprintf(var2,"%d\n",s);
                                        fprintf(var3,"%d \tDo you want to enter a Function [1/yes or 0/no] \n\n",s);
                                        fclose(var2);
                                        fclose(var3);
	
				}
				else if(c==71||c==70)
				{
					s=1;
				}
				else
				{
					s=0;
				}
				
				((CgnsZone->ZoneBc+BC_index)->Data_Set+(count-1))->Bc_Type=c;
// 				printf("\n BCTYPE=%d",c);                                
				
				if(s==0)
				{
	
					((CgnsZone->ZoneBc+BC_index)->Data_Set+(count-1))->Data_Array=(double *)User_Enter(base,((CgnsZone->ZoneBc+BC_index)->Data_Set+(count-1))->Bc_Type,(((CgnsZone->ZoneBc+BC_index)->Data_Set+(count-1))->Name_of_variable[0]));
					((CgnsZone->ZoneBc+BC_index)->Data_Set+(count-1))->funct=NULL;
					((CgnsZone->ZoneBc+BC_index)->Data_Set+(count-1))->funct_beta=NULL;
					((CgnsZone->ZoneBc+BC_index)->Data_Set+(count-1))->funct_gama=NULL;
					*(((CgnsZone->ZoneBc+BC_index)->Data_Set+(count-1))->identity)=0;
					*(((CgnsZone->ZoneBc+BC_index)->Data_Set+(count-1))->identity+1)=0;
					*(((CgnsZone->ZoneBc+BC_index)->Data_Set+(count-1))->identity+2)=0;
					*(((CgnsZone->ZoneBc+BC_index)->Data_Set+(count-1))->identity+3)=0;
				}
	
				else if(s==1)
				{	                                  
					*(((CgnsZone->ZoneBc+BC_index)->Data_Set+(count-1))->identity)=1;
	
					//system("gvim Functions.dat");
					if(c==41||c==40||c==71||c==70)
					{
						printf("\nDO YOU WANT TO ENTER  A FUNCTION IN SPACE CO-ORDINATES FOR ALPHA [1/yes or 0/no]\n");
						do 
						{
							scanf("%d",&rc);
							test=check_function(rc);
						}while(test!=1);                                        
                                                
                                                var2=fopen("Settings/coefficients.dat","a");
                                                var3=fopen("Settings/ProgramReadInputValues.dat","a");                                        
                                                fprintf(var2,"%d\n",rc);
                                                fprintf(var3,"%d \tDO YOU WANT TO ENTER  A FUNCTION IN SPACE CO-ORDINATES FOR ALPHA [1/yes or 0/no] \n\n",rc);
                                                fclose(var2);
                                                fclose(var3);
	
						if(rc==1)
						{
							*(((CgnsZone->ZoneBc+BC_index)->Data_Set+(count-1))->identity+1)=1;
							printf("\nENTER THE RESPECTIVE INDEX OF FUNCTION OF YOUR CHOICE FOR ALPHA\n");
							scanf("%d",&choice);
                                                        
                                                        var2=fopen("Settings/coefficients.dat","a");
                                                        var3=fopen("Settings/ProgramReadInputValues.dat","a");                                        
                                                        fprintf(var2,"%d\n",choice);
                                                        fprintf(var3,"%d \tENTER THE RESPECTIVE INDEX OF FUNCTION OF YOUR CHOICE FOR ALPHA \n\n",choice);
                                                        fclose(var2);
                                                        fclose(var3);
                                                
							*(((CgnsZone->ZoneBc+BC_index)->Data_Set+(count-1))->FunctionIndex+0)=choice;
							((CgnsZone->ZoneBc+BC_index)->Data_Set+(count-1))->funct=(double*)malloc((maxn)*sizeof(double));
	
							if(((CgnsZone->ZoneBc+BC_index)->Data_Set+(count-1))->funct==NULL)
							{
								printf("MEMORY ALLOCATION FAILED");
								exit(1);
							}
							Call_user_func(choice,((CgnsZone->ZoneBc+BC_index)->Data_Set+(count-1))->funct,CgnsZone,(CgnsZone->ZoneBc+BC_index),0.0);
	
						}
						else
						{
	
							printf("\n ENTER THE CONSTANT VALUE OF ALPHA\n");
							if(ts==0){((CgnsZone->ZoneBc+BC_index)->Data_Set+(count-1))->Data_Array=(double*)malloc(3*sizeof(double));}
							scanf("%lf",&dum);
                                                        
                                                        var2=fopen("Settings/coefficients.dat","a");
                                                        var3=fopen("Settings/ProgramReadInputValues.dat","a");                                        
                                                        fprintf(var2,"%lf\n",dum);
                                                        fprintf(var3,"%lf \tENTER THE CONSTANT VALUE OF ALPHA \n\n",dum);
                                                        fclose(var2);
                                                        fclose(var3);
                                                        
							*(((CgnsZone->ZoneBc+BC_index)->Data_Set+(count-1))->Data_Array)=dum;
							ts=ts+1;
							((CgnsZone->ZoneBc+BC_index)->Data_Set+(count-1))->funct=NULL;
							*(((CgnsZone->ZoneBc+BC_index)->Data_Set+(count-1))->identity+1)=0;
						}	
						printf("\nDO YOU WANT TO ENTER  A FUNCTION IN SPACE CO-ORDINATES FOR BETA [1/yes or 0/no]\n");
						do 
						{
							scanf("%d",&rc);
							test=check_function(rc);
						}while(test!=1);                                        
                                                
                                                var2=fopen("Settings/coefficients.dat","a");
                                                var3=fopen("Settings/ProgramReadInputValues.dat","a");                                        
                                                fprintf(var2,"%d\n",rc);
                                                fprintf(var3,"%d \tDO YOU WANT TO ENTER  A FUNCTION IN SPACE CO-ORDINATES FOR BETA [1/yes or 0/no] \n\n",rc);
                                                fclose(var2);
                                                fclose(var3);
	
						if(rc==1)
						{
							*(((CgnsZone->ZoneBc+BC_index)->Data_Set+(count-1))->identity+2)=1;
							printf("\nENTER THE RESPECTIVE INDEX OF FUNCTION OF YOUR CHOICE FOR BETA\n");
							scanf("%d",&choice);                                        
                                                        
                                                        var2=fopen("Settings/coefficients.dat","a");
                                                        var3=fopen("Settings/ProgramReadInputValues.dat","a");
                                                        fprintf(var2,"%d\n",choice);
                                                        fprintf(var3,"%d \tENTER THE RESPECTIVE INDEX OF FUNCTION OF YOUR CHOICE FOR BETA \n\n",choice);
                                                        fclose(var2);
                                                        fclose(var3);
        
							*(((CgnsZone->ZoneBc+BC_index)->Data_Set+(count-1))->FunctionIndex+1)=choice;
							((CgnsZone->ZoneBc+BC_index)->Data_Set+(count-1))->funct_beta=(double*)malloc((maxn)*sizeof(double));
							if(((CgnsZone->ZoneBc+BC_index)->Data_Set+(count-1))->funct_beta==NULL)
							{
								printf("MEMORY ALLOCATION FAILED");
								exit(1);
							}
							Call_user_func(choice,((CgnsZone->ZoneBc+BC_index)->Data_Set+(count-1))->funct_beta,CgnsZone,(CgnsZone->ZoneBc+BC_index),0.0);
	
						}
						else
						{
							printf("\n ENTER THE CONSTANT VALUE OF BETA\n");
							if(ts==0){((CgnsZone->ZoneBc+BC_index)->Data_Set+(count-1))->Data_Array=(double*)malloc(3*sizeof(double));}
							scanf("%lf",&dum);
                                                        
                                                        
                                                        var2=fopen("Settings/coefficients.dat","a");
                                                        var3=fopen("Settings/ProgramReadInputValues.dat","a");
                                                        fprintf(var2,"%lf\n",dum);
                                                        fprintf(var3,"%lf \tENTER THE CONSTANT VALUE OF BETA \n\n",dum);
                                                        fclose(var2);
                                                        fclose(var3);
                                                        
							*(((CgnsZone->ZoneBc+BC_index)->Data_Set+(count-1))->Data_Array+1)=dum;
							ts=ts+1;
							((CgnsZone->ZoneBc+BC_index)->Data_Set+(count-1))->funct_beta=NULL;
							*(((CgnsZone->ZoneBc+BC_index)->Data_Set+(count-1))->identity+2)=0;
						}	
						printf("\nDO YOU WANT TO ENTER  A FUNCTION IN SPACE CO-ORDINATES FOR GAMA [1/yes or 0/no]\n");
						do 
						{
							scanf("%d",&rc);
							test=check_function(rc);
						}while(test!=1);                                        
                                                
	                                        var2=fopen("Settings/coefficients.dat","a");
                                                var3=fopen("Settings/ProgramReadInputValues.dat","a");
                                                fprintf(var2,"%d\n",rc);
                                                fprintf(var3,"%d \tDO YOU WANT TO ENTER  A FUNCTION IN SPACE CO-ORDINATES FOR GAMA [1/yes or 0/no] \n\n",rc);
                                                fclose(var2);
                                                fclose(var3);
                                                        
						if(rc==1)
						{
							*(((CgnsZone->ZoneBc+BC_index)->Data_Set+(count-1))->identity+3)=1;
							printf("\nENTER THE RESPECTIVE INDEX OF FUNCTION OF YOUR CHOICE FOR GAMA\n");
							scanf("%d",&choice);
                                                        
                                                        var2=fopen("Settings/coefficients.dat","a");
                                                        var3=fopen("Settings/ProgramReadInputValues.dat","a");
                                                        fprintf(var2,"%d\n",choice);
                                                        fprintf(var3,"%d \tENTER THE RESPECTIVE INDEX OF FUNCTION OF YOUR CHOICE FOR GAMA \n\n",choice);
                                                        fclose(var2);
                                                        fclose(var3);
                                                        
							*(((CgnsZone->ZoneBc+BC_index)->Data_Set+(count-1))->FunctionIndex+2)=choice;
							((CgnsZone->ZoneBc+BC_index)->Data_Set+(count-1))->funct_gama=(double*)malloc((maxn)*sizeof(double));
							if(((CgnsZone->ZoneBc+BC_index)->Data_Set+(count-1))->funct_gama==NULL)
							{
								printf("MEMORY ALLOCATION FAILED");
								exit(1);
							}
							Call_user_func(choice,((CgnsZone->ZoneBc+BC_index)->Data_Set+(count-1))->funct_gama,CgnsZone,(CgnsZone->ZoneBc+BC_index),0.0);
	
						}
						else
						{
							printf("\n ENTER THE CONSTANT VALUE OF GAMA\n");
							if(ts==0){((CgnsZone->ZoneBc+BC_index)->Data_Set+(count-1))->Data_Array=(double*)malloc(3*sizeof(double));}
							scanf("%lf",&dum);                                        
                                                                                                                
                                                        var2=fopen("Settings/coefficients.dat","a");
                                                        var3=fopen("Settings/ProgramReadInputValues.dat","a");
                                                        fprintf(var2,"%lf\n",dum);
                                                        fprintf(var3,"%lf \tENTER THE CONSTANT VALUE OF GAMA \n\n",dum);
                                                        fclose(var2);
                                                        fclose(var3);
                                                        
							*(((CgnsZone->ZoneBc+BC_index)->Data_Set+(count-1))->Data_Array+2)=dum;
							((CgnsZone->ZoneBc+BC_index)->Data_Set+(count-1))->funct_gama=NULL;
							*(((CgnsZone->ZoneBc+BC_index)->Data_Set+(count-1))->identity+3)=0;
	
						}	
					}
					else
					{	*(((CgnsZone->ZoneBc+BC_index)->Data_Set+(count-1))->identity+1)=1;
					*(((CgnsZone->ZoneBc+BC_index)->Data_Set+(count-1))->identity+2)=0;
					*(((CgnsZone->ZoneBc+BC_index)->Data_Set+(count-1))->identity+3)=0;
					((CgnsZone->ZoneBc+BC_index)->Data_Set+(count-1))->funct_beta=NULL;
					((CgnsZone->ZoneBc+BC_index)->Data_Set+(count-1))->funct_gama=NULL;
	                                printf("\n Read BC data from readu4bc.dat or User defined function from APKuser_defined_funcs.c file? [readu4bc.dat 1/APKuser_defined_funcs.c 0] \n");
					scanf("%d",&fof);                                        
                                        
                                        var2=fopen("Settings/coefficients.dat","a");
                                        var3=fopen("Settings/ProgramReadInputValues.dat","a");
                                        fprintf(var2,"%d\n",fof);
                                        fprintf(var3,"%d \tBC data from readu4bc.dat or User defined function from APKuser_defined_funcs.c file? [readu4bc.dat 1/APKuser_defined_funcs.c 0\n\n",fof);
                                        fclose(var2);
                                        fclose(var3);
        
					if(fof==1)
                                            { 
					  choice=20;
                                                                          
                                            }
					else{
					printf("\nENTER THE RESPECTIVE INDEX OF FUNCTION OF YOUR CHOICE\n");
					scanf("%d",&choice);                                                                                               
                                        
                                        var2=fopen("Settings/coefficients.dat","a");
                                        var3=fopen("Settings/ProgramReadInputValues.dat","a");
                                        fprintf(var2,"%d\n",choice);
                                        fprintf(var3,"%d \tENTER THE RESPECTIVE INDEX OF FUNCTION OF YOUR CHOICE \n\n",choice);
                                        fclose(var2);
                                        fclose(var3);
                                     
					}
					((CgnsZone->ZoneBc+BC_index)->Data_Set+(count-1))->funct=(double*)malloc((maxn)*sizeof(double));
					Call_user_func(choice,((CgnsZone->ZoneBc+BC_index)->Data_Set+(count-1))->funct,CgnsZone,(CgnsZone->ZoneBc+BC_index),0.0);
	
					}                                	        
				}                                
	
			}
			else
			{	//Added by sushil for multiphase flow variables
				//Allocate memory for Data Array for vartags 31,32,33. and initialize its value to zero.
				((CgnsZone->ZoneBc+BC_index)->Data_Set+(count-1))->Data_Array=(double *)malloc(3*sizeof(double));
// 				((CgnsZone->ZoneBc+BC_index)->Data_Set+(count-1))->Bc_Type = 2;
				*(((CgnsZone->ZoneBc+BC_index)->Data_Set+(count-1))->Data_Array)=0;
				*(((CgnsZone->ZoneBc+BC_index)->Data_Set+(count-1))->Data_Array+1)=0;
				*(((CgnsZone->ZoneBc+BC_index)->Data_Set+(count-1))->Data_Array+2)=0;
				
			}
			
			//Added by sushil
			
// 			printf("\nBCindex=%d",*(zbc+BC_index));
			//update_boundary( base, nv, CgnsZone, zbc+BC_index);
			if((base->flagMPFV==1)&&(base->flagTurb==0))//changed from  base->flagMPFV==1 to (base->flagMPFV==1)&&(base->flagTurb==0) to account only for laminar MPFV--chandrasekhar
			{
				if(base->VARTAG[count-1]==41)
				{
					count1=find_nv(base,31);
					((CgnsZone->ZoneBc+BC_index)->Data_Set+(count1))->Bc_Type=c;
				}
				if(base->VARTAG[count-1]==42)
				{
					count1=find_nv(base,32);	
					((CgnsZone->ZoneBc+BC_index)->Data_Set+(count1))->Bc_Type=c;
				}
				if(base->VARTAG[count-1]==43)
				{
					count1=find_nv(base,33);
					((CgnsZone->ZoneBc+BC_index)->Data_Set+(count1))->Bc_Type=c;
				}
				
				if ((zbc+BC_index)->Direction==1||(zbc+BC_index)->Direction==2)
					size=(nk-1)*(nj-1);
				if ((zbc+BC_index)->Direction==3||(zbc+BC_index)->Direction==4)
					size=(nk-1)*(ni-1);
				if ((zbc+BC_index)->Direction==5||(zbc+BC_index)->Direction==6)
					size=(nj-1)*(ni-1);
				
				domain = MakeDomainArray(CgnsZone,nv);

				slab   = FindSlabPtrBC((zbc+BC_index),nv);
				
				boundary_d_t_s(base,(zbc+BC_index)->Direction,domain,nv,CgnsZone,size,((zbc+BC_index)->Data_Set+find_nv(base,nv)),slab);
			}
			
			//Added by chandrasekhar --- to assign the same boundary tags for Global momentum as that of phasic momentum for Turbulent MPFV
			if((base->flagMPFV==1)&&(base->flagTurb==1))
			{
				if(base->VARTAG[count-1]==41)
				{
					count1=find_nv(base,31);
					((CgnsZone->ZoneBc+BC_index)->Data_Set+(count1))->Bc_Type=c;
				}
				if(base->VARTAG[count-1]==42)
				{
					count1=find_nv(base,32);	
					((CgnsZone->ZoneBc+BC_index)->Data_Set+(count1))->Bc_Type=c;
				}
				if(base->VARTAG[count-1]==43)
				{
					count1=find_nv(base,33);
					((CgnsZone->ZoneBc+BC_index)->Data_Set+(count1))->Bc_Type=c;
				}
			}	
			//End--Chandrasekhar
			
			//Added by adhi+shash
			if((base->flagParticle==1)||(base->flagDroplet==1))
			{
				Find_Direction(&((zbc+BC_index)->Direction),(zbc+BC_index)->Point_Range);
				if ((zbc+BC_index)->Direction==1||(zbc+BC_index)->Direction==2)
					size=(nk-1)*(nj-1);
				if ((zbc+BC_index)->Direction==3||(zbc+BC_index)->Direction==4)
					size=(nk-1)*(ni-1);
				if ((zbc+BC_index)->Direction==5||(zbc+BC_index)->Direction==6)
					size=(nj-1)*(ni-1);
				
				domain = MakeDomainArray(CgnsZone,nv);
				
				slab   = FindSlabPtrBC((zbc+BC_index),nv);
				if((nv==62)||(nv==63)||(nv==64)||(nv==65)||(nv==710))//shash //||(nv==67)-- number density equation is no longer sovled for
				{
					printf("\nI am entering boundary_d_t_s in InputBoundarycondition and nv=%d",nv);
					boundary_d_t_s(base,(zbc+BC_index)->Direction,domain,nv,CgnsZone,size,((zbc+BC_index)->Data_Set+find_nv(base,nv)),slab);
				}
			}
			ts=0;
		}
		
		printf("\n");
		for(nv=0;nv<NS;nv++)
			printf("nv = %d\tName of Variable = %s\n", nv,((CgnsZone->ZoneBc+BC_index)->Data_Set+nv)->Name_of_variable);
	}
	
	if((base->flagMPFV==1)&&(base->flagTurb==0))//changed from  base->flagMPFV==1 to (base->flagMPFV==1)&&(base->flagTurb==0) to account only for laminar MPFV--chandrasekhar
	{
		for(BC_index=0;BC_index<(CgnsZone->Number_of_BC_Patches);BC_index++)
			update_boundaryDMPFV(base,CgnsZone,zbc+BC_index);
	}

//The following lines are commented as the momentum equations of dispersed-phase are being solved for velocity components and not momentum components
/*	if((base->flagParticle==1)||(base->flagDroplet==1))
	{
		for(BC_index=0;BC_index<(CgnsZone->Number_of_BC_Patches);BC_index++)
			update_boundaryDisperse(base,CgnsZone,zbc+BC_index);
	}	
*/	
	return (0);
}
