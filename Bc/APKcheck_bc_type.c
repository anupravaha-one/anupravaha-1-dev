#include "APKcheck_bc_type.h"

/**
The present function check_bc_type() checks whether the  BC type entered by the user is one of the option 
provided by the solver or not it gives an error message in the later case if the conditions do not satisfyies. 
*/


FILE *var2,*var3;
int check_bc_type(Base *base,ZoneBoundary *zBc,int *bc,int vartag)
{
	int test=0;
	int c;
	
        var2=fopen("Settings/coefficients.dat","a");	
        var3=fopen("Settings/ProgramReadInputValues.dat","a");
        do
	{
		if(base->flagTurbModel == 1 || base->flagTurbModel == 3)
		{
			
			if(base->flagMPFV==0)//Turbulent single phase flow--chandrasekhar
			{
				if(zBc->Boundary_Tag == 22)
				{
					//printf("entered inside of the loop vartag=%d\n",vartag);
					//printf("I am entering the Boundary tag 22 in single phase flow(Turbulent) loop in the check_bc_type\n");
					switch(vartag)
					{
						case 1:
						case 2:
						case 3:
							scanf("%d", &c);        
                                                        fprintf(var2,"%d\n",c);        
                                                        fprintf(var3,"%d \t \n\n",c);
							//printf("scanned value is %d\n",c);
							if(c == 2 || c == 3 || c == 7)
								test = 1;
							else
							{
								printf("Only Dirichlet Boundary Condition is allowed ... \n");
								//printf("\nENTER \n 2 for Dirichlet..\n 3 for Homogenous Dirichlet\n 7 for TimeDependent \n");
								printf("\nFor Dirichlet Boundary Conditions...............[Enter 2]");
								printf("\nFor Homogenous Dirichlet Boundary Conditions....[Enter 3]");
								printf("\nFor TimeDependent Boundary Conditions...........[Enter 7]\n");
								test = 0;
							}
							break;
						case 6:
							scanf("%d", &c);                                                        
                                                        fprintf(var2,"%d\n",c);        
                                                        fprintf(var3,"%d \t \n\n",c);
							//printf("scanned value is %d\n",c);
							if(c==1 || c == 2 || c == 3 || c == 7)
								test = 1;
							else
							{
								printf("\nINCORRECT INPUT\n");
								//printf("\nENTER \n1 for  Neumann\n 2 for Dirichlet\n 3 for Homogenous Dirichlet\n 7 for TimeDependent \n");
								printf("\nFor Neumann Boundary Conditions.................[Enter 1]");
								printf("\nFor Dirichlet Boundary Conditions...............[Enter 2]");
								printf("\nFor Homogenous Dirichlet Boundary Conditions....[Enter 3]");
								printf("\nFor TimeDependent Boundary Conditions...........[Enter 7]\n");
								printf("\nPLEASE REFER THE ABOVE OPTIONS AND ENTER THE RIGHT CHOICE\n");
								test=0;
							}
							break;
						case 600:
							scanf("%d", &c);                                                        
                                                        fprintf(var2,"%d\n",c);        
                                                        fprintf(var3,"%d \t \n\n",c);
							if(c == 1)
								test = 1;
							else
							{
								printf("Only Homogeneous Neumann (Zero-gradient) Boundary Condition is allowed ... \n");
								//printf("\nENTER  1 for  Neumann\n");
								printf("\nFor Neumann Boundary Conditions.................[Enter 1]\n");
								test = 0;
							}
							break;
						case 601:
							c = 3;
							test = 1;
							printf("\nDIRICHLET BOUNDARY CONDITION is applied for Epsilon.(Ref:See Wall Function Treatment in Anupravaha Report)\n");
							//printf("\nHOMOGENEOUS DIRICHLET BOUNDARY CONDITION is applied for Epsilon.(Ref:See Wall Function Treatment in Project Report)\n");
							break;
						case 14://Shyam
							scanf("%d", &c);                                                        
                                                        fprintf(var2,"%d\n",c);        
                                                        fprintf(var3,"%d \t \n\n",c);
							if(c == 1)
								test = 1;
							else
							{
								printf("Only Homogeneous Neumann (Zero-gradient) Boundary Condition is allowed ... \n");
								//printf("\nENTER  1 for  Neumann\n");
								printf("\nFor Neumann Boundary Conditions.................[Enter 1]\n");
								test = 0;
							}
							break;
						case 15://Shyam
							scanf("%d", &c);                                                        
                                                        fprintf(var2,"%d\n",c);        
                                                        fprintf(var3,"%d \t \n\n",c);
							if(c == 1)
								test = 1;
							else
							{
								printf("Only Homogeneous Neumann (Zero-gradient) Boundary Condition is allowed ... \n");
								//printf("\nENTER  1 for  Neumann\n");
								printf("\nFor Neumann Boundary Conditions.................[Enter 1]\n");
								test = 0;
							}
							break;
						case 7: //Shyam
							scanf("%d", &c);                                                        
                                                        fprintf(var2,"%d\n",c);        
                                                        fprintf(var3,"%d \t \n\n",c);
							if(c == 1&&base->flagTurbComb==1) //Shyam, for enthalpy wall treatment is not there, only insulated wall BC is applicable at present
								test = 1;
							else
							{
								printf("Only Homogeneous Neumann (Zero-gradient) Boundary Condition is allowed ... \n");
								//printf("\nENTER  1 for  Neumann\n");
								printf("\nFor Neumann Boundary Conditions.................[Enter 1]\n");
								test = 0;
							}
							break;
						case 4:
							scanf("%d", &c);                                                        
                                                        fprintf(var2,"%d\n",c);        
                                                        fprintf(var3,"%d \t \n\n",c);
							//printf("scanned value is %d\n",c);
							if(c==0||c==1||c==2||c==3||c==4||c==5||c==7)
								test=1;
							else
							{
								test=0;
								printf("\nINCORRECT INPUT\n");
								printf("\nPLEASE REFER THE ABOVE OPTIONS AND ENTER THE RIGHT CHOICE\n");
							}
							break;
						default:
							scanf("%d", &c);                                                        
                                                        fprintf(var2,"%d\n",c);        
                                                        fprintf(var3,"%d \t \n\n",c);
							if(c==0||c==1||c==2||c==3||c==4||c==5||c==7)
								test=1;
							else
							{
								test=0;
								printf("\nINCORRECT INPUT\n");
								printf("\nPLEASE REFER THE ABOVE OPTIONS AND ENTER THE RIGHT CHOICE\n");
							}
							//printf("Should not be there in check_bctype\n exiting......\n");
							break;
					}
				}
				else 
				{
					//printf("I am entering the general turbulent flag for single phase flow\n");
					scanf("%d", &c);                                                        
                                        fprintf(var2,"%d\n",c);        
                                        fprintf(var3,"%d \t \n\n",c);
					printf("scanned value is %d\n",c);
					if(c==0||c==1||c==2||c==3||c==4||c==5||c==7)
						test=1;
					else
					{
						test=0;
						printf("\nINCORRECT INPUT\n");
						printf("\nPLEASE REFER THE ABOVE OPTIONS AND ENTER THE RIGHT CHOICE\n");
					}
				}
			}
			else if(base->flagMPFV==1)//Turbulent multiphase flow(bubbly)--chandrasekhar
			{
				if(zBc->Boundary_Tag == 22)
				{
					//printf("entered inside of the loop vartag=%d\n",vartag);
					printf("I am entering the Boundary tag 22 in MPFV(Multi phase turbulent) loop in the check_bc_type\n");
					switch(vartag)
					{
						case 31:
						case 32:
						case 33:
						case 41:
						case 42:
						case 43:
							scanf("%d", &c);                                                        
                                                        fprintf(var2,"%d\n",c);        
                                                        fprintf(var3,"%d \t \n\n",c);
							//printf("scanned value is %d\n",c);
							if(c == 2 || c == 3)
								test = 1;
							else
							{
								printf("Only Dirichlet Boundary Condition is allowed ... \n");
								//printf("\nENTER \n 2 for Dirichlet..\n 3 for Homogenous Dirichlet\n ");
								printf("\nFor Dirichlet Boundary Conditions...............[Enter 2]");
								printf("\nFor Homogenous Dirichlet Boundary Conditions....[Enter 3]\n");
								test = 0;
							}
							break;
						case 44:
						case 45:
						case 46:
							scanf("%d", &c);                                                        
                                                        fprintf(var2,"%d\n",c);        
                                                        fprintf(var3,"%d \t \n\n",c);
							if(c == 2 || c == 3)
								test = 1;
							else
							{
								printf("Only Dirichlet Boundary Condition is allowed ... \n");
								//printf("\nENTER \n 2 for Dirichlet..\n 3 for Homogenous Dirichlet\n ");
								printf("\nFor Dirichlet Boundary Conditions...............[Enter 2]");
								printf("\nFor Homogenous Dirichlet Boundary Conditions....[Enter 3]\n");
								test = 0;
							}
							break;
						case 600:
							scanf("%d", &c);                                                        
                                                        fprintf(var2,"%d\n",c);        
                                                        fprintf(var3,"%d \t \n\n",c);
							if(c == 1)
								test = 1;
							else
							{
								printf("Only Homogeneous Neumann (Zero-gradient) Boundary Condition is allowed ... \n");
								//printf("\nENTER  1 for  Neumann\n");
								printf("\nFor Neumann Boundary Conditions.................[Enter 1]\n");
								test = 0;
							}
							break;
						case 601:
							c = 3;
							test = 1;
							printf("\nHOMOGENEOUS DIRICHLET BOUNDARY CONDITION is applied for Epsilon.(Ref:See Wall Function Treatment in Project Report)\n");
							break;
							
						case 4:
							scanf("%d", &c);                                                        
                                                        fprintf(var2,"%d\n",c);        
                                                        fprintf(var3,"%d \t \n\n",c);
							//printf("scanned value is %d\n",c);
							if(c==0||c==1||c==2||c==3)
								test=1;
							else
							{
								test=0;
								printf("\nINCORRECT INPUT\n");
								printf("\nPLEASE REFER THE ABOVE OPTIONS AND ENTER THE RIGHT CHOICE\n");
							}
							break;
						default:
							scanf("%d", &c);                                                        
                                                        fprintf(var2,"%d\n",c);        
                                                        fprintf(var3,"%d \t \n\n",c);
							if(c==0||c==1||c==2||c==3)
								test=1;
							else
							{
								test=0;
								printf("\nINCORRECT INPUT\n");
								printf("\nPLEASE REFER THE ABOVE OPTIONS AND ENTER THE RIGHT CHOICE\n");
							}
							//printf("Should not be there in check_bctype\n exiting......\n");
							break;
					}
				}
				else 
				{
					printf("I am entering the general turbulent flag in MPFV\n");
					scanf("%d", &c);                                                        
                                        fprintf(var2,"%d\n",c);        
                                        fprintf(var3,"%d \t \n\n",c);
					printf("scanned value is %d\n",c);
					if(c==0||c==1||c==2||c==3)
						test=1;
					else
					{
						test=0;
						printf("\nINCORRECT INPUT\n");
						printf("\nPLEASE REFER THE ABOVE OPTIONS AND ENTER THE RIGHT CHOICE\n");
					}
				}
			}
				
		}
		//Added by sushil for multiphase flow
		
		else if((base->flagMPFV==1)&&(base->flagTurb==0))//changed from (base->flagMPFV==1) to (base->flagMPFV==1)&&(base->flagTurb==0) to acc. only for laminar MPFV--chandrasekhar
		{	
			printf("Am entering the MPFV(Laminar) in the check_bc_type\n");
			scanf("%d", &c);                                                        
                        fprintf(var2,"%d\n",c);        
                        fprintf(var3,"%d \t \n\n",c);
			printf("scanned value is %d\n",c);
			if(c==0||c==1||c==2||c==3)
			test=1;
			else
			{
				test=0;
				printf("\nINCORRECT INPUT\n");
				printf("\nPLEASE REFER THE ABOVE OPTIONS AND ENTER THE RIGHT CHOICE\n");
			}
		
		}
		
		//END
		else 
		{	
			//printf("I am entering the general(Laminar single phase flow) case in check_bc_type\n");
			scanf("%d", &c);                                                        
                        fprintf(var2,"%d\n",c);        
                        fprintf(var3,"%d \t \n\n",c);
			printf("scanned value is %d\n",c);
			if(c==0||c==1||c==2||c==3||c==4||c==5||c==7||c==9||c==90)
				test=1;
			else
			{
				test=0;
				printf("\nINCORRECT INPUT\n");
				printf("\nPLEASE REFER THE ABOVE OPTIONS AND ENTER THE RIGHT CHOICE\n");
			}
		}
		printf("test=%d\n",test);
	}while(test!=1);

	*bc = c;
        fclose(var2);
        fclose(var3);
return(test);
}
