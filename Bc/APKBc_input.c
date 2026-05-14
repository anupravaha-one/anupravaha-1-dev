#include "APKBc_input.h"

/**
This function Bc_input() acts as an interface between the various modules of the boundary condition entry 
and saving, recovery and editing of the entered boundary conditions. This module is a calling function
which calls the various other modules of the boundary condition routine according to choice of    
user for example `InputBoundarycondition' for entering boundary condition related information 
or `edit_bc' for editing some of the information entered or `Save_bc' for saving the boundary    
condition entered by the user or maybe for recovering the boundary condition entered 

This function first of calls initialize() function to initialize all the variables.
Bc_input() offers three modes of BC input for the user 
1) For all fresh input for boundary conditions, InputBoundarycondition() is called for this purpuse 
   functsave() is then called to record these boundary conditions so that they could be used later
	
2) Solving with previously stored Boundary conditions. Retrive_bc() is called for this purpose.
     
3) Solving with editing some of the  previously stored boundary conditions. edit_bc() along with Retrive_bc is 
called for this purpose the freshly edited BC is stored again by calling functsave(). 
*/

int Bc_input(Root CgnsRoot)
{
	FILE *bc,*var2,*var3;
	FILE *chfile;
	int base_ind,zone_ind,No_of_variable=0,bc_index,nv,Variables=0,RadiationVariable=0;
	int choice,ai,bi,ci,di,bci;	
        
	/*The looping over all the bases starts here*/
	for(base_ind=0;base_ind<CgnsRoot.Number_of_bases;base_ind++)
	{
		//The following  module initialises the variables to be solved and allocates them sufficient memory
		if((CgnsRoot.CgnsBase+base_ind)->BaseParameters->flagRadiation==0)
		{	
			
			if((CgnsRoot.CgnsBase+base_ind)->flagGroundwater==1)
			{
				No_of_variable=initializeGroundwater(CgnsRoot.CgnsBase+base_ind);
			}
			else
			{
				No_of_variable=initialize(CgnsRoot.CgnsBase+base_ind);
			}
			
			if((CgnsRoot.CgnsBase+base_ind)->BaseParameters->flagflowRadiation==1)
			{
				printf("\nGathering information for Radiation Variable");
				printf("\nDo u want to specify Radiative Boundary Condition in terms of Radiative flux or Temperature");
				printf("\n0->RadiativeFlux-1->Temperature\n");
				scanf("%d",&(CgnsRoot.CgnsBase->BaseParameters->flagRadiativeTemperature));                             
                                var2=fopen("Settings/coefficients.dat","a");        
                                var3=fopen("Settings/ProgramReadInputValues.dat","a");
                                fprintf(var2,"%d\n",CgnsRoot.CgnsBase->BaseParameters->flagRadiativeTemperature);	
                                fprintf(var3,"%d \tDo u want to specify Radiative Boundary Condition in terms of Radiative flux or Temperature 0->RadiativeFlux-1->Temperature\n\n",CgnsRoot.CgnsBase->BaseParameters->flagRadiativeTemperature);
                                fclose(var2);
                                fclose(var3);
			
				Variables=initializeRadiation(CgnsRoot.CgnsBase+base_ind, No_of_variable);
				RadiationVariable=Variables-No_of_variable;
				CgnsRoot.CgnsBase->Radvar = RadiationVariable;
			}
			else
			CgnsRoot.CgnsBase->Radvar = 0;
			
			//printf("\nIf you want to solve with ENTIRELY NEW BOUNDARY CONDITIONS:ENTER...'0' : ");
			/*printf("\nIF you want to solve with OLD BOUNDARY CONDITIONS:ENTER...'1'\n");
			printf("\nIF you want to solve with OLD BOUNDARY CONDITIONS while EDITING some of them :ENTER...'2'\n");*/
		
			//scanf("%d",&choice);           
			choice = 0;           
                        //var2=fopen("Settings/coefficients.dat","a");   	
                        //var3=fopen("Settings/ProgramReadInputValues.dat","a");  
                        //fprintf(var2,"%d\n",choice);        
                        //fprintf(var3,"%d \tIf you want to solve with ENTIRELY NEW BOUNDARY CONDITIONS:ENTER...'0' \n\n",choice);
                        //fclose(var2);
                        //fclose(var3);
			/*if (choice==1||choice==2)
			{
	// 			This file lists all the previously saved boundary condition
				//system("gvim Settings/file_bc.dat");
			}*/
			if(choice==0)
			{
				/*The boundary conditions entered by the user are saved in this binary file*/
				bc=fopen("Settings/BC.dat","w");
				/*The listing of the previously entered Boundary conditions are done in this file*/
				chfile=fopen("Settings/file_bc.dat","w");
				/*The looping over all zones starts here*/
				for(zone_ind=0;zone_ind<(CgnsRoot.CgnsBase+base_ind)->Number_of_zones;zone_ind++)
				{
					fprintf(chfile,"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
					fprintf(chfile,"\nZone Name=%s\n",((CgnsRoot.CgnsBase+base_ind)->CgnsZone+zone_ind)->zonename);
	
					//The following module does acts as a basic interface for entry of Boundary condition
					/*if(CgnsRoot.CgnsBase->flagTurbModel == 1 || CgnsRoot.CgnsBase->flagTurbModel == 3)	//Standard and RNG k-epsilon model
						InputBoundaryconditionTurb(CgnsRoot.CgnsBase,No_of_variable,(CgnsRoot.CgnsBase+base_ind)->CgnsZone+zone_ind,(CgnsRoot.CgnsBase+base_ind)->BaseParameters->flag2d);
					else*/
					InputBoundarycondition(CgnsRoot.CgnsBase,No_of_variable,(CgnsRoot.CgnsBase+base_ind)->CgnsZone+zone_ind,(CgnsRoot.CgnsBase+base_ind)->BaseParameters->flag2d);
					
					if((CgnsRoot.CgnsBase+base_ind)->BaseParameters->flagflowRadiation==1)
					InputRadiationBoundarycondition(CgnsRoot.CgnsBase+base_ind,No_of_variable,RadiationVariable,(CgnsRoot.CgnsBase+base_ind)->CgnsZone+zone_ind,(CgnsRoot.CgnsBase+base_ind)->BaseParameters->flag2d);
					/*This module saves all the boundary condition entered*/
					//functsave(bc,No_of_variable,(CgnsRoot.CgnsBase+base_ind)->CgnsZone+zone_ind);
						
					for(bc_index=0;bc_index<((CgnsRoot.CgnsBase+base_ind)->CgnsZone+zone_ind)->Number_of_BC_Patches;bc_index++)
					{
						fprintf(chfile,"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
						fprintf(chfile,"\nBcPatch name=%s\n",(((CgnsRoot.CgnsBase+base_ind)->CgnsZone+zone_ind)->ZoneBc+bc_index)->BCName);	
						for(nv=0;nv<No_of_variable;nv++)
						{  
	 					//	printf("\n\nnv=%d\tvalue=%f",nv,*(((((CgnsRoot.CgnsBase+base_ind)->CgnsZone+zone_ind)->ZoneBc+bc_index)->Data_Set+nv)->Data_Array));//extra by sushil
							fprintf(chfile,"\nVariable=%c,BcType=%d\n",((((CgnsRoot.CgnsBase+base_ind)->CgnsZone+zone_ind)->ZoneBc+bc_index)->Data_Set+nv)->Name_of_variable[0],((((CgnsRoot.CgnsBase+base_ind)->CgnsZone+zone_ind)->ZoneBc+bc_index)->Data_Set+nv)->Bc_Type);
							ai=*(((((CgnsRoot.CgnsBase+base_ind)->CgnsZone+zone_ind)->ZoneBc+bc_index)->Data_Set+nv)->identity);
							bi=*(((((CgnsRoot.CgnsBase+base_ind)->CgnsZone+zone_ind)->ZoneBc+bc_index)->Data_Set+nv)->identity+1);
							ci=*(((((CgnsRoot.CgnsBase+base_ind)->CgnsZone+zone_ind)->ZoneBc+bc_index)->Data_Set+nv)->identity+2);
							di=*(((((CgnsRoot.CgnsBase+base_ind)->CgnsZone+zone_ind)->ZoneBc+bc_index)->Data_Set+nv)->identity+3);
							bci=((((CgnsRoot.CgnsBase+base_ind)->CgnsZone+zone_ind)->ZoneBc+bc_index)->Data_Set+nv)->Bc_Type;
							if(bci!=4 || bci!=9 || bci!=90)
							{
								if(ai==0)
								{	
								
									
	
									
									fprintf(chfile,"\nUser entered =%lf\n",*(((((CgnsRoot.CgnsBase+base_ind)->CgnsZone+zone_ind)->ZoneBc+bc_index)->Data_Set+nv)->Data_Array));
								}
                                                               
								else
								{
									fprintf(chfile,"\nUser entered Function\n");
								}
							}
							else if(bci!=9 || bci!=90)
							{
								if(bi==0)
								{
									fprintf(chfile,"\nUser entered Alpha value=%lf\n",*(((((CgnsRoot.CgnsBase+base_ind)->CgnsZone+zone_ind)->ZoneBc+bc_index)->Data_Set+nv)->Data_Array));
								}
								else
								{
									fprintf(chfile,"\nUser entered a Function for alpha\n");
								}
								if(ci==0)
								{
									fprintf(chfile,"\nUser entered Beta value=%lf\n",*(((((CgnsRoot.CgnsBase+base_ind)->CgnsZone+zone_ind)->ZoneBc+bc_index)->Data_Set+nv)->Data_Array+1));
								}
								else
								{
									fprintf(chfile,"\nUser entered a Function for Beta\n");
								}
								if(di==0)
								{
									fprintf(chfile,"\nUser entered Gamma value=%lf\n",*(((((CgnsRoot.CgnsBase+base_ind)->CgnsZone+zone_ind)->ZoneBc+bc_index)->Data_Set+nv)->Data_Array+2));
								}
								else
								{
									fprintf(chfile,"\nUser entered a Function for Gamma\n");
								}
								
								
							}
							if(bci==9 || bci==90)
							  {
                                                          	fprintf(chfile,"\nUser entered zero efflux condition\n");
                                                          }
						}
						fprintf(chfile,"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
					}	
					fprintf(chfile,"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
				}
				fclose(bc);
				fclose(chfile);
			}
		}
		else if((CgnsRoot.CgnsBase+base_ind)->BaseParameters->flagRadiation==1)
		{
			printf("\nGathering information for Radiation Variable");
			printf("\nDo u want to specify Radiative Boundary Condition in terms of Radiative flux or Temperature");
			printf("\n0->RadiativeFlux-1->Temperature\n");
			scanf("%d",&(CgnsRoot.CgnsBase->BaseParameters->flagRadiativeTemperature));                            
                        var2=fopen("Settings/coefficients.dat","a");      	
                        var3=fopen("Settings/ProgramReadInputValues.dat","a");                           
                        fprintf(var2,"%d\n",CgnsRoot.CgnsBase->BaseParameters->flagRadiativeTemperature);         
                        fprintf(var3,"%d \tDo u want to specify Radiative Boundary Condition in terms of Radiative flux or Temperature 0->RadiativeFlux-1->Temperature\n\n",CgnsRoot.CgnsBase->BaseParameters->flagRadiativeTemperature);       
                        fclose(var2);
                        fclose(var3);
	
			Variables=initializeRadiation(CgnsRoot.CgnsBase+base_ind, No_of_variable);

			RadiationVariable=Variables-No_of_variable;
			CgnsRoot.CgnsBase->Radvar = RadiationVariable;
			for(zone_ind=0;zone_ind<(CgnsRoot.CgnsBase+base_ind)->Number_of_zones;zone_ind++)
			InputRadiationBoundarycondition(CgnsRoot.CgnsBase+base_ind,No_of_variable,RadiationVariable,(CgnsRoot.CgnsBase+base_ind)->CgnsZone+zone_ind,(CgnsRoot.CgnsBase+base_ind)->BaseParameters->flag2d);

		}
	}
	
	return (No_of_variable+RadiationVariable);
}
