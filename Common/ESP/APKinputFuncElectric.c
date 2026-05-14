#include "APKinputFuncElectric.h"
/**
This function, inputFuncElectric() is for deciding different options available in the solver for Electric field. The selection 
of options decides type of coupling, type of scheme, type of solver etc. that well suits your problem. 
*/
int inputFuncElectric(Root CgnsRoot)
{
        FILE *var2,*var3;
	int scheme=0,flagT=0,Etag=0;
        
        var2=fopen("Settings/coefficients.dat","a");	
        var3=fopen("Settings/ProgramReadInputValues.dat","a");
	do
	{
		//printf("\nEnter whether you want solve an Electric Field problem[Enter->1] or No Electric Field[Enter->0]");
		//scanf("%d",&flagT);
		flagT = 0;
	}while(!(flagT == 1 || flagT == 0));                
        //fprintf(var2,"%d\n",flagT);        
        //fprintf(var3,"%d \tEnter whether you want solve an Electric Field problem[Enter->1] or No Electric Field[Enter->0] \n\n",flagT);
	CgnsRoot.CgnsBase->flagEfield=flagT;
	
	if(CgnsRoot.CgnsBase->flagEfield==1)
	{
		flagT = 0;
		do
		{
			printf("\nEnter whether you want solve an Electrically coupled Fluid Flow problem[Enter->1] or only Electric Field[Enter->2]");
			scanf("%d",&flagT);
		}while(!(flagT == 1 || flagT == 2));        
                fprintf(var2,"%d\n",flagT);        
                fprintf(var3,"%d \tEnter whether you want solve an Electrically coupled Fluid Flow problem[Enter->1] or only Electric Field[Enter->2] \n\n",flagT);
		CgnsRoot.CgnsBase->flagEfield=flagT;
		
		Etag = flagT; //Implies we are solving for electric field
		
		if(CgnsRoot.CgnsBase->flagEfield==1)
		{
			flagT = 0;
			do
			{
				printf("\nPlease Enter the type of coupling between Flow and Electric Fields:-\nOne-Way Coupling->1:\nTwo-Way Coupling->2:");
				scanf("%d",&flagT);
			}while(!(flagT == 1 || flagT == 2));        
                        fprintf(var2,"%d\n",flagT);        
                        fprintf(var3,"%d \tPlease Enter the type of coupling between Flow and Electric Fields:One-Way Coupling->1: Two-Way Coupling->2: \n\n",flagT);
			CgnsRoot.CgnsBase->FtoEcouplingType=flagT;
		}
		
		printf("\nThe following information is specifically intended for ELECTRIC FIELD ONLY\n");
		
		do
		{
			printf("\nPlease Enter the solution method for ELECTRIC FIELD Equations:-\nCOUPLED->0:\nSEMICOUPLED->1:\nSEGREGATED->2:");
			scanf("%d",&scheme);
			CgnsRoot.CgnsBase->EcoupTag=scheme; 
		}while(!(scheme == 0 || scheme == 1 || scheme == 2));        
                fprintf(var2,"%d\n",scheme);        
                fprintf(var3,"%d \tPlease Enter the solution method for ELECTRIC FIELD Equations:-COUPLED->0: SEMICOUPLED->1: SEGREGATED->2: \n\n",scheme);
		
		if(scheme==0 || scheme==1)
		{
			printf("\nEnter Solver Scheme for ELECTRIC FIELD Equations:\n");
			printf("If Implicit-0\n");
			printf("If SemiImplicit-1\n");
			scanf("%d",&scheme);        
                        fprintf(var2,"%d\n",scheme);        
                        fprintf(var3,"%d \tSolver Scheme Implicit-0 SemiImplicit-1 \n\n",scheme);
			CgnsRoot.CgnsBase->EsolScheme = scheme; 
		}
		else
		{
			printf("\nEnter Solver Scheme for ELECTRIC FIELD Equations:\n");
			printf("If Implicit-0\n");
			printf("If SemiImplicit-1\n");
			printf("If PseudoCompressiblity-\n");
			printf("If (a)one T and one Tau step-21 \n");
			printf("If (a)Impose Continuity at each step-22\n");
			scanf("%d",&scheme);        
                        fprintf(var2,"%d\n",scheme);        
                        fprintf(var3,"%d \tSolver Scheme Implicit-0 SemiImplicit-1 PseudoCompressiblity (a)one T and one Tau step-21 (a)Impose Continuity at each step-22\n\n",scheme);
			CgnsRoot.CgnsBase->EsolScheme = scheme; 
		}
		
		
	}
	else
		Etag = 0; //Implies we are not solving for electric field
	
        fclose(var2);
        fclose(var3);
	return(Etag);
}

