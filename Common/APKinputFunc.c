#include "APKinputFunc.h"
/**
This function, inputFunc() is for deciding different options available in the solver. The selection of options decides type of flow,
type of scheme, type of solver etc. that well suits your problem. 
*/

int inputFunc(Root CgnsRoot)
{
        FILE *var2,*var3;
        
	int open,scheme,restart,varRho,flagT,flagLogKEps=0;
	int flagsoli,flagTP;//bhans..anji
	int densityvar,flagGr,flagRadiative,nt,np,numzones,iZone;
	int flagLMMHD;// naveen
 	int flaginduced;//rajendra
	//int flagMPFV;//sushil//shash+adhi
	double C,number;
	Zone *zone;
	int flagMultiPhase,Multiflowtype,flagevaporation,flagdropfuelmassfrac,constant,Dragflag;//shash+adhi
	int flagDispersion;
        int flagWf;
	numzones = CgnsRoot.CgnsBase->Number_of_zones;
        
	var2=fopen("Settings/coefficients.dat","a");
	if (var2==NULL){printf("File opening if failed in APKinputFunc.c\nCan not open Settings/coefficients.dat to append...........\n");}
        fclose(var2);
        
//	printf("Open the Summary File (1 yes/0 no):");
//	scanf("%d",&open);
//	if(open==1)	
//	system("gvim Settings/Summary &");

	do
	{
		//printf("\nDo you want to solve only Radiative Model [1] or Flow solution with or without Radiation [0] : ");
		//scanf("%d",&flagRadiative);
		flagRadiative = 0;
	}while(!(flagRadiative == 1 || flagRadiative == 0));
		
        //var2=fopen("Settings/coefficients.dat","a");
	//var3=fopen("Settings/ProgramReadInputValues.dat","a");
        //fprintf(var2,"%d\n",flagRadiative);
        //fprintf(var3,"%d \tDo you want to solve only Radiative Model [1] or Flow solution with or without Radiation [0] \n\n",flagRadiative);
        //fclose(var2);
        //fclose(var3);
	
	if(flagRadiative==0){
		do
		{
			//printf("\nFlow Solution with Radiation [1] or without Radiation [0] : ");
			//scanf("%d",&(CgnsRoot.CgnsBase->BaseParameters->flagflowRadiation));
			CgnsRoot.CgnsBase->BaseParameters->flagflowRadiation = 0;
		}while(!(CgnsRoot.CgnsBase->BaseParameters->flagflowRadiation==1||CgnsRoot.CgnsBase->BaseParameters->flagflowRadiation==0));
                	
                //var2=fopen("Settings/coefficients.dat","a");
                //var3=fopen("Settings/ProgramReadInputValues.dat","a");
                //fprintf(var2,"%d\n",CgnsRoot.CgnsBase->BaseParameters->flagflowRadiation);
                //fprintf(var3,"%d \tFlow Solution with Radiation [1] or without Radiation [0] \n\n",CgnsRoot.CgnsBase->BaseParameters->flagflowRadiation);
                //fclose(var2);
                //fclose(var3);
        }
        
	if(flagRadiative==1||CgnsRoot.CgnsBase->BaseParameters->flagflowRadiation==1)
		{
		printf("\nEnter The Number of Discrete Cells in Theta Direction : ");
		scanf("%d",&nt);	
                
                var2=fopen("Settings/coefficients.dat","a");
                var3=fopen("Settings/ProgramReadInputValues.dat","a");
                fprintf(var2,"%d\n",nt);
                fprintf(var3,"%d \tEnter The Number of Discrete Cells in Theta Direction \n\n",nt);
                fclose(var2);
                fclose(var3);
        
		printf("\nEnter The Number of Discrete Cells in Phi Direction : ");
		scanf("%d",&np);	
                	
                var2=fopen("Settings/coefficients.dat","a");
                var3=fopen("Settings/ProgramReadInputValues.dat","a");
                fprintf(var2,"%d\n",np);
                fprintf(var3,"%d \tEnter The Number of Discrete Cells in Phi Direction \n\n",np);
                fclose(var2);
                fclose(var3);
        
		CgnsRoot.CgnsBase->BaseParameters->ntheta=nt;
		CgnsRoot.CgnsBase->BaseParameters->nphi=np;
		}
		
	CgnsRoot.CgnsBase->BaseParameters->flagRadiation=flagRadiative;
	
	if(flagRadiative==0||CgnsRoot.CgnsBase->BaseParameters->flagflowRadiation==1)
	{
	//Call the electric field settings
	flagT = inputFuncElectric(CgnsRoot); //Returns the tag whether we are solving for electric field or not	
	
	do
	{
		printf("\nEnter whether you want to solve a ");
		printf("\nTurbulent Flow problem----[Enter->1] ");
		printf("\nLaminar flow problem------[Enter->0] : ");
		scanf("%d",&flagT);
	}while(!(flagT == 1 || flagT == 0));	
        	
        var2=fopen("Settings/coefficients.dat","a");
	var3=fopen("Settings/ProgramReadInputValues.dat","a");
        fprintf(var2,"%d\n",flagT);
        fprintf(var3,"%d \tEnter whether you want to solve a Turbulent Flow problem[Enter->1] or Laminar flow problem[Enter->0] \n\n",flagT);
        fclose(var2);
        fclose(var3);
	
	CgnsRoot.CgnsBase->BaseParameters->Turb = (ParametersTurb *)calloc(1,sizeof(ParametersTurb));
	CgnsRoot.CgnsBase->BaseParameters->Turb->numTurbScalars = 0;
	CgnsRoot.CgnsBase->BaseParameters->Turb->numDualTimeIter = 1;//Do not tamper with this ... Very important statement <<< Narsimha.
	CgnsRoot.CgnsBase->BaseParameters->Turb->flagLimKEps = 0;
	CgnsRoot.CgnsBase->BaseParameters->Turb->flagLimTVisc = 0;
	CgnsRoot.CgnsBase->BaseParameters->Turb->flagLimTimeScale = 0;
	CgnsRoot.CgnsBase->BaseParameters->Turb->flagLimProduction = 0;
	CgnsRoot.CgnsBase->flagTurb=flagT;
	CgnsRoot.CgnsBase->flagTurbModel=flagT;
	CgnsRoot.CgnsBase->flagLogKEps = 0;	//This is neccessary for laminar flow problems which if not written causes problem in Dirichlet Boundary Condition Implementation ... Narsimha.
	CgnsRoot.CgnsBase->flagTurbComb = 0;
        CgnsRoot.CgnsBase->flagSolidification = 0;
	CgnsRoot.CgnsBase->flagTwoPhase = 0;
	CgnsRoot.CgnsBase->flagPreNonPre = 0;       //Shyam, Initializing all flags with zeroes to avoid taking some garbage values
	CgnsRoot.CgnsBase->flagDiffKinetCont = 0;
	CgnsRoot.CgnsBase->flagtypeVariableDensity = 0;
	CgnsRoot.CgnsBase->flagVariableDensityProblem = 0;
	CgnsRoot.CgnsBase->flagMPFV=0;
	CgnsRoot.CgnsBase->flagParticle=0;
	CgnsRoot.CgnsBase->flagDroplet=0;
	
	do
	{
		//printf("\nEnter whether you want to solve groundwater flow and contaminant transport (0 ->no / 1->yes): ");
		//scanf ("%d",&flagGr);
		flagGr = 0;
	}while(!(flagGr == 1 || flagGr == 0));	
                
        //var2=fopen("Settings/coefficients.dat","a");
	//var3=fopen("Settings/ProgramReadInputValues.dat","a");
        //fprintf(var2,"%d\n",flagGr);
        //fprintf(var3,"%d \tEnter whether you want to solve groundwater flow and contaminant transport (0 ->no / 1->yes): \n\n",flagGr);
        //fclose(var2);
        //fclose(var3);
        
	CgnsRoot.CgnsBase->flagGroundwater=flagGr;
	

	//added by sushil for MULTIPHSE FLOW VARIABLE(MPFV) problem.
	/*******************************************************************************************************************/		
		
		//adhi+shash
	do
	{
		//printf("\nSpecify whether you want to solve  multi phase  flow problem (0 ->no / 1->yes): ");
		//scanf ("%d",&flagMultiPhase);
		flagMultiPhase = 0;
	}while(!(flagMultiPhase == 1 || flagMultiPhase == 0));	
                
        //var2=fopen("Settings/coefficients.dat","a");
	//var3=fopen("Settings/ProgramReadInputValues.dat","a");
        //fprintf(var2,"%d\n",flagMultiPhase);
        //fprintf(var3,"%d \tSpecify whether you want to solve  multi phase  flow problem (0 ->no / 1->yes): \n\n",flagMultiPhase);
        //fclose(var2);
        //fclose(var3);
	
	if(flagMultiPhase == 1 )
	{
            do
            {
		printf("\nIf you are solving for bubbly flow enter 1; for Particles enter 2; for droplets enter 3: \n");
		scanf("%d",&Multiflowtype);
            }while(!(Multiflowtype==1 || Multiflowtype==2 || Multiflowtype==3));
            
            var2=fopen("Settings/coefficients.dat","a");
            var3=fopen("Settings/ProgramReadInputValues.dat","a");
            fprintf(var2,"%d\n",Multiflowtype);	
            fprintf(var3,"%d \tIf you are solving for bubbly flow enter 1; for Particles enter 2; for droplets enter 3: \n\n",Multiflowtype);
            fclose(var2);
            fclose(var3);

		if(Multiflowtype==1)
		{	
			CgnsRoot.CgnsBase->flagMPFV = 1 ; 
			printf("\nDo you want to include interphase momentum transfer term(i.e,Drag) If YES Enter-> 1  No Enter-> 0 : ");
			scanf("%d",&Dragflag);            
                        	
                        var2=fopen("Settings/coefficients.dat","a");
                        var3=fopen("Settings/ProgramReadInputValues.dat","a");
                        fprintf(var2,"%d\n",Dragflag);
                        fprintf(var3,"%d \tDo you want to include interphase momentum transfer term(i.e,Drag) If YES Enter-> 1  No Enter-> 0 \n\n",Dragflag);
                        fclose(var2);
                        fclose(var3);
			
			CgnsRoot.CgnsBase->BaseParameters->Dragflag=Dragflag;
			if(CgnsRoot.CgnsBase->BaseParameters->Dragflag==1)
			{	
				printf("\nFor the  drag calculation, If value of Drag coefficient(CD) is constant Enter 0\n\nOtherwise for user defined function for Drag coefficient(CD) Enter 1 : ");
				scanf("%d",&constant);            
                                	
                                var2=fopen("Settings/coefficients.dat","a");
                                var3=fopen("Settings/ProgramReadInputValues.dat","a");
                                fprintf(var2,"%d\n",constant);
                                fprintf(var3,"%d \tFor Drag coefficient(CD) is constant Enter 0, for user defined function Enter 1\n\n",constant);
                                fclose(var2);
                                fclose(var3);
        
				if(constant==0)
				{
					printf("\nEnter the constant  value of drag coefficinet(CD) : ");
					scanf("%lf",&C);            
                                        	
                                        var2=fopen("Settings/coefficients.dat","a");
                                        var3=fopen("Settings/ProgramReadInputValues.dat","a");
                                        fprintf(var2,"%lf\n",C);
                                        fprintf(var3,"%lf \tEnter the constant  value of drag coefficinet(CD) \n\n",C);
                                        fclose(var2);
                                        fclose(var3);
        
					CgnsRoot.CgnsBase->BaseParameters->ConstantCD = C;	
				}
					
				else
				{	
					printf("\n It will directly go to user defined function for CD calculation\n");
					printf("\nUser have freedom to change  user defined function in (user/APKuserDefinedVarCoefficients/userDefinedCD)\n\n");
					
					CgnsRoot.CgnsBase->constant = constant;
				}
				
				printf("\nDo you want to use normal formulation If YES  Enter 0\n\n If u want to use Small bubble approximation then Enter 1 : ");
				scanf("%d",&constant);            
                                	
                                var2=fopen("Settings/coefficients.dat","a");
                                var3=fopen("Settings/ProgramReadInputValues.dat","a");
                                fprintf(var2,"%d\n",constant);
                                fprintf(var3,"%d \tFor normal formulation  Enter 0, Small bubble approximation Enter 1 \n\n",constant);
                                fclose(var2);
                                fclose(var3);
        
				CgnsRoot.CgnsBase->BaseParameters->Bubbleapprox = constant;
			}
			//Added by chandrasekhar to include turbulent dispersion term in volume fraction equation
			if(CgnsRoot.CgnsBase->flagTurb==1)
			{
			printf("\nDo you want to include Turbulent Dispersion Term in Volumefraction equation If YES Enter-> 1  No Enter-> 0 : ");
			scanf("%d",&flagDispersion);            
                        	
                        var2=fopen("Settings/coefficients.dat","a");
                        var3=fopen("Settings/ProgramReadInputValues.dat","a");
                        fprintf(var2,"%d\n",flagDispersion);
                        fprintf(var3,"%d \tDo you want to include Turbulent Dispersion Term in Volumefraction equation If YES Enter-> 1  No Enter-> 0 \n\n",flagDispersion);
                        fclose(var2);
                        fclose(var3);
        
			CgnsRoot.CgnsBase->flagDispersion=flagDispersion;
			
				if(CgnsRoot.CgnsBase->flagDispersion==1)
				{
				printf("\nEnter the value of turbulent prandtl number : ");
				scanf("%lf",&number);            
                                	
                                var2=fopen("Settings/coefficients.dat","a");
                                var3=fopen("Settings/ProgramReadInputValues.dat","a");
                                fprintf(var2,"%lf\n",number);
                                fprintf(var3,"%lf \tEnter the value of turbulent prandtl number \n\n",number);
                                fclose(var2);
                                fclose(var3);
        
				CgnsRoot.CgnsBase->number=number;
				
				}
				
				if(CgnsRoot.CgnsBase->BaseParameters->Dragflag==1)
				{
					do
					{
						printf("\nDo you want to include Liftforce also as Interaction term(0 ->no / 1->yes): ");
						scanf("%d",&(CgnsRoot.CgnsBase->flagLift));
					}while(!(CgnsRoot.CgnsBase->flagLift == 1 || CgnsRoot.CgnsBase->flagLift == 0));            
                                        	
                                        var2=fopen("Settings/coefficients.dat","a");
                                        var3=fopen("Settings/ProgramReadInputValues.dat","a");
                                        fprintf(var2,"%d\n",CgnsRoot.CgnsBase->flagLift);
                                        fprintf(var3,"%d \tDo you want to include Liftforce also as Interaction term(0 ->no / 1->yes):  \n\n",CgnsRoot.CgnsBase->flagLift);
                                        fclose(var2);
                                        fclose(var3);
				}
			}
				//End--Chandrasekhar		
		}	
		else if(Multiflowtype==2)
		{
			CgnsRoot.CgnsBase->flagParticle = 1 ;  
			CgnsRoot.CgnsBase->BaseParameters->flagVariableDia = 0 ;
			CgnsRoot.CgnsBase->BaseParameters->flagDropFuelMassFrac =0;
			printf("\nFor the  drag calculation, If value of Drag coefficient(CD) is constant Enter 0\n\nOtherwise for user defined function for Drag coefficient(CD) Enter 1 : ");
			scanf("%d",&constant);                                
                        	
                        var2=fopen("Settings/coefficients.dat","a");
                        var3=fopen("Settings/ProgramReadInputValues.dat","a");
                        fprintf(var2,"%d\n",constant);
                        fprintf(var3,"%d \tConstant Drag coefficient(CD) Enter 0, User defined function Enter 1\n\n",constant);
                        fclose(var2);
                        fclose(var3);
        
			if(constant==0)
			{
				printf("\nEnter the constant  value of drag coefficinet(CD) : ");
                                scanf("%lf",&C);                                
                                	
                                var2=fopen("Settings/coefficients.dat","a");
                                var3=fopen("Settings/ProgramReadInputValues.dat","a");
                                fprintf(var2,"%lf\n",C);
                                fprintf(var3,"%lf \tEnter the constant  value of drag coefficinet(CD) \n\n",C);
                                fclose(var2);
                                fclose(var3);
        
				CgnsRoot.CgnsBase->BaseParameters->ConstantCD = C;	
			}
			else
			{
			printf("\n It will directly go to user defined function for CD calculation\n");
			printf("\nUser have freedom to change  user defined function in (user/APKuserDefinedVarCoefficients/userDefinedCD)\n\n");
			CgnsRoot.CgnsBase->constant = constant;
			}
		}
		else if(Multiflowtype==3)
		{
			
			CgnsRoot.CgnsBase->flagDroplet = 1 ;  
                   	printf("\nIf you want to consider evaporation of droplets enter 1 \n else enter 0: ");
		   	scanf("%d",&flagevaporation);                                
                        	
                        var2=fopen("Settings/coefficients.dat","a");
                        var3=fopen("Settings/ProgramReadInputValues.dat","a");
                        fprintf(var2,"%d\n",flagevaporation);
                        fprintf(var3,"%d \tIf you want to consider evaporation of droplets enter 1 \n else enter 0: \n\n",flagevaporation);
                        fclose(var2);
                        fclose(var3);
        
		   	CgnsRoot.CgnsBase->BaseParameters->flagVariableDia = flagevaporation ;
		   
		   	printf("\n If you want to solve evaporated fuel mass fraction equation enter 1 \n else enter 0: \n(If you have considered evaporation of the droplets then solve fuel mass fraction)\n");
		   	scanf("%d",&flagdropfuelmassfrac);                                
                        	
                        var2=fopen("Settings/coefficients.dat","a");
                        var3=fopen("Settings/ProgramReadInputValues.dat","a");
                        fprintf(var2,"%d\n",flagdropfuelmassfrac);
                        fprintf(var3,"%d \tIf you want to solve evaporated fuel mass fraction equation enter 1 else enter 0: \n\n",flagdropfuelmassfrac);
                        fclose(var2);
                        fclose(var3);
        
		   	CgnsRoot.CgnsBase->BaseParameters->flagDropFuelMassFrac = flagdropfuelmassfrac;
			
			printf("\nFor the  drag calculation, If value of Drag coefficient(CD) is constant Enter 0\n\nOtherwise for user defined function for Drag coefficient(CD) Enter 1\n");
			scanf("%d",&constant);                                
                        	
                        var2=fopen("Settings/coefficients.dat","a");
                        var3=fopen("Settings/ProgramReadInputValues.dat","a");
                        fprintf(var2,"%d\n",constant);
                        fprintf(var3,"%d \tConstant Drag coefficient(CD) Enter 0, User defined function Enter 1 \n\n",constant);
                        fclose(var2);
                        fclose(var3);
        
			if(constant==0)
			{
				printf("\n Enter the constant  value of drag coefficinet(CD)");			
                                scanf("%lf",&C);                                
                                                               
                                var2=fopen("Settings/coefficients.dat","a");
                                var3=fopen("Settings/ProgramReadInputValues.dat","a");
                                fprintf(var2,"%lf\n",C);
                                fprintf(var3,"%lf \tEnter the constant  value of drag coefficinet(CD) \n\n",C);
                                fclose(var2);
                                fclose(var3);
                                
				CgnsRoot.CgnsBase->BaseParameters->ConstantCD = C;	
			}
			else
			{
			printf("\n It will directly go to user defined function for CD calculation\n");
			printf("\nUser have freedom to change  user defined function in (user/APKuserDefinedVarCoefficients/userDefinedCD)\n\n");
			CgnsRoot.CgnsBase->constant = constant;
			}
			
		}
		else
		{
			CgnsRoot.CgnsBase->BaseParameters->flagVariableDia = 0;
		  	CgnsRoot.CgnsBase->BaseParameters->flagDropFuelMassFrac =0; 
		}
	}

	/**********************************************************/	
	

	if(CgnsRoot.CgnsBase->flagTurb==0&&CgnsRoot.CgnsBase->flagGroundwater==0)
	{
		//added by bhans for solidification problem START
		do
		{
			//printf("\nSpecify whether the problem Solidification (or) melting problem(0 ->no / 1->yes): ");
			//scanf ("%d",&flagsoli);
			flagsoli = 0;
		}while(!(flagsoli == 1 || flagsoli == 0));                                
               	
                //var2=fopen("Settings/coefficients.dat","a");
                //var3=fopen("Settings/ProgramReadInputValues.dat","a");
                //fprintf(var2,"%d\n",flagsoli);
                //fprintf(var3,"%d \tSpecify whether the problem Solidification (or) melting problem(0 ->no / 1->yes): \n\n",flagsoli);
                //fclose(var2);
                //fclose(var3);
        
		CgnsRoot.CgnsBase->flagSolidification=flagsoli;
		//added by bhans for solidification problem END
		
		//added by anji for twophase (flow boiling problem) START
		do
		{
			//printf("\nSpecify whether the problem is related to Phase change (flow boiling)Two-phase problem(0 ->no / 1->yes): ");
			//scanf ("%d",&flagTP);
			flagTP = 0;
		}while(!(flagTP == 1 || flagTP == 0));                                
                	
                //var2=fopen("Settings/coefficients.dat","a");
                //var3=fopen("Settings/ProgramReadInputValues.dat","a");
                //fprintf(var2,"%d\n",flagTP);
                //fprintf(var3,"%d \tSpecify whether the problem is related to Phase change (flow boiling)Two-phase problem(0 ->no / 1->yes): \n\n",flagTP);
                //fclose(var2);
                //fclose(var3);
        
		CgnsRoot.CgnsBase->flagTwoPhase=flagTP;
		//added by anji for twophase (flow boiling problem) END
		
		//added by naveen for LMMHD problem START
		do
		{
			printf("\nSpecify whether you want to solve a liquid metal MHD flow problem (0 ->no / 1->yes): ");
	                scanf ("%d",&flagLMMHD);
		}while(!(flagLMMHD == 1 || flagLMMHD == 0));
			
                var2=fopen("Settings/coefficients.dat","a");
                var3=fopen("Settings/ProgramReadInputValues.dat","a");
                fprintf(var2,"%d\n",flagLMMHD);
                fprintf(var3,"%d \tSpecify whether you want to solve a liquid metal MHD flow problem (0 ->no / 1->yes): \n\n",flagLMMHD);
                fclose(var2);
                fclose(var3);
        
		CgnsRoot.CgnsBase->flagLMMHD=flagLMMHD;
		//Bydefault the tag is made zero.
		CgnsRoot.CgnsBase->flaginduced = 0;
		CgnsRoot.CgnsBase->flagWf = 0;
		//added by naveen for LMMHD problem END
		if(flagLMMHD==1)
			{
				do
				{
					printf("\nDo you want to solve for the induced magnetic field also (0 ->no / 1->yes):");
					scanf("%d",&flaginduced);
				}while(!(flaginduced == 1 || flaginduced == 0));                                
                                	
                                var2=fopen("Settings/coefficients.dat","a");
                                var3=fopen("Settings/ProgramReadInputValues.dat","a");
                                fprintf(var2,"%d\n",flaginduced);
                                fprintf(var3,"%d \tDo you want to solve for the induced magnetic field also (0 ->no / 1->yes): \n\n",flaginduced);
                                fclose(var2);
                                fclose(var3);
        
				CgnsRoot.CgnsBase->flaginduced=flaginduced;
			if(CgnsRoot.CgnsBase->BaseParameters->flag2d==0 && flaginduced==0)
			{
				do
				{	printf("\nDo you want to apply Wall-functions (0 ->no / 1->yes):");
					printf("\nNote that the boundary needs to be tagged as Wall Viscous to apply wall-functions:");
					scanf("%d",&flagWf);
				}while(!(flagWf == 1 || flagWf == 0));                                
                                	
                                var2=fopen("Settings/coefficients.dat","a");
                                var3=fopen("Settings/ProgramReadInputValues.dat","a");
                                fprintf(var2,"%d\n",flagWf);
                                fprintf(var3,"%d \tDo you want to apply Wall-functions (0 ->no / 1->yes): \n\n",flagWf);
                                fclose(var2);
                                fclose(var3);
        
				CgnsRoot.CgnsBase->flagWf=flagWf;

			}
			}

               	//added by rajendra for LMMHD problem END
		
		if(CgnsRoot.CgnsBase->flagMPFV == 1)
		{
			for(iZone=0;iZone<numzones;iZone++)
			{
				zone = CgnsRoot.CgnsBase->CgnsZone+iZone;
				
				zone->FlagMPFV = flagMultiPhase;
			}
			CgnsRoot.CgnsBase->flagVariableDensityProblem = 0;
		}
	
		if((CgnsRoot.CgnsBase->flagParticle==1)||(CgnsRoot.CgnsBase->flagDroplet==1))
		{
			//printf("\nPlease choose whether you are solving for a\nVariable Density case (Enter 1)\nOr a Constant Fluid Properties case (Enter 0) for the carrier phase: \n");
			
			//scanf("%d",&varRho);                                
			varRho = 0;                                
                        	
                        //var2=fopen("Settings/coefficients.dat","a");
                        //var3=fopen("Settings/ProgramReadInputValues.dat","a");
                        //fprintf(var2,"%d\n",varRho);
                        //fprintf(var3,"%d \tVariable Density case (Enter 1) or a Constant Fluid Properties case (Enter 0) for the carrier phase: \n\n",varRho);
                        //fclose(var2);
                        //fclose(var3);
        
			CgnsRoot.CgnsBase->flagVariableDensityProblem = varRho;
			
			printf("\n If this is a Constant Density problem with disperse flow then please choose SEMICOUPLED solver only: \n");
		}
		else
		{
			//printf("\nPlease choose whether you are solving for a\nVariable Density case (Enter 1)\nOr a Constant Fluid Properties case (Enter 0):");
			//scanf("%d",&varRho);                                
			varRho = 0;                                
                        	
                        //var2=fopen("Settings/coefficients.dat","a");
                        //var3=fopen("Settings/ProgramReadInputValues.dat","a");
                        //fprintf(var2,"%d\n",varRho);
                        //fprintf(var3,"%d \tVariable Density case (Enter 1) or a Constant Fluid Properties case (Enter 0): \n\n",varRho);
                        //fclose(var2);
                        //fclose(var3);
        
			CgnsRoot.CgnsBase->flagVariableDensityProblem = varRho;
		}
		
		if(varRho==0)
		{
			printf("\nPlease Enter the SOLUTION SCHEME:-\nSEMICOUPLED->[Enter 1]:\nSEGREGATED-->[Enter 2]:");
			scanf("%d",&scheme);                                
                        	
                        var2=fopen("Settings/coefficients.dat","a");
                        var3=fopen("Settings/ProgramReadInputValues.dat","a");
                        fprintf(var2,"%d\n",scheme);
                        fprintf(var3,"%d \tPlease Enter the SOLUTION SCHEME:-SEMICOUPLED->[Enter 1]: SEGREGATED->[Enter 2]: \n\n",scheme);
                        fclose(var2);
                        fclose(var3);
        
			CgnsRoot.CgnsBase->couplingTag=scheme;
		}
		else
		{
			printf("\nPlease Enter the SOLUTION SCHEME:-\nSEMICOUPLED->[Enter 1]:\nSEGREGATED-->[Enter 2]:");
			printf("\nSince this is a Variable Density problem, please choose SEMICOULED solver only :");
			scanf("%d",&scheme);                        
                        	
                        var2=fopen("Settings/coefficients.dat","a");
                        var3=fopen("Settings/ProgramReadInputValues.dat","a");
                        fprintf(var2,"%d\n",scheme);
                        fprintf(var3,"%d \tPlease Enter the SOLUTION SCHEME:-SEMICOUPLED->[Enter 1]: SEGREGATED->[Enter 2]: \n\n",scheme);
                        fclose(var2);
                        fclose(var3);
        
			CgnsRoot.CgnsBase->couplingTag=scheme;
			if(CgnsRoot.CgnsBase->flagMPFV==0)//Added by sushil
			{
				printf("\nIf Density variation is due to fluid mixing Enter 1\n");
				printf("If Density variation is due to Temperature Enter 0:");
				scanf("%d",&densityvar);                        
                                	
                                var2=fopen("Settings/coefficients.dat","a");
                                var3=fopen("Settings/ProgramReadInputValues.dat","a");
                                fprintf(var2,"%d\n",densityvar);
                                fprintf(var3,"%d \tIf Density variation is due to fluid mixing Enter 1, due to Temperature Enter 0:\n\n",densityvar);
                                fclose(var2);
                                fclose(var3);
        
				CgnsRoot.CgnsBase->flagtypeVariableDensity=densityvar;
				scheme = 1; //Currently we are solving Variable density problems using semicoupled solver only. We are explicitly making the solver SemiCoupled.
				CgnsRoot.CgnsBase->couplingTag=scheme;
			}
		}
		
		if(scheme==0 || scheme==1)
		{
			printf("\nEnter Solver Scheme if:");
			printf("\nIf Implicit-----[Enter 0]");
			printf("\nIf SemiImplicit-[Enter 1] : ");
			scanf("%d",&scheme);                        
                        	
                        var2=fopen("Settings/coefficients.dat","a");
                        var3=fopen("Settings/ProgramReadInputValues.dat","a");
                        fprintf(var2,"%d\n",scheme);
                        fprintf(var3,"%d \tEnter Solver Scheme if Implicit-0: SemiImplicit-1 \n\n",scheme);
                        fclose(var2);
                        fclose(var3);
        
			CgnsRoot.CgnsBase->SolverScheme = scheme;
		}
		else
		{
			printf("\nEnter Solver Scheme if:");
			printf("\nIf Implicit--------------------------[Enter  0]");
			printf("\nIf SemiImplicit----------------------[Enter  1]");
			printf("\nIf PseudoCompressiblity");
			printf("\nIf (a)one T and one Tau step---------[Enter 21]");
			printf("\nIf (b)Impose Continuity at each step-[Enter 22] : ");
			scanf("%d",&scheme);                        
                        	
                        var2=fopen("Settings/coefficients.dat","a");
                        var3=fopen("Settings/ProgramReadInputValues.dat","a");
                        fprintf(var2,"%d\n",scheme);
                        fprintf(var3,"%d \tEnter Solver Scheme Implicit-0 SemiImplicit-1 PseudoCompressiblity (a)one T and one Tau step-21 (b)Impose Continuity at each step-22 \n\n",scheme);
                        fclose(var2);
                        fclose(var3);
        
			CgnsRoot.CgnsBase->SolverScheme = scheme;
		}
	}
	else if(CgnsRoot.CgnsBase->flagTurb==0&&CgnsRoot.CgnsBase->flagGroundwater==1)
	{
		printf("\nPlease Enter the SOLUTION SCHEME:-\nSEMICOUPLED->[Enter 1]:\nSEGREGATED-->[Enter 2]:");
		printf("\nFor Groundwater, please choose SEMICOUPLED solver only :");
		scanf("%d",&scheme);                        
                	
                var2=fopen("Settings/coefficients.dat","a");
                var3=fopen("Settings/ProgramReadInputValues.dat","a");
                fprintf(var2,"%d\n",scheme);
                fprintf(var3,"%d \tPlease Enter the SOLUTION SCHEME:-SEMICOUPLED->[Enter 1]: SEGREGATED->[Enter 2]: \n\n",scheme);
                fclose(var2);
                fclose(var3);
        
		CgnsRoot.CgnsBase->couplingTag=scheme;

		printf("\nEnter Solver Scheme if:");
		printf("\nIf Implicit-----[Enter 0]");
		printf("\nIf SemiImplicit-[Enter 1] : ");
		scanf("%d",&scheme);                        
                	
                var2=fopen("Settings/coefficients.dat","a");
                var3=fopen("Settings/ProgramReadInputValues.dat","a");
                fprintf(var2,"%d\n",scheme);
                fprintf(var3,"%d \tEnter Solver Scheme Implicit-0 SemiImplicit-1 \n\n",scheme);
                fclose(var2);
                fclose(var3);
        
		CgnsRoot.CgnsBase->SolverScheme = scheme;
		
		
		printf("\n For Groundwater problem Fully implicit scheme with Semi Coupled Solver\n is applied by default\n");
		CgnsRoot.CgnsBase->SolverScheme=0; /*Implicit Solver*/
		CgnsRoot.CgnsBase->couplingTag=1; /*Semi Coupled Solver*/
		printf(" ******************************************************************************\n");
		
		printf("\nPlease choose whether you are solving for \n a Rigid Porous medium or Incompressible fluid(Enter 0)\n or for Deformable porous medium or Compressible flow a case (Enter 1):");
		scanf("%d",&flagGr);                        
                	
                var2=fopen("Settings/coefficients.dat","a");
                var3=fopen("Settings/ProgramReadInputValues.dat","a");
                fprintf(var2,"%d\n",flagGr);
                fprintf(var3,"%d \tRigid Porous medium or Incompressible fluid(Enter 0) or for Deformable porous medium or Compressible flow a case (Enter 1): \n\n",flagGr);
                fclose(var2);
                fclose(var3);
        
		CgnsRoot.CgnsBase->flagdeformable = flagGr;
		
		if(CgnsRoot.CgnsBase->BaseParameters->flag2d==0)
		{
			 CgnsRoot.CgnsBase->BaseParameters->flagdat =1;
		}
		else
		{
			printf("\nenter `0'if the plane is along the datum plane or `1'if in other plane ");
			scanf("%d",&flagGr);                        
                        	
                        var2=fopen("Settings/coefficients.dat","a");
                        var3=fopen("Settings/ProgramReadInputValues.dat","a");
                        fprintf(var2,"%d\n",flagGr);
                        fprintf(var3,"%d \tEnter `0'if the plane is along the datum plane or `1'if in other plane \n\n",flagGr);
                        fclose(var2);
                        fclose(var3);
        
			CgnsRoot.CgnsBase->BaseParameters->flagdat=flagGr;
		}
		
		printf("\nEnter `0' if the problem is of Saturated flow");
		printf("\nEnter `1' if the problem is of UnSaturated flow or Variably saturated flow");
		scanf ("%d",&flagGr);                        
                	
                var2=fopen("Settings/coefficients.dat","a");
                var3=fopen("Settings/ProgramReadInputValues.dat","a");
                fprintf(var2,"%d\n",flagGr);
                fprintf(var3,"%d \tSaturated flow Enter `0': UnSaturated flow or Variably saturated flow Enter `1'\n\n",flagGr);
                fclose(var2);
                fclose(var3);
        
		CgnsRoot.CgnsBase->flagSaturation=flagGr;
	
	}
	else if(CgnsRoot.CgnsBase->flagTurb==1)//Turbulence area, START
	{
		
		//added by bhans for solidification problem START
		do
		{
			//printf("\nSpecify whether the problem Solidification (or) melting problem(0 ->no / 1->yes): ");
			//scanf ("%d",&flagsoli);
			flagsoli = 0;
		}while(!(flagsoli == 1 || flagsoli == 0));                        
                	
                //var2=fopen("Settings/coefficients.dat","a");
                //var3=fopen("Settings/ProgramReadInputValues.dat","a");
                //fprintf(var2,"%d\n",flagsoli);
                //fprintf(var3,"%d \tSpecify whether the problem Solidification (or) melting problem(0 ->no / 1->yes): \n\n",flagsoli);
                //fclose(var2);
                //fclose(var3);
        
		CgnsRoot.CgnsBase->flagSolidification=flagsoli;
		//added by bhans for solidification problem END
		
		
		//added by anji for twophase (flow boiling problem) START
		do
		{
			//printf("\nSpecify whether the problem is related to Phase change (flow boiling)Two-phase problem(0 ->no / 1->yes): ");
			//scanf ("%d",&flagTP);
			flagTP = 0;
		}while(!(flagTP == 1 || flagTP == 0));                        
                	
                //var2=fopen("Settings/coefficients.dat","a");
                //var3=fopen("Settings/ProgramReadInputValues.dat","a");
                //fprintf(var2,"%d\n",flagTP);
                //fprintf(var3,"%d \tSpecify whether the problem is related to Phase change (flow boiling)Two-phase problem(0 ->no / 1->yes): \n\n",flagTP);
                //fclose(var2);
                //fclose(var3);
        
		CgnsRoot.CgnsBase->flagTwoPhase=flagTP;
		//added by anji for twophase (flow boiling problem) END
		
		//added by chandrasekhar for Turbulent MPFV START
		if(CgnsRoot.CgnsBase->flagMPFV == 1)
		{
			for(iZone=0;iZone<numzones;iZone++)
			{
				zone = CgnsRoot.CgnsBase->CgnsZone+iZone;
				
				zone->FlagMPFV = flagMultiPhase;
			}
			CgnsRoot.CgnsBase->flagVariableDensityProblem = 0;
			
		}
		//added by chandrasekhar for Turbulent MPFV END

                do//Shyam
		{
			//printf("\nSpecify whether the problem is Turbulent Combustion(0 ->no / 1->yes): ");
			//scanf("%d",&(CgnsRoot.CgnsBase->flagTurbComb));
			CgnsRoot.CgnsBase->flagTurbComb = 0;
                }while(!(CgnsRoot.CgnsBase->flagTurbComb == 1 || CgnsRoot.CgnsBase->flagTurbComb == 0));                        
                	
                //var2=fopen("Settings/coefficients.dat","a");
                //var3=fopen("Settings/ProgramReadInputValues.dat","a");
                //fprintf(var2,"%d\n",CgnsRoot.CgnsBase->flagTurbComb);
                //fprintf(var3,"%d \tSpecify whether the problem is Turbulent Combustion(0 ->no / 1->yes): \n\n",CgnsRoot.CgnsBase->flagTurbComb);
                //fclose(var2);
                //fclose(var3);
		//Shyam

                if(CgnsRoot.CgnsBase->flagTurbComb==1)//Shyam, START
                {
// 			do
// 			{
// 			printf("Enter the following options for Turbulent Reactive Flows:\n");
// 			printf("1 - for Premixed Turbulent Reactive Flow\n");
// 			printf("2 - for Non-Premixed(Diffusion) Turbulent Reactive Flow\n\n");
// 			scanf("%d",&(CgnsRoot.CgnsBase->flagPreNonPre));
// 			}while(!(CgnsRoot.CgnsBase->flagPreNonPre==1||CgnsRoot.CgnsBase->flagPreNonPre==2));	
//                    var2=fopen("Settings/coefficients.dat","a");
//                    var3=fopen("Settings/ProgramReadInputValues.dat","a");
//                    fprintf(var2,"%d\n",CgnsRoot.CgnsBase->flagPreNonPre);
//                    fprintf(var3,"%d \t \n\n",CgnsRoot.CgnsBase->flagPreNonPre);
//                    fclose(var2);
//                    fclose(var3);
// 			
// 			if(CgnsRoot.CgnsBase->flagPreNonPre==1)
// 			{
// 				CgnsRoot.CgnsBase->flagDiffKinetCont=2;
// 				printf("Kinetically Controlled Reactions Method is used for Premixed Reactive Flows\n");
// 			}
			CgnsRoot.CgnsBase->flagPreNonPre = 2;
			if(CgnsRoot.CgnsBase->flagPreNonPre==2)
			{
				do
				{
					printf("Enter the following options for Non-Premixed Turbulent Reactive Flows:\n");
					printf("1 - for Diffusion Controlled\n");
					printf("2 - for Kinetically Controlled\n\n");
					scanf("%d",&(CgnsRoot.CgnsBase->flagDiffKinetCont));
				}while(!(CgnsRoot.CgnsBase->flagDiffKinetCont==1||CgnsRoot.CgnsBase->flagDiffKinetCont==2));                
                                	
                                var2=fopen("Settings/coefficients.dat","a");
                                var3=fopen("Settings/ProgramReadInputValues.dat","a");
                                fprintf(var2,"%d\n",CgnsRoot.CgnsBase->flagDiffKinetCont);
                                fprintf(var3,"%d \t1 - for Diffusion Controlled: 2 - for Kinetically Controlled Non-Premixed Turbulent Reactive Flows \n\n",CgnsRoot.CgnsBase->flagDiffKinetCont);
                                fclose(var2);
                                fclose(var3);
			}
		}//Shyam, END
		
		do
		{
			//printf("\nPlease choose whether you are solving for a\nVariable Density case (Enter 1)\nOr a Constant Fluid Properties case (Enter 0):");
			//scanf("%d",&varRho);
			varRho = 0;
		}while(!(varRho == 1 || varRho == 0));                                
                	
                //var2=fopen("Settings/coefficients.dat","a");
                //var3=fopen("Settings/ProgramReadInputValues.dat","a");
                //fprintf(var2,"%d\n",varRho);
                //fprintf(var3,"%d \tVariable Density case (Enter 1) or a Constant Fluid Properties case (Enter 0): \n\n",varRho);
                //fclose(var2);
                //fclose(var3);
        
		CgnsRoot.CgnsBase->flagVariableDensityProblem = varRho;
		//printf("\nCgnsRoot.CgnsBase->flagVariableDensityProblem=%d",CgnsRoot.CgnsBase->flagVariableDensityProblem);

		if (CgnsRoot.CgnsBase->flagTurb==0) {
			printf("\nPlease Enter the SOLUTION SCHEME:-\nSEMICOUPLED->[Enter 1]:\nSEGREGATED-->[Enter 2]:");
			//printf("\nFor Turbulence problems, please choose SEMICOULED solver only :");
			//printf("\nFor Variable Density and/or Turbulence problems, please choose SEMICOULED solver only :");
			scanf("%d",&scheme);                
                	
                	var2=fopen("Settings/coefficients.dat","a");
                	var3=fopen("Settings/ProgramReadInputValues.dat","a");
                	fprintf(var2,"%d\n",scheme);
                	fprintf(var3,"%d \tPlease Enter the SOLUTION SCHEME:-SEMICOUPLED->[Enter 1]: SEGREGATED->[Enter 2]: \n\n",scheme);
                	fclose(var2);
                	fclose(var3);
		} else { scheme = 1; }
        
		CgnsRoot.CgnsBase->couplingTag=scheme;	
		printf("\nEnter Solver Scheme if:");
		printf("\nIf Implicit-----[Enter 0]");
		printf("\nIf SemiImplicit-[Enter 1] : ");
		scanf("%d",&scheme);                
                	
                var2=fopen("Settings/coefficients.dat","a");
                var3=fopen("Settings/ProgramReadInputValues.dat","a");
                fprintf(var2,"%d\n",scheme);
                fprintf(var3,"%d \tEnter Solver Scheme Implicit-0: SemiImplicit-1 \n\n",scheme);
                fclose(var2);
                fclose(var3);
        
		CgnsRoot.CgnsBase->SolverScheme = scheme;
		
		if(varRho==1)
		{
			do
			{
				printf("\nIf Density variation is due to fluid mixing Enter 1\n");
				printf("If Density variation is due to Temperature Enter 0\n");
				printf("NOTE: Choose only '0' for Turbulent Combustion problem:");
				scanf("%d",&densityvar);
			}while(!(densityvar == 1 || densityvar == 0));                
                        	
                        var2=fopen("Settings/coefficients.dat","a");
                        var3=fopen("Settings/ProgramReadInputValues.dat","a");
                        fprintf(var2,"%d\n",densityvar);
                        fprintf(var3,"%d \tIf Density variation is due to fluid mixing Enter 1: due to Temperature Enter 0 \n\n",densityvar);
                        fclose(var2);
                        fclose(var3);
        
			CgnsRoot.CgnsBase->flagtypeVariableDensity=densityvar;
		}
		
		//NOTE: Use Semicoupled solver for variable density problem.
		printf("\n ******************************************************************************");
		printf("\n Since Turbulence problem is choosen, fully implicit scheme with Semi Coupled Solver\n are applied by default even if the user enters other solvers\n");
		CgnsRoot.CgnsBase->SolverScheme=0; //Implicit Solver
		CgnsRoot.CgnsBase->couplingTag=1; //Semi Coupled Solver
		printf(" ******************************************************************************\n");
		
		do
		{
                	//printf("NOTE: Please choose only standard K-Epsilon model for Turbulent Combustion problem \n");
			//printf("NOTE: Please choose Standard K-Epsilon model for MPFV(Liquid-Gas) or Gas-Particle flow problem \n");
			printf("Choose the Turbulence Model:\n");
			
			//printf("\n############### ZERO-EQUATION TURBULENCE MODEL(S) ###############\n");
			//printf("Baldwin-Lomax Model (1978) [Enter-9]\n");
			
			printf("\n############### TWO-EQUATION TURBULENCE MODEL(S) ###############\n");
			printf("Standard K-Epsilon [Enter-1] \n");
			printf("Standard K-Omega [Enter-2]\n");
			printf("RNG K-Epsilon [Enter-3]\n");
			printf("Low Reynolds Number K-Epsilon: Launder and Sharma [Enter-4]\n");
			printf("Low Reynolds Number K-Epsilon: Lam - Bremhorst [Enter-5]\n");
			printf("Shear Stress Transport (SST) Model [Enter-6]\n");
			
			printf("\n############### ONE-EQUATION TURBULENCE MODEL(S) ###############\n");
			printf("Spallart-Almaras (ONE-equation) Model [Enter-8]\n");
			
			//printf("\n############### FOUR-EQUATION TURBULENCE MODEL(S) ###############\n");
			//printf("k-epsilon-v2-f (FOUR-equation) Model [Enter-7]\n");
			
			scanf("%d",&flagT);
			printf("You have entered %d\n", flagT);
		}while(!((flagT >= 1 && flagT <= 6) || flagT == 8));// || flagT == 9));
		//while(!(flagT == 1 || flagT == 2 || flagT == 3 || flagT == 4 || flagT == 5 || flagT == 6 || flagT == 8 || flagT == 9));                
                	
                var2=fopen("Settings/coefficients.dat","a");
                var3=fopen("Settings/ProgramReadInputValues.dat","a");
                fprintf(var2,"%d\n",flagT);
                fprintf(var3,"%d \tChoice of Turbulence Model \n\n",flagT);
                fclose(var2);
                fclose(var3);
        
		CgnsRoot.CgnsBase->flagTurbModel=flagT;
		//printf("\n CgnsRoot.CgnsBase->flagTurbModel=%d",CgnsRoot.CgnsBase->flagTurbModel);
		if(CgnsRoot.CgnsBase->flagTurbModel==7)		//v2f turbulence model
			CgnsRoot.CgnsBase->BaseParameters->Turb->numTurbScalars = 4;
		else if(CgnsRoot.CgnsBase->flagTurbModel==8)	//Spalart-Allmaras turbulence model
			CgnsRoot.CgnsBase->BaseParameters->Turb->numTurbScalars = 1;
		else if(CgnsRoot.CgnsBase->flagTurbModel==9)	//Baldwin-Lomax zero equation turbulence model
			CgnsRoot.CgnsBase->BaseParameters->Turb->numTurbScalars = 0;
		else
			CgnsRoot.CgnsBase->BaseParameters->Turb->numTurbScalars = 2;	//Other turbulence models
		
		if(CgnsRoot.CgnsBase->flagTurbModel == 1)	//If turbulence model choosen is K-Epsilon Model
		{
			do
			{
                        	printf("Do you want to solve K-Epsilon equations in their LOGARITHMIC form [1->YES/0->NO]:\n");
				//printf("NOTE: Please don't choose LOGARITHMIC form for Turbulent Combustion problem or MPFV(Liquid-Gas) problem or Gas-Particle problem\n");
				scanf("%d", &flagLogKEps);
			}while(!(flagLogKEps == 1 || flagLogKEps == 0));                
                        	
                        var2=fopen("Settings/coefficients.dat","a");
                        var3=fopen("Settings/ProgramReadInputValues.dat","a");
                        fprintf(var2,"%d\n",flagLogKEps);
                        fprintf(var3,"%d \tDo you want to solve K-Epsilon equations in their LOGARITHMIC form [1->YES/0->NO]: \n\n",flagLogKEps);
                        fclose(var2);
                        fclose(var3);
        
			CgnsRoot.CgnsBase->flagLogKEps = flagLogKEps;
			printf("********** CgnsRoot.CgnsBase->flagLogKEps = %d ***********\n", CgnsRoot.CgnsBase->flagLogKEps);
		
// 			if(CgnsRoot.CgnsBase->flagMPFV == 1)
// 			{
// 				printf("Enter Type of model you want to use for Liquid-gas flow:\nEnetr 0 for Mixture turbulence model and Enetr 1 for Dispersed turbulence model\n");
// 				scanf("%d", &flagMultiPhase);
// 				for(iZone=0;iZone<numzones;iZone++)
// 				{
// 					zone = CgnsRoot.CgnsBase->CgnsZone+iZone;
// 				
// 					zone->flagTurbMPFVmodel = flagMultiPhase;
// 				}
// 				CgnsRoot.CgnsBase->flagTurbMPFVmodel =flagMultiPhase ;
// 			}
			if(CgnsRoot.CgnsBase->flagParticle  == 1)
			{
				printf("\n Do you want to include Turbulence Modulation term:[1->YES/0->NO]\n ");
				scanf("%d", &flagMultiPhase);                
                                	
                                var2=fopen("Settings/coefficients.dat","a");
                                var3=fopen("Settings/ProgramReadInputValues.dat","a");
                                fprintf(var2,"%d\n",flagMultiPhase);
                                fprintf(var3,"%d \tDo you want to include Turbulence Modulation term:[1->YES/0->NO] \n\n",flagMultiPhase);
                                fclose(var2);
                                fclose(var3);
        
				CgnsRoot.CgnsBase->flagModulation =flagMultiPhase ;
			}
		}
		if(CgnsRoot.CgnsBase->flagTurbModel == 3)
		{
			if(CgnsRoot.CgnsBase->flagParticle  == 1)
			{
				printf("\n Do you want to include Turbulence Modulation term:[1->YES/0->NO]\n ");
				scanf("%d", &flagMultiPhase);                
                                	
                                var2=fopen("Settings/coefficients.dat","a");
                                var3=fopen("Settings/ProgramReadInputValues.dat","a");
                                fprintf(var2,"%d\n",flagMultiPhase);
                                fprintf(var3,"%d \tDo you want to include Turbulence Modulation term:[1->YES/0->NO] \n\n",flagMultiPhase);
                                fclose(var2);
                                fclose(var3);
        
				CgnsRoot.CgnsBase->flagModulation =flagMultiPhase ;
			}
		}	
		if(CgnsRoot.CgnsBase->BaseParameters->Turb->numTurbScalars >= 2)
		{
			printf("Enter the number of dual time step iterations to be performed for solving the governing equations for turbulence variables\n");
			do
			{
				scanf("%d", &CgnsRoot.CgnsBase->BaseParameters->Turb->numDualTimeIter);
			}while(CgnsRoot.CgnsBase->BaseParameters->Turb->numDualTimeIter < 1);                
                        	
                        var2=fopen("Settings/coefficients.dat","a");
                        var3=fopen("Settings/ProgramReadInputValues.dat","a");
                        fprintf(var2,"%d\n",CgnsRoot.CgnsBase->BaseParameters->Turb->numDualTimeIter);
                        fprintf(var3,"%d \tEnter the number of dual time step iterations to be performed for solving the governing equations for turbulence variables \n\n",CgnsRoot.CgnsBase->BaseParameters->Turb->numDualTimeIter);
                        fclose(var2);
                        fclose(var3);
			//CgnsRoot.CgnsBase->BaseParameters->Turb->numDualTimeIter = 10;
		}
		else
		{
			CgnsRoot.CgnsBase->BaseParameters->Turb->numDualTimeIter = 1;
		}
		if(CgnsRoot.CgnsBase->flagTurbModel == 1)
		{
			printf("Choose 1->YES/0->NO for the following limiters used standard k-epsilon model\n\n");
			printf("Do you want to apply limiters on turbulence time-scale (Durbin's Realizability constraints) (k/epsilon) [0->NO / 1-YES]: \n");
			do
			{
				scanf("%d", &flagT);
			}while(!(flagT==0 || flagT==1));                                
                        	
                        var2=fopen("Settings/coefficients.dat","a");
                        var3=fopen("Settings/ProgramReadInputValues.dat","a");
                        fprintf(var2,"%d\n",flagT);
                        fprintf(var3,"%d \tDo you want to apply limiters on turbulence time-scale (Durbin's Realizability constraints) (k/epsilon) [0->NO / 1-YES]: \n\n",flagT);
                        fclose(var2);
                        fclose(var3);
        
			CgnsRoot.CgnsBase->BaseParameters->Turb->flagLimTimeScale = flagT;
			
			printf("Do you want to limit Turbulent Viscosity [0->NO / 1-YES]: \n");
			do
			{
				scanf("%d", &flagT);
			}while(!(flagT==0 || flagT==1));                        
                        	
                        var2=fopen("Settings/coefficients.dat","a");
                        var3=fopen("Settings/ProgramReadInputValues.dat","a");
                        fprintf(var2,"%d\n",flagT);
                        fprintf(var3,"%d \tDo you want to limit Turbulent Viscosity [0->NO / 1-YES]: \n\n",flagT);
                        fclose(var2);
                        fclose(var3);
        
			CgnsRoot.CgnsBase->BaseParameters->Turb->flagLimTVisc = flagT;
			
			printf("Do you want to limit Turbulent Production [0->NO / 1-YES]: \n");
			do
			{
				scanf("%d", &flagT);
			}while(!(flagT==0 || flagT==1));                        
                        	
                        var2=fopen("Settings/coefficients.dat","a");
                        var3=fopen("Settings/ProgramReadInputValues.dat","a");
                        fprintf(var2,"%d\n",flagT);
                        fprintf(var3,"%d \tDo you want to limit Turbulent Production [0->NO / 1-YES]: \n\n",flagT);
                        fclose(var2);
                        fclose(var3);
        
			CgnsRoot.CgnsBase->BaseParameters->Turb->flagLimProduction = flagT;
			
			printf("Do you want to limit k, Epsilon with minimum and maximum values [0->NO / 1-YES]: \n");
			do
			{
				scanf("%d", &flagT);
			}while(!(flagT==0 || flagT==1));                        
                        	
                        var2=fopen("Settings/coefficients.dat","a");
                        var3=fopen("Settings/ProgramReadInputValues.dat","a");
                        fprintf(var2,"%d\n",flagT);
                        fprintf(var3,"%d \tDo you want to limit k, Epsilon with minimum and maximum values [0->NO / 1-YES]: \n\n",flagT);
                        fclose(var2);
                        fclose(var3);
        
			CgnsRoot.CgnsBase->BaseParameters->Turb->flagLimKEps = flagT;
			if(flagT == 1)
			{
				printf("Enter the minimum limit for turbulent kinetic energy (k):\n\n");
				scanf("%lf", &CgnsRoot.CgnsBase->BaseParameters->Turb->kmin);                        
                                	
                                var2=fopen("Settings/coefficients.dat","a");
                                var3=fopen("Settings/ProgramReadInputValues.dat","a");
                                fprintf(var2,"%lf\n",CgnsRoot.CgnsBase->BaseParameters->Turb->kmin);
                                fprintf(var3,"%lf \tEnter the minimum limit for turbulent kinetic energy (k): \n\n",CgnsRoot.CgnsBase->BaseParameters->Turb->kmin);

				
				printf("Enter the minimum limit for dissipation rate of turbulent kinetic energy (epsilon):\n\n");
				scanf("%lf", &CgnsRoot.CgnsBase->BaseParameters->Turb->epsmin);                        
                                fprintf(var2,"%lf\n",CgnsRoot.CgnsBase->BaseParameters->Turb->epsmin);        
                                fprintf(var3,"%lf \tEnter the minimum limit for dissipation rate of turbulent kinetic energy (epsilon): \n\n",CgnsRoot.CgnsBase->BaseParameters->Turb->epsmin);
				
				printf("Enter the reference velocity for limiting the maximum turbulent kinetic energy:\n\n");
				printf("NOTE that 20 percent of mean specific kinetic energy based on this reference velocity is used as the limiting value\n"); 
				scanf("%lf", &CgnsRoot.CgnsBase->BaseParameters->Turb->Uref);                        
                                fprintf(var2,"%lf\n",CgnsRoot.CgnsBase->BaseParameters->Turb->Uref);           
                                fprintf(var3,"%lf \tEnter the reference velocity for limiting the maximum turbulent kinetic energy: \n\n",CgnsRoot.CgnsBase->BaseParameters->Turb->Uref);                             
                                
                                fclose(var2);
                                fclose(var3);
			}
		}
	} //Turbulence area, END
	}//Without Radiation, End


/*	do
	{
		printf("\nIF YOU WANT TO START WITH FRESH INITIAL CONDITION------------: [ENTER 0]");
		printf("\nELSE IF YOU WANT TO START WITH PARTIALLY CONVERGED SOLUTIONS : [ENTER 1] : ");
		scanf("%d",&restart);
	}while(!(restart == 0 || restart == 1));                        
           	
        var2=fopen("Settings/coefficients.dat","a");
	var3=fopen("Settings/ProgramReadInputValues.dat","a");
        fprintf(var2,"%d\n",restart);
        fprintf(var3,"%d \tIF YOU WANT TO START WITH FRESH INITIAL CONDITION:ENTER 0, PARTIALLY CONVERGED SOLUTIONS : ENTER:1 \n\n",restart);
        fclose(var2);
        fclose(var3);

	printf("\nCGNS Version Number  : %f\n",CgnsRoot.version);
	printf("Number of Bases      : %d\n",CgnsRoot.Number_of_bases);
	CgnsRoot.CgnsBase->BaseParameters->restart = restart;       */

	return(restart);
}

