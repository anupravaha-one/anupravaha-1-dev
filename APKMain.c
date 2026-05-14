#include "APKMain.h"

/**
This function is the root function where the program execution starts and sequentially carries out various operations
like reading grid, calculating geometrical parameters, asking for boundary conditions and finally actually solving the problem.

Upon execution this function reads the grids, using the function ReadCgnsFile() which returns the Root pointer. The Root pointer
is supplied to the function inputFunc(), which asks the user to select the solver, restart options etc. and returns an integer tag
``restart". restart is set to zero, to start the computation afresh and set to one to start the computations from a previous solution.
initializeSolidLiquid() is used to ask the user about the continuum properties of the domain (used for Conjugate Heat transfer problems).
UpdateInterfaceStruct() is called to Calculate and store the API parameters inside the inerface structure. Cell surface areas, cell volumes
and  diffusions lenghts are calculated by CalculateGeometry(). RefreshSurfaceVolumes() copies the near inerface cell volumes of the neighboring zones
into the volume slabs of the all the Zone s in a Base. Bc_input() asks and stores the problem parameters and boundary conditions for starting a
computation afresh, which in tern returns the total number of variables to be solved. Bc_inputRestart() allows to continue the computation from
a previous solution by reading the previously stored problem parameters. ReadCgnsFileRestart() reads the solution used to restart the computation
from restart.cgns. TimeDepend(), RefreshSurfaceCells() are called for u, v, w and p for updating the time dependent boundaries and other
boundary and interface slabs. massfluxes are initializes using massflux() and massfluxstar() depending on whether the solution is started
afresh or restarted from a previous solution. solverCoupled(), solverCoupled2d(), solverSemiCoupled() and  solverSegregated() are segregated using
case statements and the choice is made by the user depending on the problem. After the solution converges (for steady state problems) the solution
is written by WriteCgnsFile()*/

int main(int argc,char *argv[])
{
	FILE *nscon,*var2,*var3;
        FILE *fpwd;
        char Working_directory[1035];
	int No_of_variable,restart,vartag,nv,read=0;
 	//int nop,nophase;
	double time;
	Root CgnsRoot;

	char *filename,*filename2,*filename3;
	char filename4[] = "SOLUTIONS/initialdata.cgns";//DO NOT TAMPER WITH THIS NAME

	filename =(char *)malloc(1040*sizeof(char));
	filename2 =(char *)malloc(1040*sizeof(char));
	filename3 =(char *)malloc(1040*sizeof(char));

	char strtemporary1[1040],strtemporary2[1040],strtemporary3[1040],strtemporary4[1040];
	
	if(argc!=3)
	{
		printf("\nWrong usage...\nfollow the accompanying syntax...\n./run input_cgns_mesh_file output_cgns_file\n");
			//printf("\nWrong usage...\nfollow the accompanying syntax...\n./run GRIDS/input_file SOLUTIONS/output_file\n");
		exit(0);
	}

        printf("\n\n\t\tANUPRAVAHA-I\n");
        printf("Educational Institution Version 2018 \n");
        printf("Multi-Block Finite-Volume-Based Navier-Stokes Solver \n");

        /* Find out Working directory */
        fpwd = popen("pwd", "r");
        if (fpwd == NULL) {
          printf("Failed to run command popen\n" );
          //exit(1);
        }
        /* Read the output a line at a time - output it. */
        while (fgets(Working_directory, sizeof(Working_directory), fpwd) != NULL) {
            printf("\n\nWorking directory is = %s", Working_directory);
        }
        pclose(fpwd);
        
	var2=fopen("Settings/coefficients.dat","w");
	var3=fopen("Settings/ProgramReadInputValues.dat","w");	
        if ( var2==NULL ){  printf("\nError opening File Settings/coefficients.dat for writing \n");}
        if ( var3==NULL ){  printf("\nError opening File Settings/ProgramReadInputValues.dat for writing \n");	}
        if (( var2==NULL ) && ( var3==NULL ) ) {  
	   printf("\n--- trying to create Settings directory.... \n");	
	   system("mkdir Settings");
	   var2=fopen("Settings/coefficients.dat","w");
	   var3=fopen("Settings/ProgramReadInputValues.dat","w");	
           if (( var2==NULL ) && ( var3==NULL ) ) {  exit(1);} else { printf("...Created 'Settings' directory\n");}
	}        
        fclose(var2);
        fclose(var3);
        
	//if arguments are proper copy them into input and output filename strings.
	strcpy(strtemporary1,"");
	strcpy(strtemporary2,argv[1]);
	filename = strcat(strtemporary1,strtemporary2);
	strcpy(strtemporary3,"SOLUTIONS/");
	strcpy(strtemporary4,argv[2]);
	filename3 = strcat(strtemporary3,strtemporary4);
	strcpy(filename2,"SOLUTIONS/restart.cgns");//DO NOT TAMPER WITH THIS NAME
	
	CgnsRoot = ReadCgnsFile(filename);//calling the read function

	file_refresh(&CgnsRoot,filename,filename3);// Storing input and output filenames in the Root structure.
	
	restart = inputFunc(CgnsRoot);
	
	/**Geometry related function**/
	
	UpdateInterfaceStruct(CgnsRoot.CgnsBase); /*This func Updates Interface structure by deriving Transformation Matrix  and also assigns directions to Interfaces ** */

	find_zone_face_info(CgnsRoot.CgnsBase);//This function is added by bhanesh to improve the efficiency by removing FindInterfaceBCPtrMod
	
	initializeSolidLiquid(CgnsRoot.CgnsBase);
	
	CalculateGeometry(CgnsRoot.CgnsBase);
	
	FindMinMaxCoordinatesOfZones(CgnsRoot.CgnsBase);
	
	/**End of Geometry functions**/
	
//	InputParticle(CgnsRoot);
	
	RefreshSurfaceVolumes(CgnsRoot.CgnsBase);
	
//	printf("\nEnter 0 for GS and 1 for Multigrid\n");
//	scanf("%d",&(CgnsRoot.CgnsBase->Multigridflag));
//      var2=fopen("Settings/coefficients.dat","a");
//	var3=fopen("Settings/ProgramReadInputValues.dat","a");
//      fprintf(var2,"%d\n",CgnsRoot.CgnsBase->Multigridflag);
//      fprintf(var3,"%d \tEnter 0 for GS and 1 for Multigrid \n\n",CgnsRoot.CgnsBase->Multigridflag);
//      fclose(var2);
//      fclose(var3);

	CgnsRoot.CgnsBase->MultiGridflag = 0; 
	
	CgnsRoot.CgnsBase->masstag=0;//Added for case of Turbulent MPFV to implement boundary conidtions 
	
	if((CgnsRoot.CgnsBase->MultiGridflag)==1)
	{
		printf("Enter 0 for V cycle & 1 for W cycle\n");
		scanf("%d",&(CgnsRoot.CgnsBase->VWcycleflag));
                
                var2=fopen("Settings/coefficients.dat","a");
                var3=fopen("Settings/ProgramReadInputValues.dat","a");
                fprintf(var2,"%d\n",CgnsRoot.CgnsBase->VWcycleflag);
                fprintf(var3,"%d \tEnter 0 for V cycle & 1 for W cycle \n\n",CgnsRoot.CgnsBase->VWcycleflag);
                fclose(var2);
                fclose(var3);
        
		GridParameters1(CgnsRoot);
	}
	else // Not MultiGrid
	{
		GridParameters2(CgnsRoot);
	}
	CalcIFACE_data(CgnsRoot.CgnsBase);

	do
	{
		printf("\nIs this solution restart?");
		printf("\nNote: This will look for the solutions file written by Anupravaha in SOLUTIONS/RESTART_DATA/solution.dat");
		printf("\nIf you have solution in CGNS format - use next option");
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
	CgnsRoot.CgnsBase->BaseParameters->restart = restart;
	CgnsRoot.CgnsBase->BaseParameters->read_initialdata_cgns = 0;

	if(restart==0)
	{
		do
		{
			printf("\n\nDo you want to give initial flow field from a CGNS file [0->No/1->Yes]");//Can be opted only for single phase flow problems.
			printf("\nNOTE: If so, keep the initial flow field solution with the name \"initialdata.cgns\" in SOLUTIONS directory");
			printf("\nSOLUTIONS/restart.cgns file from previous run can be copied as initial flow field to restart a simulation : ");
			scanf("%d", &read);
		}while(!(read==0 || read==1));	
		CgnsRoot.CgnsBase->BaseParameters->read_initialdata_cgns = read;
                
                var2=fopen("Settings/coefficients.dat","a");
                var3=fopen("Settings/ProgramReadInputValues.dat","a");
                fprintf(var2,"%d\n",read);        
                fprintf(var3,"%d \tDo you want to give initial flow field from a CGNS file [0->No/1->Yes] \n\n",read);
                fclose(var2);        
                fclose(var3);
		
		nscon=fopen("Settings/converg.dat","w");
		fprintf(nscon,"0");
		fclose(nscon);
		No_of_variable=Bc_input(CgnsRoot);
		if(CgnsRoot.CgnsBase->flagTurbComb==1)//Shyam
			InitializeNonSolveCombVariables(CgnsRoot.CgnsBase);
			
		if(read == 1)
			ReadInitialData(filename4,CgnsRoot);
	}
	
	else if(restart==1)
	{
		No_of_variable=Bc_input(CgnsRoot);
		
		if(CgnsRoot.CgnsBase->flagTurbComb==1)//Shyam
			InitializeNonSolveCombVariables(CgnsRoot.CgnsBase);
		
		//ReadCgnsFileRestart(filename2,CgnsRoot);
		Read_Solution(CgnsRoot.CgnsBase);
	}

	time=0.0;

	CallMakeSlabPtrArray(CgnsRoot.CgnsBase);//Position is shifted to accommodate for MakeSlabPtrArray_com.c

	if(((CgnsRoot.CgnsBase->flagTurbModel==1)||(CgnsRoot.CgnsBase->flagTurbModel==3))&&(CgnsRoot.CgnsBase->CgnsZone->Flowsol->T!=NULL))
	{
		CalcTurbViscosity(CgnsRoot.CgnsBase);
		CalcMuEff(CgnsRoot.CgnsBase,1);
		//if(CgnsRoot.CgnsBase->flagTurbComb!=1) //For Combustion problems temperature equation is not solved but obtained from the enthalpy <<<<< Narsimha
			//CalcTurbconductivity(CgnsRoot.CgnsBase);
	}
	
	printTimeStepResidOneline(CgnsRoot.CgnsBase);

	
	time=0.0;
	for(nv=0;nv<No_of_variable;nv++)
	{
		vartag = (CgnsRoot.CgnsBase)->VARTAG[nv];
		TimeDepend(time,vartag,CgnsRoot.CgnsBase);
		RefreshSurfaceCells(vartag,CgnsRoot.CgnsBase);
	}
	//exit(12);
       
	if(CgnsRoot.CgnsBase->flagGroundwater == 0)
	{
		EvaluateVarCoefficients(CgnsRoot.CgnsBase);
		for(vartag=700;vartag<=726;vartag++)//shash
			RefreshVariableCoefficients(vartag,CgnsRoot.CgnsBase);
		
// 		if(CgnsRoot.CgnsBase->flagMPFV == 1)
// 		{	
// 			nop = CgnsRoot.CgnsBase->nphase;
// 			for(nophase=0;nophase<nop;nophase++)
// 				RefreshVariableCoefficients(720+nophase,CgnsRoot.CgnsBase);
// 		}
	}
	
	if(CgnsRoot.CgnsBase->flagSolidification==1)
		fillprevdensity(CgnsRoot.CgnsBase);
	
	if(CgnsRoot.CgnsBase->flagParticle_tracking == 0)
	{
		switch(CgnsRoot.CgnsBase->BaseParameters->flagRadiation)
		{
			
			case 0://Flow Solution with No Radiation
		
				if(restart==0)
				{
					if(CgnsRoot.CgnsBase->BaseParameters->flag2d==0)
					{
						volumeflux(CgnsRoot.CgnsBase);
						if(CgnsRoot.CgnsBase->flaginduced==1)
							Magneticflux(CgnsRoot.CgnsBase);//rajendra
					}
					else
					{
						volumeflux2d(CgnsRoot.CgnsBase);
						if(CgnsRoot.CgnsBase->flaginduced==1)
							Magneticflux2d(CgnsRoot.CgnsBase);//rajendra
					}
				}
				else if(restart==1)
				{
					if(CgnsRoot.CgnsBase->flagMPFV==0)
					{	
						RefreshSurfaceCells(100,CgnsRoot.CgnsBase);
						RefreshSurfaceCells(200,CgnsRoot.CgnsBase);
					}
					else
					{
						RefreshSurfaceCells(3100,CgnsRoot.CgnsBase);
						RefreshSurfaceCells(3200,CgnsRoot.CgnsBase);
						
					}	
					
					if((CgnsRoot.CgnsBase->BaseParameters->flag2d==0)&&(CgnsRoot.CgnsBase->flagMPFV==0))
						RefreshSurfaceCells(300,CgnsRoot.CgnsBase);
					else if((CgnsRoot.CgnsBase->BaseParameters->flag2d==0)&&(CgnsRoot.CgnsBase->flagMPFV==1))	
						RefreshSurfaceCells(3300,CgnsRoot.CgnsBase);
				}		
	
				switch(CgnsRoot.CgnsBase->flagEfield)
				{
					case 0://No Electric Field
						switch(CgnsRoot.CgnsBase->couplingTag)
						{
							case 0:
								if(CgnsRoot.CgnsBase->BaseParameters->flag2d==0)
									solverCoupled(CgnsRoot);
								else
									solverCoupled2d(CgnsRoot);
								break;
							case 1:
								if(CgnsRoot.CgnsBase->BaseParameters->flag2d==0)
									solverSemiCoupled(CgnsRoot);
								else
									solverSemiCoupled2d(CgnsRoot);
								break;
							case 2:
							 
								solverSegregated(CgnsRoot);
								break;
						}
						break;
						case 1://Both Electric Field and Flow Field
							SolverElectric(CgnsRoot);
	
							switch(CgnsRoot.CgnsBase->couplingTag)
							{
								case 0:
									if(CgnsRoot.CgnsBase->BaseParameters->flag2d==0)
										solverCoupled(CgnsRoot);
									else
										solverCoupled2d(CgnsRoot);
									break;
								case 1:
									if(CgnsRoot.CgnsBase->BaseParameters->flag2d==0)
										solverSemiCoupled(CgnsRoot);
									else
										solverSemiCoupled2d(CgnsRoot);
									break;
								case 2:
									solverSegregated(CgnsRoot);
									break;
							}
							break;
							case 2://Only Electric Field
								SolverElectric(CgnsRoot);
								break;
				}
				break;
				case 1://Only Radiation
					if(CgnsRoot.CgnsBase->BaseParameters->flag2d==0)
						RadiationSolverFiniteVolume3D(CgnsRoot);
					else
						RadiationSolverFiniteVolume2D(CgnsRoot);
	
					break;
	
		}
	}
	else
	{
		ParticleTracking(CgnsRoot);
	}
	
	//PressureForce(CgnsRoot.CgnsBase);
	//printf("I am before printsolver\n");
	
	printsolver(CgnsRoot.CgnsBase);
	//printu4bc(CgnsRoot.CgnsBase);
	//printslab(CgnsRoot.CgnsBase);
	
	WriteCgnsFile(filename3,CgnsRoot);
	
	WriteCgnsFile_vertex("SOLUTIONS/vertex.cgns",CgnsRoot);
	printf("Solution is written to file %s \n",filename3);
	return 0;
}



