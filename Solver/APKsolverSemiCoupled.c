#include "APKsolverSemiCoupled.h"
/**
The function solverSemiCoupled() is the solver for three dimensional geometries which solves all the variables
in a semicoupled approach. This function is a control function which calls different other functions in an order
as given in the semicoupled algorithm. See Project Report No.1, December 2005, Section 3.10.0.3 Page 48.
Brief about Algo:
This function calls the massvelocityoneIter() to calculate massvelocities and converges them simultaneously.
Calculates vilumeflux from the calculated value of massvelocities.
Finds out Pressure using the function PressureIter().
Repeats the above steps till Volume Flux converges.
Computes velocities by using the function CalcVariableOneIter() and converges them simultaneously.
Computes all other solvable Scalars like Temperature etc and converges simultaneusly.
Computes variables coefficients like Density, Viscosity etc. if they are userdefined functions,
using the function EvaluateVarCoefficients().
Saves the current timestep data if the user chooses so, or saves after a desired no of time steps as choosen by the user
using the function SaveTimeStepData()
Returns to massvelocity calculation till the total residual of all the solvable variables are less than their
individual convergence criteria.
If the user chooses to stop the computation at any point(by changing the value to 0 from 1 in the file "iwrite.dat" in the
directory "Settings"), then the solution at that time level is written using the function WriteCgnsFile().
*/
int solverSemiCoupled(Root CgnsRoot)
{
    int timetag,vartag,CounterTstep=0;
    int total,index_var;
    double time=0.0,dt;
    FILE *save,*odt,*Q,*fp,*fpTstep,*omega,*convcri, *var2, *var3;
    char filename2[] = "SOLUTIONS/restart.cgns"; //This value is last somewhere while calling stop - needed to re-assign later
    int timeStepToken,variableIterTag,cntVariable,stop=1;
    double variableResidual[11];
    double ErrorTimeStep[22]; //Increased the size to 16, to accomodate two combustion variables, Shyam
    //adhi+shash
    //Increased the size to 21 to accomodate the disperse phase variables//22 by shash
    double residMassFlux;
    int tempTurbTag=0, changeTimeStepTag=0;
    double turbTime=0.0,r=1.0,dummyViscosity=0.0;
    int dummyViscosityTag=0,TurbIntroTag=0,numTurbScalars = 0,initial_numTurbScalars=0;
    int fluxcount;//adhi+shash
    //for Solidification and Mleting problems BHANS
    int solidification_token=0,Iter_counter=0;
    double solidification_rms;
    //Radiative Variable
    double *Kp,Irrad,*Omega,*Int,sigma=5.67e-8,*TempP,*TempN;
    int i,j,k,l,m,ni,nj,nk,nt,np,nzone,Z_ind,INDEX,T_Bicg=400;
    int nop,nophase,numzones;
    Zone *zone;
    int massflux_iter;

    nt=CgnsRoot.CgnsBase->BaseParameters->ntheta;
    np=CgnsRoot.CgnsBase->BaseParameters->nphi;
    nzone=CgnsRoot.CgnsBase->Number_of_zones;

    FILE *countBiCg;

    countBiCg=fopen("Settings/Switch_Iter.dat","r");

    fscanf(countBiCg,"%d",&T_Bicg);

    fclose(countBiCg);


    if(CgnsRoot.CgnsBase->flagTurb == 1)
    {
        printf("Do you want to add DummyViscosity to stabilize the turbulent flow problem [0->No/1->Yes]:\n");
        do
        {
            scanf("%d", &dummyViscosityTag);
        }
        while(!(dummyViscosityTag == 0 || dummyViscosityTag == 1));

	var2=fopen("Settings/coefficients.dat","a");
	var3=fopen("Settings/ProgramReadInputValues.dat","a");
	fprintf(var2,"%d\n",dummyViscosityTag);
	fprintf(var3,"%d \tDo you want to add DummyViscosity to stabilize the turbulent flow problem [0->No/1->Yes]: \n\n",dummyViscosityTag);
	fclose(var2);
	fclose(var3);

        if(dummyViscosityTag == 1)
        {
            printf("Enter the dummyViscosity: \n");
            scanf("%lf", &dummyViscosity);

            fp=fopen("Settings/DummyViscocity.dat","w");
            fprintf(fp,"%0.9lf", dummyViscosity);
            fclose(fp);
            fp=fopen("Settings/DummyViscCofactor.dat","w");
            fprintf(fp,"0.99999");
            fclose(fp);
            printf("\n**********************************************************************************\n");
            printf("* This dummyViscosity decreases every time step with a cofactor of 0.99999.       *\nTo make it vanish more rapidly/slowly, increase/decrease the value in \"Settings/DummyViscCofactor.dat\" \n");
            printf("**********************************************************************************\n");

	    var2=fopen("Settings/coefficients.dat","a");
	    var3=fopen("Settings/ProgramReadInputValues.dat","a");
	    fprintf(var2,"%lf\n",dummyViscosity);
	    fprintf(var3,"%lf \tEnter the dummyViscosity: \n\n",dummyViscosity);
	    fclose(var2);
	    fclose(var3);
        }

        numTurbScalars = CgnsRoot.CgnsBase->BaseParameters->Turb->numTurbScalars;

        initial_numTurbScalars = numTurbScalars;

        if(CgnsRoot.CgnsBase->flagTurbModel == 6 || CgnsRoot.CgnsBase->flagTurbModel == 5 || CgnsRoot.CgnsBase->flagTurbModel == 8 || CgnsRoot.CgnsBase->flagTurbModel == 9)	//SST turbulence model
        {
            printf("\n************ Calculating Minimum Distance to Wall ***********\n");
            Distance2Wall(CgnsRoot.CgnsBase);
        }

        for(vartag=600; vartag<600+numTurbScalars; vartag++)
        {
            printf("Implementing boundary conditions for the first time, vartag=%d,\ttime=%lf\n", vartag,time);
            TimeDepend(time,vartag,CgnsRoot.CgnsBase);
            RefreshSurfaceCells(vartag,CgnsRoot.CgnsBase);
        }
        FillNewsolToPrevsol(CgnsRoot.CgnsBase);

        if(!(CgnsRoot.CgnsBase->flagTurbModel==4 || CgnsRoot.CgnsBase->flagTurbModel==5))
            CalcUtau(CgnsRoot.CgnsBase);

        if(CgnsRoot.CgnsBase->BaseParameters->Turb->flagAutoInit == 0)
        {
            CalcGamma(CgnsRoot.CgnsBase);
            CalcTurbViscosity(CgnsRoot.CgnsBase);
            freeGamma(CgnsRoot.CgnsBase);
        }
        //Added by chandrasekhar --- to caluculate diffusion coefficients for 31(Global momentum),35(Disperse phase Volume fraction),44(Disperse phase momentum)
        if(CgnsRoot.CgnsBase->flagMPFV==0)
        {
            CalcMuEff(CgnsRoot.CgnsBase,1);
            if(CgnsRoot.CgnsBase->flagParticle==1)
            {
                CalcMuEff(CgnsRoot.CgnsBase,62);//Diff coeff for disperse phase momentum in case of gas-particle flow
                CalcMuEff(CgnsRoot.CgnsBase,65);//Diff coeff for volume fraction of disperse phase in case of gas-particle flow
            }
        }
        else if(CgnsRoot.CgnsBase->flagMPFV==1)
        {
            CalcMuEff(CgnsRoot.CgnsBase,31);
            CalcMuEff(CgnsRoot.CgnsBase,35);
            CalcMuEff(CgnsRoot.CgnsBase,44);
        }
        //End --- chandrasekhar

        printf("Do you want to introduce Turbulence after few time steps?  [0->No/1->Yes]:\n");
        do
        {
            scanf("%d", &TurbIntroTag);
        }
        while(!(TurbIntroTag == 0 || TurbIntroTag == 1));

	var2=fopen("Settings/coefficients.dat","a");
	var3=fopen("Settings/ProgramReadInputValues.dat","a");
	fprintf(var2,"%d\n",TurbIntroTag);
	fprintf(var3,"%d \tDo you want to introduce Turbulence after few time steps?  [0->No/1->Yes] \n\n",TurbIntroTag);
	fclose(var2);
	fclose(var3);


        if(TurbIntroTag == 1)
        {
            printf("Enter after what TIME turbulence should be introduced?  \n");
            scanf("%lf", &turbTime);

            fpTstep=fopen("Settings/TurbulenceStartTime.dat","w");
            fprintf(fpTstep,"%0.9lf", turbTime);
            fclose(fpTstep);

	    var2=fopen("Settings/coefficients.dat","a");
	    var3=fopen("Settings/ProgramReadInputValues.dat","a");
	    fprintf(var2,"%lf\n",turbTime);
	    fprintf(var3,"%lf \tEnter after what TIME turbulence should be introduced?\n\n",turbTime);
	    fclose(var2);
	    fclose(var3);

            if(CgnsRoot.CgnsBase->flagTurb == 1)
                tempTurbTag = CgnsRoot.CgnsBase->flagTurb;	//This is used to introduce turbulence after few laminar time iterations

            CgnsRoot.CgnsBase->flagTurb = 0;
            numTurbScalars = 0;
        }

        if(CgnsRoot.CgnsBase->flagTurbComb==1)
        {
            ComputeCombustionVariables(CgnsRoot.CgnsBase);
        }
    }

    Q=fopen("SOLUTIONS/MassConservation.dat","w");
    dt=CgnsRoot.CgnsBase->BaseParameters->TimeStep;

    printf("\n I am in solverSemiCoupled \n");


    timetag = CgnsRoot.CgnsBase->BaseParameters->TimeStepingSchemeTag;

    //CgnsRoot.CgnsBase->BaseParameters->Time=0.0;

    total = CgnsRoot.CgnsBase->total_var;



    if(CgnsRoot.CgnsBase->flagGroundwater==0)
    {

        CDcoeffMemAlloc(CgnsRoot.CgnsBase);

        DiffusionCoeffCalculation(CgnsRoot.CgnsBase);

        //DiffusionCoeffCalculationForVariables(CgnsRoot.CgnsBase);

        if(CgnsRoot.CgnsBase->flagSteadyState==0)
        {
            if(CgnsRoot.CgnsBase->BaseParameters->flag_Potentialflow==1)
            {
                ///This is the solver for potential flow solution.  Giving potential flow solution as a initial conditions for the solver improves the continuity guess and pressure convergence.  If the user chooses to solve for the potential flow then the following functions will be called.
                printf("Solving for potential flow.....\n\n");

                //if(CgnsRoot.CgnsBase->BiCgflag==0)
                CalcConjugate(CgnsRoot.CgnsBase,990);
                //else //if(CgnsRoot.CgnsBase->BiCgflag==1)
                //	Conj_CalcConjugate(CgnsRoot.CgnsBase,990);
                //else if(CgnsRoot.CgnsBase->BiCgflag==2)
                //ModiCalcConjugate(CgnsRoot.CgnsBase,990);

                find_velocities(CgnsRoot.CgnsBase);
                RefreshSurfaceCells(1,CgnsRoot.CgnsBase);
                RefreshSurfaceCells(2,CgnsRoot.CgnsBase);
                RefreshSurfaceCells(3,CgnsRoot.CgnsBase);

                volumeflux(CgnsRoot.CgnsBase);

                printf("The mass flux difference is =%e\n",MassFluxConv(CgnsRoot.CgnsBase));
                //strcpy(filename2,"SOLUTIONS/Potential_Flow_Solution.cgns");
                //WriteCgnsFile(filename2,CgnsRoot);
		WriteCgnsFile("SOLUTIONS/potentialflow.cgns",CgnsRoot);

// 				if(CgnsRoot.CgnsBase->CgnsZone->Flowsol->T!=NULL && CgnsRoot.CgnsBase->flagTurbComb==0)
// 				{
// 					///For steady state problems giving pure diffusion equation solution as initial conditions for solver will improve the convergence rate.  This case really helps in case of problems where CP values are extremely high and the chosen time step is low because of the flow.
// 					printf("Solving diffusion equation for Temperature .....\n\n");
// 					CallCalcCoefficientsForVariables(CgnsRoot.CgnsBase,0,6);
// 					CalcConjugate(CgnsRoot.CgnsBase,6);
// 					WriteCgnsFile(filename2,CgnsRoot);
// 				}
                printf("The Potential flow solution with other variables is written to potentialflow.cgns\n\n");
            }
        }

        do/*time stepping loop*/
        {
            //DiffusionCoeffCalculationForVariables(CgnsRoot.CgnsBase);//Must include the header file before using this function <<< Narsimha.
#ifdef FILE_OPEN
            {
                convcri=fopen("Settings/convergenceCriteria.dat","r");
                for(cntVariable=0; cntVariable<12; cntVariable++)
                    fscanf(convcri,"%lf",(CgnsRoot.CgnsBase->BaseParameters->convergenceCriterion+cntVariable));
                fclose(convcri);

                omega=fopen("Settings/omega.dat","r");
                for(cntVariable=0; cntVariable<14; cntVariable++)
                    fscanf(omega,"%lf",(CgnsRoot.CgnsBase->BaseParameters->OMEGA+cntVariable));
                fclose(omega);

// 			fp=fopen("SOLUTIONS/No_of_TimeSteps","w");
// 			fprintf(fp,"%d\n%lf", CounterTstep,time);
// 			fclose(fp);

                save=fopen("Settings/iwrite.dat","r");
                fscanf(save,"%d",&stop);
                fclose(save);

                if(dummyViscosityTag == 1)
                {
                    fp=fopen("Settings/DummyViscCofactor.dat","r");
                    fscanf(fp,"%lf", &r);
                    fclose(fp);

                    fp=fopen("Settings/DummyViscocity.dat","r");
                    fscanf(fp,"%lf", &dummyViscosity);
                    fclose(fp);

                    dummyViscosity = r * dummyViscosity;

                    fp=fopen("Settings/DummyViscocity.dat","w");
                    fprintf(fp,"%0.9lf", dummyViscosity);
                    fclose(fp);
                }

                fpTstep = fopen("Settings/ChangeTimeStep.dat","r");
                fscanf(fpTstep,"%d", &changeTimeStepTag);
                if(changeTimeStepTag == 1)
                {
                    fscanf(fpTstep,"%lf", &dt);
                    printf("\n ******** Changing the Time step to %lf ******** \n", dt);
                    CgnsRoot.CgnsBase->BaseParameters->TimeStep = dt;
                    CgnsRoot.CgnsBase->BaseParameters->Turb->DualTimeStep = dt/CgnsRoot.CgnsBase->BaseParameters->Turb->numDualTimeIter;
                }
                fclose(fpTstep);
            }
#endif
            if(tempTurbTag == 1)
            {
#ifdef FILE_OPEN
                {
                    fpTstep=fopen("Settings/TurbulenceStartTime.dat","r");
                    fscanf(fpTstep,"%lf", &turbTime);
                    fclose(fpTstep);
                }
#endif
                if(time >=turbTime)
                {
                    CgnsRoot.CgnsBase->flagTurb = 1;
                    tempTurbTag = 0;
                    numTurbScalars = initial_numTurbScalars;

                    printf("######## Turbulence is introduced after Time = %lf #########\n", time);
                    for(vartag=600; vartag<600+numTurbScalars; vartag++)
                    {
                        printf("Implementing boundary conditions for Scalars 1 and 2 , vartag=%d,\ttime=%lf\n", vartag,time);
                        TimeDepend(time,vartag,CgnsRoot.CgnsBase);
                        RefreshSurfaceCells(vartag,CgnsRoot.CgnsBase);
                    }
                    FillNewsolToPrevsolKepsilon(CgnsRoot.CgnsBase);
                    if(!(CgnsRoot.CgnsBase->flagTurbModel==4 || CgnsRoot.CgnsBase->flagTurbModel==5))
                        CalcUtau(CgnsRoot.CgnsBase);
                    CalcTurbViscosity(CgnsRoot.CgnsBase);
                    CalcGamma(CgnsRoot.CgnsBase);
                    CalcMuEff(CgnsRoot.CgnsBase,1);
                    freeGamma(CgnsRoot.CgnsBase);
                }
            }

            if (timetag==1 && time==0) /*Even the scheme is CN, the first timestep is done by BTCS scheme*/
                CgnsRoot.CgnsBase->BaseParameters->TimeStepingSchemeTag=0;
            else if (timetag==1 && time>0)
                CgnsRoot.CgnsBase->BaseParameters->TimeStepingSchemeTag=1;

            time = time+dt;
            CounterTstep=CounterTstep+1;


            if(stop==0)
            {
		printf("Writing the solution data ...\n");
		strcpy(filename2,"SOLUTIONS/restart.cgns");
		//printf("CGNS File Name : %s\n",filename2);
                WriteCgnsFile(filename2,CgnsRoot);
                WriteCgnsFile_vertex("SOLUTIONS/vertex.cgns",CgnsRoot);
                printsolver(CgnsRoot.CgnsBase);
                dt=CgnsRoot.CgnsBase->BaseParameters->TimeStep;
                odt=fopen("Settings/old_dt.dat","w");
                fprintf(odt,"%12.10lf %12.10lf",dt,CgnsRoot.CgnsBase->BaseParameters->Time);
                fclose(odt);
                printslab(CgnsRoot.CgnsBase);
                //printu4bc(CgnsRoot.CgnsBase);

                Write_Solution(CgnsRoot.CgnsBase,"SOLUTIONS/RESTART_DATA/solution.dat");
		//Write_gnuplot_style_solution(CgnsRoot.CgnsBase,"SOLUTIONS/gnuplot/gnu_solution.dat");

                exit(1);
            }

            CgnsRoot.CgnsBase->BaseParameters->Time=CgnsRoot.CgnsBase->BaseParameters->Time+dt;

            if(CounterTstep%1 == 0)
                printf("\nNo of Time Steps: %d\tTime: %lf\n", CounterTstep,CgnsRoot.CgnsBase->BaseParameters->Time);//<----This has been shifted to sovlerCoupled.c

            FillNewsolToPrevsol(CgnsRoot.CgnsBase);

            if(CgnsRoot.CgnsBase->flagLMMHD==1)//rajendra
            {
                if(CgnsRoot.CgnsBase->flaginduced==1)
                {
                    InducedLMMHD_Source_Calculation(CgnsRoot.CgnsBase);
                }
                else
                {

                    LMMHDmass_Source_Calculation(CgnsRoot.CgnsBase);
                    //  LMMHD_Source_Calculation(CgnsRoot.CgnsBase);
                }

            }


            if((CgnsRoot.CgnsBase->flagParticle==1)|| (CgnsRoot.CgnsBase->flagDroplet==1))
            {

                if((CgnsRoot.CgnsBase->BaseParameters->flagVariableDia) == 1)
                {
                    //  CalcVarDropDia(CgnsRoot.CgnsBase);
                    CalcEvapoConst(CgnsRoot.CgnsBase);
                    //  Evaporate(CgnsRoot.CgnsBase);
                }

                volumefluxDisperse2d(CgnsRoot.CgnsBase);

                Drag(CgnsRoot.CgnsBase);
            }
            if(CgnsRoot.CgnsBase->flagMPFV==0)
                CallCalcCoefficientsForVariables(CgnsRoot.CgnsBase,0,1);
            else
            {
                CallCalcCoefficientsForVariables(CgnsRoot.CgnsBase,0,31);
                CallCalcCoefficientsForVariables(CgnsRoot.CgnsBase,0,4);
            }
            massflux_iter=0;
            if(CgnsRoot.CgnsBase->SolverScheme==0)//Fully - Implicit Solver
            {
                do/*mass flux loop*/
                {
                    savevolumeflux(CgnsRoot.CgnsBase);
                    if(CgnsRoot.CgnsBase->flagMPFV==1)
                    {
                        Conv_Coeff_Calculation(CgnsRoot.CgnsBase,0,31);

                        CDcoeffCalculation(CgnsRoot.CgnsBase,31);

                        CalcSource_opti(CgnsRoot.CgnsBase,3100);
                        massVelocity(CgnsRoot.CgnsBase, 31);

                        CalcSource_opti(CgnsRoot.CgnsBase,3200);
                        massVelocity(CgnsRoot.CgnsBase, 32);

                        CalcSource_opti(CgnsRoot.CgnsBase,3300);
                        massVelocity(CgnsRoot.CgnsBase, 33);

                        volumefluxstar(CgnsRoot.CgnsBase);

                        if(CgnsRoot.CgnsBase->BiCgflag==0)
                            PressureIter(CgnsRoot.CgnsBase);//Calculate new pressure from new Massflux
                        else if(CgnsRoot.CgnsBase->BiCgflag==1)
                            Conj_PressureIter(CgnsRoot.CgnsBase);
                        else if(CgnsRoot.CgnsBase->BiCgflag==2)
                            ModiPressureIter(CgnsRoot.CgnsBase);

                        CalcNewVolumeFlux(CgnsRoot.CgnsBase);

                        residMassFlux=FindVolumeFluxDifference(CgnsRoot.CgnsBase);
                    }
                    else if(CgnsRoot.CgnsBase->flagDroplet==1 && CgnsRoot.CgnsBase->BaseParameters->flagVariableDia==1)
                    {
                        puts("\nInside Fully Implicit Solver for Evaporating Droplets");

                        savevolumeflux(CgnsRoot.CgnsBase);

                        Conv_Coeff_Calculation(CgnsRoot.CgnsBase,0,1);
                        CDcoeffCalculation(CgnsRoot.CgnsBase,1);

                        CalcSource_opti(CgnsRoot.CgnsBase,100);
                        massVelocity(CgnsRoot.CgnsBase, 1);

                        CalcSource_opti(CgnsRoot.CgnsBase,200);
                        massVelocity(CgnsRoot.CgnsBase, 2);

                        CalcSource_opti(CgnsRoot.CgnsBase,300);
                        massVelocity(CgnsRoot.CgnsBase, 3);

                        volumefluxstar(CgnsRoot.CgnsBase);

                        Conv_Coeff_Calculation(CgnsRoot.CgnsBase,0,62);
                        CDcoeffCalculation(CgnsRoot.CgnsBase,62);

                        CalcSource_opti(CgnsRoot.CgnsBase,6200);
                        massVelocity(CgnsRoot.CgnsBase, 62);

                        CalcSource_opti(CgnsRoot.CgnsBase,6300);
                        massVelocity(CgnsRoot.CgnsBase, 63);

                        CalcSource_opti(CgnsRoot.CgnsBase,6400);
                        massVelocity(CgnsRoot.CgnsBase, 64);

                        DisperseVolumeFluxStar(CgnsRoot.CgnsBase);

                        if(CgnsRoot.CgnsBase->BiCgflag==0)
                            PressureIter(CgnsRoot.CgnsBase);//Calculate new pressure from new Massflux
                        else if(CgnsRoot.CgnsBase->BiCgflag==1)
                            Conj_PressureIter(CgnsRoot.CgnsBase);
                        else if(CgnsRoot.CgnsBase->BiCgflag==2)
                            ModiPressureIter(CgnsRoot.CgnsBase);


                        CalcNewVolumeFlux(CgnsRoot.CgnsBase);
                        CalcNewVolumeFluxDisperse(CgnsRoot.CgnsBase);

                        residMassFlux=FindDisperseVolumeFluxDifference(CgnsRoot.CgnsBase);
                    }
                    else
                    {
                        Conv_Coeff_Calculation(CgnsRoot.CgnsBase,0,1);
                        CDcoeffCalculation(CgnsRoot.CgnsBase,1);

                        CalcSource_opti(CgnsRoot.CgnsBase,100);
                        massVelocity(CgnsRoot.CgnsBase, 1);

                        CalcSource_opti(CgnsRoot.CgnsBase,200);
                        massVelocity(CgnsRoot.CgnsBase, 2);

                        CalcSource_opti(CgnsRoot.CgnsBase,300);
                        massVelocity(CgnsRoot.CgnsBase, 3);

                        volumefluxstar(CgnsRoot.CgnsBase);
                        residMassFlux=FindVolumeFluxDifference(CgnsRoot.CgnsBase);
                    }
                    //if(residMassFlux > 1.0e-03)
                    massflux_iter++;
                    if(massflux_iter%50==0)
                        printf("mass flux residual =%e\n",residMassFlux);
                }
                while(residMassFlux > CgnsRoot.CgnsBase->BaseParameters->MassFluxConv);


                if(CgnsRoot.CgnsBase->BiCgflag==0)
                    PressureIter(CgnsRoot.CgnsBase);//Calculate new pressure from new Massflux
                else if(CgnsRoot.CgnsBase->BiCgflag==1)
                    Conj_PressureIter(CgnsRoot.CgnsBase);
                else if(CgnsRoot.CgnsBase->BiCgflag==2)
                    ModiPressureIter(CgnsRoot.CgnsBase);

                CalcNewVolumeFlux(CgnsRoot.CgnsBase);



            }
            else if(CgnsRoot.CgnsBase->SolverScheme==1)//Semi - Implicit Solver
            {
                if(CounterTstep>1)
                {
                    if(CgnsRoot.CgnsBase->flagDroplet==1 && CgnsRoot.CgnsBase->BaseParameters->flagVariableDia ==1)
					{
							savevolumeflux(CgnsRoot.CgnsBase);

							Conv_Coeff_Calculation(CgnsRoot.CgnsBase,0,1);
							CDcoeffCalculation(CgnsRoot.CgnsBase,1);

							CalcSource_opti(CgnsRoot.CgnsBase,100);
							massVelocity(CgnsRoot.CgnsBase, 1);

							CalcSource_opti(CgnsRoot.CgnsBase,200);
							massVelocity(CgnsRoot.CgnsBase, 2);

							CalcSource_opti(CgnsRoot.CgnsBase,300);
							massVelocity(CgnsRoot.CgnsBase, 3);

							volumefluxstar(CgnsRoot.CgnsBase);

							Conv_Coeff_Calculation(CgnsRoot.CgnsBase,0,62);
							CDcoeffCalculation(CgnsRoot.CgnsBase,62);

							CalcSource_opti(CgnsRoot.CgnsBase,6200);
							massVelocity(CgnsRoot.CgnsBase, 62);

							CalcSource_opti(CgnsRoot.CgnsBase,6300);
							massVelocity(CgnsRoot.CgnsBase, 63);

							CalcSource_opti(CgnsRoot.CgnsBase,6400);
							massVelocity(CgnsRoot.CgnsBase, 64);

							DisperseVolumeFluxStar(CgnsRoot.CgnsBase);

							if(CgnsRoot.CgnsBase->BiCgflag==0)
								PressureIter(CgnsRoot.CgnsBase);//Calculate new pressure from new Massflux
							else if(CgnsRoot.CgnsBase->BiCgflag==1)
								Conj_PressureIter(CgnsRoot.CgnsBase);
							else if(CgnsRoot.CgnsBase->BiCgflag==2)
								ModiPressureIter(CgnsRoot.CgnsBase);

							CalcNewVolumeFlux(CgnsRoot.CgnsBase);
							CalcNewVolumeFluxDisperse(CgnsRoot.CgnsBase);
					}
                    else
                    {
                        savevolumeflux(CgnsRoot.CgnsBase);

                        Conv_Coeff_Calculation(CgnsRoot.CgnsBase,0,1);
                        CDcoeffCalculation(CgnsRoot.CgnsBase,1);

                        CalcSource_opti(CgnsRoot.CgnsBase,100);
                        massVelocity(CgnsRoot.CgnsBase, 1);

                        CalcSource_opti(CgnsRoot.CgnsBase,200);
                        massVelocity(CgnsRoot.CgnsBase, 2);

                        CalcSource_opti(CgnsRoot.CgnsBase,300);
                        massVelocity(CgnsRoot.CgnsBase, 3);

                        volumefluxstar(CgnsRoot.CgnsBase);

                        if(CgnsRoot.CgnsBase->BiCgflag==0)
                            PressureIter(CgnsRoot.CgnsBase);//Calculate new pressure from new Massflux
                        else if(CgnsRoot.CgnsBase->BiCgflag==1)
                            Conj_PressureIter(CgnsRoot.CgnsBase);
                        else if(CgnsRoot.CgnsBase->BiCgflag==2)
                            ModiPressureIter(CgnsRoot.CgnsBase);

                        CalcNewVolumeFlux(CgnsRoot.CgnsBase);
                    }
                }
                else
                {
                    do/*mass flux loop*/
                    {
                        if(CgnsRoot.CgnsBase->flagDroplet==1 && CgnsRoot.CgnsBase->BaseParameters->flagVariableDia ==1)
						{
							savevolumeflux(CgnsRoot.CgnsBase);

							Conv_Coeff_Calculation(CgnsRoot.CgnsBase,0,1);
							CDcoeffCalculation(CgnsRoot.CgnsBase,1);

							CalcSource_opti(CgnsRoot.CgnsBase,100);
							massVelocity(CgnsRoot.CgnsBase, 1);

							CalcSource_opti(CgnsRoot.CgnsBase,200);
							massVelocity(CgnsRoot.CgnsBase, 2);

							CalcSource_opti(CgnsRoot.CgnsBase,300);
							massVelocity(CgnsRoot.CgnsBase, 3);

							volumefluxstar(CgnsRoot.CgnsBase);

							Conv_Coeff_Calculation(CgnsRoot.CgnsBase,0,62);
							CDcoeffCalculation(CgnsRoot.CgnsBase,62);

							CalcSource_opti(CgnsRoot.CgnsBase,6200);
							massVelocity(CgnsRoot.CgnsBase, 62);

							CalcSource_opti(CgnsRoot.CgnsBase,6300);
							massVelocity(CgnsRoot.CgnsBase, 63);

							CalcSource_opti(CgnsRoot.CgnsBase,6400);
							massVelocity(CgnsRoot.CgnsBase, 64);

							DisperseVolumeFluxStar(CgnsRoot.CgnsBase);

							if(CgnsRoot.CgnsBase->BiCgflag==0)
								PressureIter(CgnsRoot.CgnsBase);//Calculate new pressure from new Massflux
							else if(CgnsRoot.CgnsBase->BiCgflag==1)
								Conj_PressureIter(CgnsRoot.CgnsBase);
							else if(CgnsRoot.CgnsBase->BiCgflag==2)
								ModiPressureIter(CgnsRoot.CgnsBase);

							CalcNewVolumeFlux(CgnsRoot.CgnsBase);
							CalcNewVolumeFluxDisperse(CgnsRoot.CgnsBase);

							residMassFlux=FindDisperseVolumeFluxDifference(CgnsRoot.CgnsBase);

							massflux_iter++;

							if(massflux_iter%50==0)
							{
								printf("\nvolume flux residual of droplet solver=%e\n",residMassFlux);
							}
						}

                        else
                        {
                        savevolumeflux(CgnsRoot.CgnsBase);

                        Conv_Coeff_Calculation(CgnsRoot.CgnsBase,0,1);
                        CDcoeffCalculation(CgnsRoot.CgnsBase,1);

                        CalcSource_opti(CgnsRoot.CgnsBase,100);
                        massVelocity(CgnsRoot.CgnsBase, 1);

                        CalcSource_opti(CgnsRoot.CgnsBase,200);
                        massVelocity(CgnsRoot.CgnsBase, 2);

                        CalcSource_opti(CgnsRoot.CgnsBase,300);
                        massVelocity(CgnsRoot.CgnsBase, 3);

                        volumefluxstar(CgnsRoot.CgnsBase);

                        if(CgnsRoot.CgnsBase->BiCgflag==0)
                            PressureIter(CgnsRoot.CgnsBase);//Calculate new pressure from new Massflux
                        else if(CgnsRoot.CgnsBase->BiCgflag==1)
                            Conj_PressureIter(CgnsRoot.CgnsBase);
                        else if(CgnsRoot.CgnsBase->BiCgflag==2)
                            ModiPressureIter(CgnsRoot.CgnsBase);

                        CalcNewVolumeFlux(CgnsRoot.CgnsBase);

                        residMassFlux=FindVolumeFluxDifference(CgnsRoot.CgnsBase);
                        massflux_iter++;
                        if(massflux_iter%50==0)
                            printf("mass flux residual =%e\n",residMassFlux);

                        }
                    }
                    while(residMassFlux > CgnsRoot.CgnsBase->BaseParameters->MassFluxConv);
                }
            }


            if(CgnsRoot.CgnsBase->flagMPFV==0)
            {
                Conv_Coeff_Calculation(CgnsRoot.CgnsBase,0,1);
                CDcoeffCalculation(CgnsRoot.CgnsBase,1);

                Pressure_Grad_Calculation(CgnsRoot.CgnsBase);//For Optimizing testing purpose..  Can be removed later.

                CalcSource_opti(CgnsRoot.CgnsBase,1);
                CalcVariable(CgnsRoot.CgnsBase,1);

                CalcSource_opti(CgnsRoot.CgnsBase,2);
                CalcVariable(CgnsRoot.CgnsBase,2);

                CalcSource_opti(CgnsRoot.CgnsBase,3);
                CalcVariable(CgnsRoot.CgnsBase,3);
            }
            else if(CgnsRoot.CgnsBase->flagMPFV==1)
            {
                nop = CgnsRoot.CgnsBase->nphase;
                Conv_Coeff_Calculation(CgnsRoot.CgnsBase,0,31);
                CDcoeffCalculation(CgnsRoot.CgnsBase,31);

                Pressure_Grad_Calculation(CgnsRoot.CgnsBase);//For Optimizing testing purpose..  Can be removed later.

                CalcSource_opti(CgnsRoot.CgnsBase,31);
                CalcVariable(CgnsRoot.CgnsBase,31);

                CalcSource_opti(CgnsRoot.CgnsBase,32);
                CalcVariable(CgnsRoot.CgnsBase,32);

                CalcSource_opti(CgnsRoot.CgnsBase,33);
                CalcVariable(CgnsRoot.CgnsBase,33);

                //Calculation of individual phase
                if(CgnsRoot.CgnsBase->BaseParameters->Dragflag==1)
                    DragMPFV(CgnsRoot.CgnsBase,CounterTstep);
                if((CgnsRoot.CgnsBase->BaseParameters->Dragflag==1)&&(CgnsRoot.CgnsBase->flagTurb == 1)&&(CgnsRoot.CgnsBase->flagLift == 1))
                    LiftMPFV(CgnsRoot.CgnsBase,CounterTstep);
                for(nophase=1; nophase<nop; nophase++) //This is done for solving disperse phase only.
                {
                    CallCalcCoefficientsForVariables(CgnsRoot.CgnsBase,0,41+3*nophase);

                    for(numzones=0; numzones<nzone; numzones++)
                    {
                        zone = CgnsRoot.CgnsBase->CgnsZone+numzones;
                        UpdateVarMPFVDiffusion(CgnsRoot.CgnsBase,zone,41+3*nophase);
                        UpdateVarMPFVDiffusion(CgnsRoot.CgnsBase,zone,41+3*nophase+1);
                        UpdateVarMPFVDiffusion(CgnsRoot.CgnsBase,zone,41+3*nophase+2);
                    }

                    CDcoeffCalculation(CgnsRoot.CgnsBase,41+3*nophase);


                    CalcSource_opti(CgnsRoot.CgnsBase,41+3*nophase);
                    CalcVariable(CgnsRoot.CgnsBase,41+3*nophase);

                    CalcSource_opti(CgnsRoot.CgnsBase,42+3*nophase);
                    CalcVariable(CgnsRoot.CgnsBase,42+3*nophase);

                    CalcSource_opti(CgnsRoot.CgnsBase,43+3*nophase);
                    CalcVariable(CgnsRoot.CgnsBase,43+3*nophase);

                    //Volume fraction calculation
                    if(CgnsRoot.CgnsBase->flagTurb==0)
                    {
                        CalcSource_opti(CgnsRoot.CgnsBase,34+nophase);//This is only for laminar flow.
                        CalcVariable(CgnsRoot.CgnsBase,34+nophase);
                    }
                    else if(CgnsRoot.CgnsBase->flagTurb==1)//Added by chandrasekhar
                    {
                        if(CgnsRoot.CgnsBase->flagDispersion==0)
                        {
                            CalcSource_opti(CgnsRoot.CgnsBase,34+nophase);
                            CalcVariable(CgnsRoot.CgnsBase,34+nophase);
                        }
                        else
                        {
                            volumefluxMPFV(CgnsRoot.CgnsBase);
                            CallCalcCoefficientsForVariables(CgnsRoot.CgnsBase,0,34+nophase);
                            CDcoeffCalculation(CgnsRoot.CgnsBase,34+nophase);

                            CalcSource_opti(CgnsRoot.CgnsBase,34+nophase);
                            CalcVariable(CgnsRoot.CgnsBase,34+nophase);
                        }
                    }
                }

                MPFV(CgnsRoot.CgnsBase);
            }

            //disperse phase....adhi+shash
			if((CgnsRoot.CgnsBase->flagParticle==1)||(CgnsRoot.CgnsBase->flagDroplet==1))
			{
				if(CgnsRoot.CgnsBase->BaseParameters->flagVariableDia==1)
				{

					Pressure_Grad_Calculation(CgnsRoot.CgnsBase);

					Conv_Coeff_Calculation(CgnsRoot.CgnsBase,0,62);
					CDcoeffCalculation(CgnsRoot.CgnsBase,62);

					CalcSource_opti(CgnsRoot.CgnsBase,62);
					CalcVariable(CgnsRoot.CgnsBase,62);

					CalcSource_opti(CgnsRoot.CgnsBase,63);
					CalcVariable(CgnsRoot.CgnsBase,63);

					CalcSource_opti(CgnsRoot.CgnsBase,64);
					CalcVariable(CgnsRoot.CgnsBase,64);
				}
				else
				{
					fluxcount=0;
					do
					{
						savevolumeflux(CgnsRoot.CgnsBase);

						Conv_Coeff_Calculation(CgnsRoot.CgnsBase,0,62);

						if(CgnsRoot.CgnsBase->flagTurb == 1)
						{
							CallCalcCoefficientsForVariables(CgnsRoot.CgnsBase,0,62);
							for(numzones=0;numzones<nzone;numzones++)
							{
								zone = CgnsRoot.CgnsBase->CgnsZone+numzones;

								UpdateVarDisperseDiffusion(CgnsRoot.CgnsBase,zone,62);

								UpdateVarDisperseDiffusion(CgnsRoot.CgnsBase,zone,63);

								UpdateVarDisperseDiffusion(CgnsRoot.CgnsBase,zone,64);
							}
						}

						CDcoeffCalculation(CgnsRoot.CgnsBase,62);

						CalcSource_opti(CgnsRoot.CgnsBase,62);
						CalcVariable(CgnsRoot.CgnsBase,62);


						CalcSource_opti(CgnsRoot.CgnsBase,63);
						CalcVariable(CgnsRoot.CgnsBase,63);

						CalcSource_opti(CgnsRoot.CgnsBase,64);
						CalcVariable(CgnsRoot.CgnsBase,64);

						volumefluxDisperse(CgnsRoot.CgnsBase);

						residMassFlux=FindMassFluxDisperseDifference(CgnsRoot.CgnsBase);
						if(fluxcount%100==0)
						{
						printf("\nresid mass flux=%12.14lf\n",residMassFlux);
						}

						fluxcount=fluxcount+1;
					}
					while(residMassFlux > CgnsRoot.CgnsBase->BaseParameters->MassFluxConv);
				}

					savevolumeflux(CgnsRoot.CgnsBase);

					Conv_Coeff_Calculation(CgnsRoot.CgnsBase,0,65);

					if(CgnsRoot.CgnsBase->flagTurb == 1)
					{
						CallCalcCoefficientsForVariables(CgnsRoot.CgnsBase,0,65);
					}

					CDcoeffCalculation(CgnsRoot.CgnsBase,65);

					CalcSource_opti(CgnsRoot.CgnsBase,65);
					CalcVariable(CgnsRoot.CgnsBase,65);

				if((CgnsRoot.CgnsBase->flagDroplet==1)&&(CgnsRoot.CgnsBase->BaseParameters->flagVariableDia==1))
				{
	 				Conv_Coeff_Calculation(CgnsRoot.CgnsBase,0,710);
					CDcoeffCalculation(CgnsRoot.CgnsBase,710);
					CalcSource_opti(CgnsRoot.CgnsBase,710);
					CalcVariable(CgnsRoot.CgnsBase,710);
				}
				/*
				if((CgnsRoot.CgnsBase->flagDroplet==1)||(CgnsRoot.CgnsBase->BaseParameters->flagVariableDia == 1))
				{
					savevolumeflux(CgnsRoot.CgnsBase);
					Conv_Coeff_Calculation(CgnsRoot.CgnsBase,0,67);

					CDcoeffCalculation(CgnsRoot.CgnsBase,67);

					CalcSource_opti(CgnsRoot.CgnsBase,67);
					CalcVariable(CgnsRoot.CgnsBase,67);
				}
			*/
				if((CgnsRoot.CgnsBase->flagDroplet==1)&&(CgnsRoot.CgnsBase->BaseParameters->flagDropFuelMassFrac == 1))
				{
					CallCalcCoefficientsForVariables(CgnsRoot.CgnsBase,0,68);
					Conv_Coeff_Calculation(CgnsRoot.CgnsBase,0,68);
					CDcoeffCalculation(CgnsRoot.CgnsBase,68);

					CalcSource_opti(CgnsRoot.CgnsBase,68);

					CalcVariable(CgnsRoot.CgnsBase,68);
				}
				if((CgnsRoot.CgnsBase->flagDroplet==1)||(CgnsRoot.CgnsBase->BaseParameters->flagVariableDia == 1))
				{
					savevolumeflux(CgnsRoot.CgnsBase);
			        	Conv_Coeff_Calculation(CgnsRoot.CgnsBase,0,69);

					CDcoeffCalculation(CgnsRoot.CgnsBase,69);

					CalcSource_opti(CgnsRoot.CgnsBase,69);
					CalcVariable(CgnsRoot.CgnsBase,69);
			        }
			}
            //##############################


// 			if(CgnsRoot.CgnsBase->CgnsZone->Flowsol->T!=NULL && CgnsRoot.CgnsBase->flagTurbComb==0)
// 			{
// 				CallCalcCoefficientsForVariables(CgnsRoot.CgnsBase,0,6);
// 				CDcoeffCalculation(CgnsRoot.CgnsBase,6);
// 			}

// 			for(iScal=0;iScal<10;iScal++)//What this loop does?? <<< Narsimha
// 			{
// 				if(CgnsRoot.CgnsBase->CgnsZone->Flowsol->Scal[iScal]!=NULL)
// 				{
// 					CallCalcCoefficientsForVariables(CgnsRoot.CgnsBase,0,600+iScal);
// 					CDcoeffCalculation(CgnsRoot.CgnsBase,600+iScal);
// 				}
// 			}

            SolveTurbVariables(CgnsRoot,variableResidual,CounterTstep);

            if(CgnsRoot.CgnsBase->flagTurbComb==1)//Shyam
            {
                CalcMuEff(CgnsRoot.CgnsBase,14);
                //Seems this statement should be used here (may be inside SolveCombustionVariables, but check which vartag is to be supplied <<< Narsimha 18th April, 2008.
                //CallCalcCoefficientsForVariables(base,0,14);
                //CDcoeffCalculation(CgnsRoot.CgnsBase,14);
                SolveCombustionVariables(CgnsRoot.CgnsBase);
                ComputeCombustionVariables(CgnsRoot.CgnsBase);
            }

            if(CgnsRoot.CgnsBase->flagLMMHD==1)//naveen
            {
                solverLMMHD(CgnsRoot.CgnsBase);//rajendra

                CalcCurrentDensity(CgnsRoot.CgnsBase);
            }

            if((CgnsRoot.CgnsBase->flagSolidification==1)||(CgnsRoot.CgnsBase->flagTwoPhase==1))
            {
                ErrorTimeStep[3]=0.0;
                solidification_token=0;
                Iter_counter=0;

                CallCalcCoefficientsForVariables(CgnsRoot.CgnsBase,0,7);
                CDcoeffCalculation(CgnsRoot.CgnsBase,7);
                //CalcSource_opti(CgnsRoot.CgnsBase,7);
                ///CalcSource_opti will be called inside of CalcVariable() for variableTag '7' and melting and solidification problems.
                do
                {
                    CalcVariable(CgnsRoot.CgnsBase,7);
                    //variableResidual[0] = CalcVariableOneIter(CgnsRoot.CgnsBase,7);
                    solidification_rms = CalculateErrorTemp(CgnsRoot.CgnsBase,7);
                    if (fabs(solidification_rms-ErrorTimeStep[3]) <=1e-07)
                        solidification_token=1;
                    else
                        ErrorTimeStep[3] = solidification_rms;
                    Iter_counter++;
                }
                while(solidification_token==0);
            }
            else
                ErrorTimeStep[3]=0.0;

            if(CgnsRoot.CgnsBase->BaseParameters->flagflowRadiation==1)
            {
                //if(time>0.01)
                RadiationSolverWithFlow3D(CgnsRoot.CgnsBase);

                for(Z_ind=0; Z_ind<nzone; Z_ind++)
                {
                    Int=(CgnsRoot.CgnsBase->CgnsZone+Z_ind)->Flowsol->RadIntensity;
                    Omega=(CgnsRoot.CgnsBase->CgnsZone+Z_ind)->ZGeom->Omega;
                    ni=*((CgnsRoot.CgnsBase->CgnsZone+Z_ind)->Size+0);
                    nj=*((CgnsRoot.CgnsBase->CgnsZone+Z_ind)->Size+1);
                    nk=*((CgnsRoot.CgnsBase->CgnsZone+Z_ind)->Size+2);
                    TempP=(CgnsRoot.CgnsBase->CgnsZone+Z_ind)->PrevFlowsol->T;
                    TempN=(CgnsRoot.CgnsBase->CgnsZone+Z_ind)->Flowsol->T;

                    Kp=(CgnsRoot.CgnsBase->CgnsZone+Z_ind)->Flowsol->kappa;

                    for(k=0; k<nk-1; k++)
                        for(j=0; j<nj-1; j++)
                            for(i=0; i<ni-1; i++)
                            {
                                INDEX=(ni-1)*(nj-1)*k+(ni-1)*j+i;
                                Irrad=0.0;
                                for(l=0; l<nt; l++)
                                    for(m=0; m<np; m++)
                                        Irrad=Irrad+*(Int+((ni-1)*(nj-1)*k+(ni-1)*j+i)*nt*np+np*l+m)*(*(Omega+np*l+m));
                                *((CgnsRoot.CgnsBase->CgnsZone+Z_ind)->Flowsol->RadEnergysource+INDEX)=*Kp*(4.0*sigma*pow(*(TempP+INDEX),3)*(*(TempN+INDEX))-Irrad);
                                //printf("%e\t",*((CgnsRoot.CgnsBase->CgnsZone+Z_ind)->Flowsol->RadEnergysource+INDEX));
                            }
                }
            }

            if(CgnsRoot.CgnsBase->flagTurbComb==0) //Shyam
            {
 		//printf("Solving Temperature equation within the time step\n");
                if(CgnsRoot.CgnsBase->CgnsZone->Flowsol->T!=NULL)
                {
                    if(CgnsRoot.CgnsBase->flagTurb == 1)
                    {
                        CalcTurbconductivity(CgnsRoot.CgnsBase);
                    }

                    EvaluateEnergysource(CgnsRoot.CgnsBase);

                    CallCalcCoefficientsForVariables(CgnsRoot.CgnsBase,0,6);
                    CDcoeffCalculation(CgnsRoot.CgnsBase,6);

                    CalcSource_opti(CgnsRoot.CgnsBase,6);

                    if(CgnsRoot.CgnsBase->flagConjugateProblem==0 && CgnsRoot.CgnsBase->flagSteadyState==0)/*Conjugate Problem and Steady State*/
                        variableResidual[0] = CalcVariableConjugate(CgnsRoot.CgnsBase,6);
                    else
                        variableResidual[0] = CalcVariable(CgnsRoot.CgnsBase,6);

                    //CalcNusseltNo(CgnsRoot.CgnsBase,6);  //Seems to give segmentaion fault for differentially heated cube test case
                }
                else
                    variableResidual[0]=0.0;
            }
            else //Shyam
            {
                variableResidual[0]=0.0;
            }

            //do//Iterates all the variables simulataneously till convergence
            {
                for(cntVariable=initial_numTurbScalars; cntVariable<10; cntVariable++)
                {
                    if(CgnsRoot.CgnsBase->CgnsZone->Flowsol->Scal[cntVariable]!=NULL)
                    {
                        CallCalcCoefficientsForVariables(CgnsRoot.CgnsBase,0,600+cntVariable);
                        CDcoeffCalculation(CgnsRoot.CgnsBase,600+cntVariable);

                        CalcSource_opti(CgnsRoot.CgnsBase,600+cntVariable);
                        CalcVariable(CgnsRoot.CgnsBase,600+cntVariable);

                        //variableResidual[cntVariable+1]=CalcVariableOneIter(CgnsRoot.CgnsBase,600+cntVariable);
                        //ShearStress(CgnsRoot.CgnsBase,600+cntVariable);
                        //printf("residual scalar[%d]=%e\n",cntVariable,variableResidual[cntVariable+1]);
                    }
                    else
                        variableResidual[cntVariable+1]=0.0;
                }
                //variableIterTag=FindVariableIterToken(variableResidual,CgnsRoot.CgnsBase);
                //printf("variableIterTag = %d\n",variableIterTag);
            }
            //while(variableIterTag==1);

            //Calculates density,viscosity,ScalCoeff[1->10] etc. ONLY INSIDE THE DOMAIN, if they are variables/functions.
            EvaluateVarCoefficients(CgnsRoot.CgnsBase);

            for(cntVariable=700; cntVariable<=713; cntVariable++)
                RefreshVariableCoefficients(cntVariable,CgnsRoot.CgnsBase);
            if(CgnsRoot.CgnsBase->flagMPFV==1)
            {
                for(nophase=0; nophase<nop; nophase++)
                    RefreshVariableCoefficients(720+nophase,CgnsRoot.CgnsBase);
            }
            FindTimeStepErrorArray(CgnsRoot.CgnsBase,ErrorTimeStep);

            timeStepToken = FindtimeStepToken(CgnsRoot.CgnsBase,ErrorTimeStep,CounterTstep);
            SaveTimeStepData(CgnsRoot,CounterTstep,CgnsRoot.CgnsBase->BaseParameters->Time);

            RefreshMassVelocities(CgnsRoot.CgnsBase);///This function is added to refresh the mass velocities at the boundary.  This call is mandatory in case of restarting the turbulence case.  Don't remove. Bhanesh.
            CgnsRoot.CgnsBase->masstag=0;//chandrasekhar
            for(index_var=0; index_var<total; index_var++)
            {
                vartag = CgnsRoot.CgnsBase->VARTAG[index_var];
                TimeDepend(CgnsRoot.CgnsBase->BaseParameters->Time,vartag,CgnsRoot.CgnsBase);//Alpha Beta Gamma values have to be updated here
                RefreshSurfaceCells(vartag,CgnsRoot.CgnsBase);
            }

            if(CgnsRoot.CgnsBase->flagTurb==1)
            {
                if(!(CgnsRoot.CgnsBase->flagTurbModel==4 || CgnsRoot.CgnsBase->flagTurbModel==5))
                    CalcUtau(CgnsRoot.CgnsBase);

                CalcGamma(CgnsRoot.CgnsBase);
                CalcTurbViscosity(CgnsRoot.CgnsBase);
                freeGamma(CgnsRoot.CgnsBase);

                //Added by chandrasekhar --- to caluculate diffusion coefficients for 31(Global momentum),35(Disperse phase Volume fraction),44(Disperse phase momentum)
                if(CgnsRoot.CgnsBase->flagMPFV==0)
                {
                    CalcMuEff(CgnsRoot.CgnsBase,1);
                    if(CgnsRoot.CgnsBase->flagParticle==1)
                    {
                        CalcMuEff(CgnsRoot.CgnsBase,62);//Diff coeff for disperse phase momentum in case of gas-particle flow
                        CalcMuEff(CgnsRoot.CgnsBase,65);//Diff coeff for volume fraction of disperse phase in case of gas-particle flow
                    }
                }
                else if(CgnsRoot.CgnsBase->flagMPFV==1)
                {
                    CalcMuEff(CgnsRoot.CgnsBase,31);
                    CalcMuEff(CgnsRoot.CgnsBase,35);
                    CalcMuEff(CgnsRoot.CgnsBase,44);
                }
                //End --- chandrasekhar
            }
            pmonitors(&CgnsRoot,CgnsRoot.CgnsBase);

            if(CounterTstep==T_Bicg)
            {
                if(CgnsRoot.CgnsBase->BiCgflag==1)
                {
                    printf("Changing to OGS\n");//can be changed to GS also.
                    CgnsRoot.CgnsBase->BiCgflag=2;
                    for(numzones=0; numzones<nzone; numzones++)
                    {
                        zone = CgnsRoot.CgnsBase->CgnsZone+numzones;

                        zone->Correction= (double *)initialize_variable(1,zone,0.0);
                        zone->AP= (double *)initialize_variable(1,zone,0.0);

                        zone->Correction_Prev= (double *)initialize_variable(1,zone,0.0);
                        zone->AP_Prev= (double *)initialize_variable(1,zone,0.0);
// 						free(zone->Residue);
                        free(zone->Residue_init);
                        free(zone->Pi);
                        free(zone->Vi);
                        free(zone->Si);
                        free(zone->Ti);
                    }
                }
            }

        }//End of time stepping loop
        while(timeStepToken==1);

        if(CgnsRoot.CgnsBase->flagLMMHD==1)
        {
            //CalcCurrentDensity(CgnsRoot.CgnsBase);
        }

        printf("Printing the solver in data files of \"SOLUTIONS\" folder ...\n");
        printsolver(CgnsRoot.CgnsBase);
    }
    else
    {
        // very important to put it
        if((CgnsRoot.CgnsBase->flagSaturation==1)&&(CgnsRoot.CgnsBase->flagVelocity==0))
        {
            UpdateWaterContent(CgnsRoot.CgnsBase);
            RefreshVariableCoefficients(719,CgnsRoot.CgnsBase);
        }

        volumeflux(CgnsRoot.CgnsBase);

        do
        {
            if (timetag==1 && time==0)
                CgnsRoot.CgnsBase->BaseParameters->TimeStepingSchemeTag=0;
            else if (timetag==1 && time>0)
                CgnsRoot.CgnsBase->BaseParameters->TimeStepingSchemeTag=1;

            time = time+dt;
            CounterTstep=CounterTstep+1;

// 			fp=fopen("SOLUTIONS/No_of_TimeSteps","w");
// 			fprintf(fp,"No of Time Steps: %d\nTime: %lf", CounterTstep,time);
// 			fclose(fp);

            convcri=fopen("Settings/convergenceCriteria.dat","r");
            for(cntVariable=0; cntVariable<12; cntVariable++)
                fscanf(convcri,"%lf",(CgnsRoot.CgnsBase->BaseParameters->convergenceCriterion+cntVariable));
            fclose(convcri);

            omega=fopen("Settings/omega.dat","r");
            for(cntVariable=0; cntVariable<14; cntVariable++)
                fscanf(omega,"%lf",(CgnsRoot.CgnsBase->BaseParameters->OMEGA+cntVariable));
            fclose(omega);

            printf("\nTime: %lf\n",time);

            save=fopen("Settings/iwrite.dat","r");
            fscanf(save,"%d",&stop);
            fclose(save);
            if(stop==0)
            {
		printf("Writing the solution data ...\n");
                WriteCgnsFile(filename2,CgnsRoot);
		WriteCgnsFile_vertex("SOLUTIONS/vertex.cgns",CgnsRoot);
                printsolver(CgnsRoot.CgnsBase);
                dt=CgnsRoot.CgnsBase->BaseParameters->TimeStep;
                odt=fopen("Settings/old_dt.dat","w");
                fprintf(odt,"%12.10lf %12.10lf",dt,CgnsRoot.CgnsBase->BaseParameters->Time);
                fclose(odt);
		Write_Solution(CgnsRoot.CgnsBase,"SOLUTIONS/RESTART_DATA/solution.dat");
		//Write_gnuplot_style_solution(CgnsRoot.CgnsBase,"SOLUTIONS/gnuplot/gnu_solution.dat");
                printslab(CgnsRoot.CgnsBase);
                exit(1);
            }

            FillNewsolToPrevsol(CgnsRoot.CgnsBase);

            savevolumeflux(CgnsRoot.CgnsBase);

            if(CgnsRoot.CgnsBase->flagVelocity==0)
            {
                PressureHeadIter(CgnsRoot.CgnsBase);
            }

            CalVelocity3d(CgnsRoot.CgnsBase);

            volumeflux(CgnsRoot.CgnsBase);

            //##########################################$$$$$$$$$$$$$$$$$$$$$$##############//////////

            for(cntVariable=0; cntVariable<10; cntVariable++)
            {
                if(CgnsRoot.CgnsBase->CgnsZone->Flowsol->Scal[cntVariable]!=NULL )
                {
                    CalDispersionCoeff3d(CgnsRoot.CgnsBase,cntVariable);
                    RefreshVariableCoefficients(700+cntVariable,CgnsRoot.CgnsBase);
                }
            }

            variableResidual[0] = 0.0; // This is for Temperature and we are not solving for Temperature.

            do
            {
                for(cntVariable=0; cntVariable<10; cntVariable++)
                {
                    if(CgnsRoot.CgnsBase->CgnsZone->Flowsol->Scal[cntVariable]!=NULL)
                    {
                        CalGradScalerXY(CgnsRoot.CgnsBase,600+cntVariable,cntVariable);
                        variableResidual[cntVariable+1]=CalcVariableOneIter_porousMedia(CgnsRoot.CgnsBase,600+cntVariable);//this subroutine is in file APKCalcVariableOneIter.c
                    }
                    else
                        variableResidual[cntVariable+1]=0.0;
                }

                variableIterTag=FindVariableIterToken(variableResidual,CgnsRoot.CgnsBase);
                printf("\n Residual Scaler=%e",variableResidual[0]);

            }
            while(variableIterTag==1);
            //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$////////////////////////////////
            CgnsRoot.CgnsBase->BaseParameters->Time=CgnsRoot.CgnsBase->BaseParameters->Time+dt;

            if(time >=(CgnsRoot.CgnsBase->BaseParameters->RequiredTime))
                timeStepToken=0;
            else
                timeStepToken =1;

            SaveTimeStepData(CgnsRoot,CounterTstep,CgnsRoot.CgnsBase->BaseParameters->Time);

        }
        while(timeStepToken==1);

    }

    fclose(Q);
    return 0;
}


