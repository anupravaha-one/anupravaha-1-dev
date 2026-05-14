#include "APKsolverSegregated.h"
/**
The function solverSegregated() is the solver for two dimensional geometries which solves all the variables
in a segregated approach. This function is a control function which calls different other functions in an order
as given in the segregated algorithm. See Project Report No.1, December 2005, Section 3.10.0.2 Page 48.

Brief about Algo:
This function calls the massvelocity() to calculate massvelocities.
Calculates volumeflux from the calculated value of massvelocities.
Finds out Pressure using the function PressureIter().
Repeats the above steps till Volume Flux converges.
Computes velocities by using the function CalcVariabler().
Returns to massvelocity calculation till the total residual of all the velocity components are less than their
individual convergence criteria.
Computes all other solvable Scalars like Temperature etc seperately till they individually converge.
Saves the current timestep data if the user chooses so, or saves after a desired no of time steps as choosen by the user
using the function SaveTimeStepData()
If the user chooses to stop the computation at any point(by changing the value to 0 from 1 in the file "iwrite.dat" in the
directory "Settings"), then the solution at that time level is written using the function WriteCgnsFile().

Note: The part of calculation of Mass velocities, Pressure Iteration and Velocities computation is done by the
functions Solutions2d() or Solutions() depending on 2d or 3d domains.
*/
int solverSegregated(Root CgnsRoot)
{
	double time=0.0,dt,CON;
	double errorTimeStep[11];
	int timeStepToken,vartag;
	FILE *save,*nscon,*odt,*file,*fp,*convcri,*omega;
	char filename2[] = "SOLUTIONS/restart.cgns";
	int stop=1,iScal,con,timetag,counter=0,cntVariable,z;
	double maxerror;
	//for Solidification and Mleting problems BHANS
	int solidification_token=0,Iter_counter=0;
	double solidification_rms;
	Zone *zone;
	//Radiative Variable
	double *Kp,Irrad,*Omega,*Int,sigma=5.67e-8,*TempP,*TempN;
	int i,j,k,l,m,ni,nj,nk,nt,np,nzone,Z_ind,INDEX,T_Bicg=350;
	
	
	nt=CgnsRoot.CgnsBase->BaseParameters->ntheta;
	np=CgnsRoot.CgnsBase->BaseParameters->nphi;
	nzone=CgnsRoot.CgnsBase->Number_of_zones;
	
	file=fopen("SOLUTIONS/MassConservation.dat","w");
	nscon=fopen("Settings/converg.dat","r");
	fscanf(nscon,"%d",&con);
	fclose(nscon);
	
	FILE *countBiCg;

	countBiCg=fopen("Settings/Switch_Iter.dat","r");
	
	fscanf(countBiCg,"%d",&T_Bicg);
	
	fclose(countBiCg);
	
	dt=CgnsRoot.CgnsBase->BaseParameters->TimeStep;
	
	timetag = CgnsRoot.CgnsBase->BaseParameters->TimeStepingSchemeTag;
	//CgnsRoot.CgnsBase->BaseParameters->Time=0.0;
	
	CDcoeffMemAlloc(CgnsRoot.CgnsBase);
	
	DiffusionCoeffCalculation(CgnsRoot.CgnsBase);
	
	//DiffusionCoeffCalculationForVariables(CgnsRoot.CgnsBase);

	if(CgnsRoot.CgnsBase->flagSteadyState==0)
	{
		

// 		if(CgnsRoot.CgnsBase->CgnsZone->Flowsol->T!=NULL)
// 		{
// 			///For steady state problems giving pure diffusion equation solution as a initial conditions for solver will improve the convergence rate.  This case really helps in case of problems where CP values are extremely high and the chosen time step is low because of the flow.
// 			printf("Solving diffusion equation for Temperature .....\n\n");
//			CallCalcCoefficientsForVariables(CgnsRoot.CgnsBase,0,6);
// 			CalcConjugate(CgnsRoot.CgnsBase,6);
// 			WriteCgnsFile(filename2,CgnsRoot);
		
// 			printf("The steady state temperature distribution with other variables is written to restart.cgns\n\n");
// 		}	

		if(CgnsRoot.CgnsBase->BaseParameters->flag_Potentialflow==1)
		{
			///This is the solver for potential flow solution.  Giving potential flow solution as a initial conditions for the solver improves the cotinuity guess and pressure convergence.  If the user chooses to solve for the potential flow then the following functions will be called.
			printf("Solving for potential flow.....\n\n");
			if(CgnsRoot.CgnsBase->BaseParameters->flag2d == 1)
			{
				if(CgnsRoot.CgnsBase->BiCgflag==0)
					CalcConjugate(CgnsRoot.CgnsBase,990);
				else if(CgnsRoot.CgnsBase->BiCgflag==1)
					Conj_CalcConjugate(CgnsRoot.CgnsBase,990);
				else if(CgnsRoot.CgnsBase->BiCgflag==2)
					ModiCalcConjugate(CgnsRoot.CgnsBase,990);
			}
			else
				CalcConjugate(CgnsRoot.CgnsBase,990);
			
			find_velocities(CgnsRoot.CgnsBase);
			RefreshSurfaceCells(1,CgnsRoot.CgnsBase);
			RefreshSurfaceCells(2,CgnsRoot.CgnsBase);
			
			if(CgnsRoot.CgnsBase->BaseParameters->flag2d==0)
			RefreshSurfaceCells(3,CgnsRoot.CgnsBase);
	
			if(CgnsRoot.CgnsBase->BaseParameters->flag2d==0)
				volumeflux(CgnsRoot.CgnsBase);
			else
			volumeflux2d(CgnsRoot.CgnsBase);
			
			//printf("The mass flux difference is =%e\n",MassFluxConv(CgnsRoot.CgnsBase));
			WriteCgnsFile("SOLUTIONS/potentialflow.cgns",CgnsRoot);
			//return(1);
			
			printf("The Potential flow solution with other variables is written to potentialflow.cgns\n\n");
		}
	}
	if(con==0)
	{
		do
		{
			if (timetag==1 && time==0) /*Even the scheme is CN, the first timestep is done by BTCS scheme*/
				CgnsRoot.CgnsBase->BaseParameters->TimeStepingSchemeTag=0;
			else if (timetag==1 && time>0)
				CgnsRoot.CgnsBase->BaseParameters->TimeStepingSchemeTag=1;
			
			time = time+dt;
			counter=counter+1;
			
			//EvaluateVarCoefficients(CgnsRoot.CgnsBase);
#ifdef FILE_OPEN
{
// 			fp=fopen("SOLUTIONS/No_of_TimeSteps","w");
// 			fprintf(fp,"No of Time Steps: %d\nTime: %lf",counter,time);
// 			fclose(fp);
			
			
			
			save=fopen("Settings/iwrite.dat","r");
			fscanf(save,"%d",&stop);
			fclose(save);

			convcri=fopen("Settings/convergenceCriteria.dat","r");
			for(cntVariable=0;cntVariable<12;cntVariable++)
			fscanf(convcri,"%lf",(CgnsRoot.CgnsBase->BaseParameters->convergenceCriterion+cntVariable));
			fclose(convcri);
			
			omega=fopen("Settings/omega.dat","r");
			for(cntVariable=0;cntVariable<14;cntVariable++)
			fscanf(omega,"%lf",(CgnsRoot.CgnsBase->BaseParameters->OMEGA+cntVariable));
			fclose(omega);
}
#endif
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
				exit(1);
			}
			CgnsRoot.CgnsBase->BaseParameters->Time=CgnsRoot.CgnsBase->BaseParameters->Time+dt;
			
			printf("\nNo of Time Steps: %d\tTime: %lf\n", counter,CgnsRoot.CgnsBase->BaseParameters->Time);//<----This has been shifted to sovlerCoupled.c
			
			FillNewsolToPrevsol(CgnsRoot.CgnsBase);
			
			CallCalcCoefficientsForVariables(CgnsRoot.CgnsBase,0,1);
			
			if(CgnsRoot.CgnsBase->BaseParameters->flag2d==1)
			{
				//make_w_zero(CgnsRoot.CgnsBase);
				Solutions2d(CgnsRoot.CgnsBase,counter);
			}
			else
			{
				Solutions(CgnsRoot.CgnsBase,counter);
			}
			
			maxerror = CalculateError(CgnsRoot.CgnsBase);
			
			//printf("maximum error in velocity time steps %12.8e\n",maxerror);
			//printf("Time=%lf\n",CgnsRoot.CgnsBase->BaseParameters->Time);

			
			if(CgnsRoot.CgnsBase->BaseParameters->TstepChoice==0)
			{
				CgnsRoot.CgnsBase->BaseParameters->TimeStep=CalculateTimeStep(CgnsRoot.CgnsBase);
				dt=CgnsRoot.CgnsBase->BaseParameters->TimeStep;
				printf("deltaT=%lf\n",CgnsRoot.CgnsBase->BaseParameters->TimeStep);
			}
			
			SaveTimeStepData(CgnsRoot,counter,CgnsRoot.CgnsBase->BaseParameters->Time);
			
			for(vartag=1;vartag<5;vartag++)
			{
				if(CgnsRoot.CgnsBase->BaseParameters->flag2d==1 && vartag!=3)
				{
					TimeDepend(CgnsRoot.CgnsBase->BaseParameters->Time,vartag,CgnsRoot.CgnsBase);/*Alpha Beta Gamma values have to be updated here*/
					RefreshSurfaceCells(vartag,CgnsRoot.CgnsBase);
				}
			}
			print_monitors(CgnsRoot.CgnsBase);

			if(counter==T_Bicg)
			{
				if(CgnsRoot.CgnsBase->BiCgflag==1)
				{
					printf("Changing to OGS\n");
					CgnsRoot.CgnsBase->BiCgflag=2;
					for(z=0;z<nzone;z++)
					{
						zone = CgnsRoot.CgnsBase->CgnsZone+z;	
						
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
			
			
		}
		while(maxerror > *(CgnsRoot.CgnsBase->BaseParameters->convergenceCriterion));
#ifdef FILE_OPEN
{

		nscon=fopen("Settings/converg.dat","w");
		fprintf(nscon,"1");
		fclose(nscon);
}
#endif
	}
	CON = MassFluxConv(CgnsRoot.CgnsBase);
	fprintf(file,"Massflux Conservation at all Boundary =%12.8e\n",CON);
	
	printf("######### VELOCITIES HAVE CONVERGED ############\n");
	//printf("I have started to check\n");
	WriteCgnsFile(filename2,CgnsRoot);
	//Write_gnuplot_style_solution(CgnsRoot.CgnsBase,"SOLUTIONS/gnuplot/gnu_solution.dat");
	printsolver(CgnsRoot.CgnsBase);
	dt=CgnsRoot.CgnsBase->BaseParameters->TimeStep;
#ifdef FILE_OPEN
{

	odt=fopen("Settings/old_dt.dat","w");
	fprintf(odt,"%12.10lf",dt);
	fclose(odt);
}
#endif
	time=0.0;
	
	Conv_Coeff_Calculation(CgnsRoot.CgnsBase,0,1);

	counter=0;
	do// loop for Solving Temperature and Species
	{
		EvaluateVarCoefficients(CgnsRoot.CgnsBase);
		
		for(cntVariable=700;cntVariable<=713;cntVariable++)
			RefreshVariableCoefficients(cntVariable,CgnsRoot.CgnsBase);
		
		//DiffusionCoeffCalculationForVariables(CgnsRoot.CgnsBase);
		
		
		
#ifdef FILE_OPEN
{
	
// 		fp=fopen("SOLUTIONS/No_of_TimeSteps","w");
// 		fprintf(fp,"No of Time Steps: %d\nTime: %lf",counter,time);
// 		fclose(fp);
		
		time = time+dt;
		
		convcri=fopen("Settings/convergenceCriteria.dat","r");
		for(cntVariable=0;cntVariable<12;cntVariable++)
			fscanf(convcri,"%lf",(CgnsRoot.CgnsBase->BaseParameters->convergenceCriterion+cntVariable));
		fclose(convcri);
			
		omega=fopen("Settings/omega.dat","r");
		for(cntVariable=0;cntVariable<14;cntVariable++)
			fscanf(omega,"%lf",(CgnsRoot.CgnsBase->BaseParameters->OMEGA+cntVariable));
		fclose(omega);
		
		save=fopen("Settings/iwrite.dat","r");
		fscanf(save,"%d",&stop);
		fclose(save);
}
#endif		
		if(stop==0)
		{
			WriteCgnsFile(filename2,CgnsRoot);
			printsolver(CgnsRoot.CgnsBase);
			
			dt=CgnsRoot.CgnsBase->BaseParameters->TimeStep;
			odt=fopen("Settings/old_dt.dat","w");
			fprintf(odt,"%12.10lf %12.10lf",dt,CgnsRoot.CgnsBase->BaseParameters->Time);
			fclose(odt);
			Write_Solution(CgnsRoot.CgnsBase,"SOLUTIONS/RESTART_DATA/solution.dat");                        
			//Write_gnuplot_style_solution(CgnsRoot.CgnsBase,"SOLUTIONS/gnuplot/gnu_solution.dat");
                        printf("Stopping AnuPravaha as per user request\n");
			exit(1);
		}
		
		CgnsRoot.CgnsBase->BaseParameters->Time=CgnsRoot.CgnsBase->BaseParameters->Time+dt;
		printf("\nNo of Time Steps: %d\tTime: %lf\n", counter,CgnsRoot.CgnsBase->BaseParameters->Time);
		
		FillNewsolToPrevsol(CgnsRoot.CgnsBase);
		
		if(CgnsRoot.CgnsBase->CgnsZone->Flowsol->T!=NULL)
		{
			EvaluateEnergysource(CgnsRoot.CgnsBase);
			
			CallCalcCoefficientsForVariables(CgnsRoot.CgnsBase,0,6);
			CDcoeffCalculation(CgnsRoot.CgnsBase,6);
			CalcSource_opti(CgnsRoot.CgnsBase,6);
			
			if(CgnsRoot.CgnsBase->flagConjugateProblem==0 && CgnsRoot.CgnsBase->flagSteadyState==0)//Conjugate Problem and Steady State
                        	CalcVariableConjugate(CgnsRoot.CgnsBase,6);

			else 
			{
				
				if(CgnsRoot.CgnsBase->BaseParameters->flagflowRadiation==1)
				{

					if(CgnsRoot.CgnsBase->BaseParameters->flag2d==1)
					RadiationSolverWithFlow2D(CgnsRoot.CgnsBase);
					else
					RadiationSolverWithFlow3D(CgnsRoot.CgnsBase);
				
					for(Z_ind=0;Z_ind<nzone;Z_ind++)
					{
						Int=(CgnsRoot.CgnsBase->CgnsZone+Z_ind)->Flowsol->RadIntensity;
						Omega=(CgnsRoot.CgnsBase->CgnsZone+Z_ind)->ZGeom->Omega;
						ni=*((CgnsRoot.CgnsBase->CgnsZone+Z_ind)->Size+0);
						nj=*((CgnsRoot.CgnsBase->CgnsZone+Z_ind)->Size+1);
						nk=*((CgnsRoot.CgnsBase->CgnsZone+Z_ind)->Size+2);
						
						TempP=(CgnsRoot.CgnsBase->CgnsZone+Z_ind)->PrevFlowsol->T;
						TempN=(CgnsRoot.CgnsBase->CgnsZone+Z_ind)->Flowsol->T;

						Kp=(CgnsRoot.CgnsBase->CgnsZone+Z_ind)->Flowsol->kappa;
						for(k=0;k<nk-1;k++)
						for(j=0;j<nj-1;j++)
						for(i=0;i<ni-1;i++)
						{
							INDEX=(ni-1)*(nj-1)*k+(ni-1)*j+i;
							Irrad=0.0;
							for(l=0;l<nt;l++)
							for(m=0;m<np;m++)
							Irrad=Irrad+*(Int+((ni-1)*(nj-1)*k+(ni-1)*j+i)*nt*np+np*l+m)*(*(Omega+np*l+m));
							*((CgnsRoot.CgnsBase->CgnsZone+Z_ind)->Flowsol->RadEnergysource+INDEX)=*Kp*(4.0*sigma*pow(*(TempP+INDEX),3)*(*(TempN+INDEX))-Irrad);
							//printf("%e\t",*((CgnsRoot.CgnsBase->CgnsZone+Z_ind)->Flowsol->RadEnergysource+INDEX));
						}
					}
				}
				
				CalcVariable(CgnsRoot.CgnsBase,6);
			}
			
			errorTimeStep[0] = CalculateErrorTemp(CgnsRoot.CgnsBase,6);
			
			printf("Temperature residual= %e\n",errorTimeStep[0]);
			
//			printf("maximum error of Temperature %12.8e\n",errorTimeStep[0]);
			TimeDepend(CgnsRoot.CgnsBase->BaseParameters->Time,6,CgnsRoot.CgnsBase);/*Alpha Beta Gamma values have to be updated here*/
			RefreshSurfaceCells(6,CgnsRoot.CgnsBase);
			//printsolver(CgnsRoot.CgnsBase);
			//exit(0);
		}//
		else if(CgnsRoot.CgnsBase->flagSolidification==1)
		{
			EvaluateEnergysource(CgnsRoot.CgnsBase);
			
			
			errorTimeStep[0]=0.0;
			solidification_token=0;
			Iter_counter=0;
			
			CallCalcCoefficientsForVariables(CgnsRoot.CgnsBase,0,7);
			CDcoeffCalculation(CgnsRoot.CgnsBase,7);
			CalcSource_opti(CgnsRoot.CgnsBase,7);
			do
			{
				Iter_counter++;
				
				CalcVariable(CgnsRoot.CgnsBase,7);
				solidification_rms = CalculateErrorTemp(CgnsRoot.CgnsBase,7);
				if (fabs(solidification_rms-errorTimeStep[0]) <=1e-05)
				solidification_token=1;
				else
				errorTimeStep[0] = solidification_rms;				
				printf("solidification rms=%lf\n",solidification_rms);
			}while(solidification_token==0);
			
			printf("maximum error of enthalpy %12.8e\n",solidification_rms);
			TimeDepend(CgnsRoot.CgnsBase->BaseParameters->Time,7,CgnsRoot.CgnsBase);/*Alpha Beta Gamma values have to be updated here*/
			RefreshSurfaceCells(7,CgnsRoot.CgnsBase);
			
		}
		else 
			errorTimeStep[0]=0.0;	
		
		for(iScal=0;iScal<10;iScal++)
		{
			if(CgnsRoot.CgnsBase->CgnsZone->Flowsol->Scal[iScal]!=NULL)
			{
				
				CallCalcCoefficientsForVariables(CgnsRoot.CgnsBase,0,600+iScal);
 				CDcoeffCalculation(CgnsRoot.CgnsBase,600+iScal);
				CalcSource_opti(CgnsRoot.CgnsBase,600+iScal);
				CalcVariable(CgnsRoot.CgnsBase,600+iScal);
				
				errorTimeStep[iScal+1] = CalculateErrorTemp(CgnsRoot.CgnsBase,600+iScal);
				//printf("maximum error in Scalar[%d] is %12.8e\n",iScal,errorTimeStep[iScal+1]);
				printf("Scalar[%d] residual= %e\n",iScal,errorTimeStep[iScal+1]);
				
				TimeDepend(CgnsRoot.CgnsBase->BaseParameters->Time,600+iScal,CgnsRoot.CgnsBase);/*Alpha Beta Gamma values have to be updated here*/
				RefreshSurfaceCells(600+iScal,CgnsRoot.CgnsBase);
			}
			else
				errorTimeStep[iScal+1] =0.0;
		}
		
		timeStepToken=FindTimeStepTokenSegr(CgnsRoot.CgnsBase,errorTimeStep);
		
		SaveTimeStepData(CgnsRoot,counter,CgnsRoot.CgnsBase->BaseParameters->Time);
		counter++;
		print_monitors(CgnsRoot.CgnsBase);
		
	}
	while(timeStepToken==1);
	
	/*if(CgnsRoot.CgnsBase->CgnsZone->Flowsol->T!=NULL)
	{
		Var = VarMassFluxConv(CgnsRoot.CgnsBase,6);
		ShearStress(CgnsRoot.CgnsBase,6);
		fprintf(file,"Energy fluxConv at all Boundary =%12.8e\n",Var);
		printf("differance in two time steps =%12.8e\n",errorTimeStep[0]);
	}
	for(iScal=0;iScal<10;iScal++)
	{
		if(CgnsRoot.CgnsBase->CgnsZone->Flowsol->Scal[iScal]!=NULL)
		{
			Iscal = VarMassFluxConv(CgnsRoot.CgnsBase,600+iScal);
			fprintf(file,"Scalar[%d] fluxConv at all Boundary =%12.8e\n",iScal,Iscal);
			ShearStress(CgnsRoot.CgnsBase,600+iScal);
		}
	}*/
	fclose(file);
	return 0;
}
//#################################################################################################
int FindTimeStepTokenSegr(Base *base,double errorTimeStep[11])
{
	int cntVar,token=0;
	for(cntVar=0;cntVar<11;cntVar++)
		if(errorTimeStep[cntVar] > *(base->BaseParameters->convergenceCriterion+cntVar+1))
			token=1;
	return(token);
}
//#################################################################################################

double CalculateError(Base *base)
{
	int numZones,cntZ,*ptrSize,NI,NJ,NK,index;
	int i,j,k,totalcells=0;
	double rmsu,rmsv,rmsw,maxerror;
	double erroru=0.0,errorv=0.0,errorw = 0.0;
	double maxu,maxv,maxw,temp;
	Zone *zone;
	
	numZones= base->Number_of_zones;
	
	maxu=0.0;
	maxv=0.0;
	maxw=0.0;
		
	for(cntZ=0;cntZ<numZones;cntZ++)
	{
		zone = (base->CgnsZone+cntZ);
		ptrSize = zone->Size;
		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);
		totalcells = totalcells+(NI-1)*(NJ-1)*(NK-1);
		index=0;
		
		for(k=0;k<NK-1;k++)   /*Storing current solution in Previous solution*/
		for(j=0;j<NJ-1;j++)
		for(i=0;i<NI-1;i++)
		{
			temp=(*(zone->Flowsol->u+index) - *(zone->PrevFlowsol->u+index));
			erroru = erroru + temp*temp;
			if(fabs(temp)>maxu)
				maxu=fabs(temp);
			
			
			temp=(*(zone->Flowsol->v+index) - *(zone->PrevFlowsol->v+index));
			errorv = errorv + temp*temp;
			if(fabs(temp)>maxv)
				maxv=fabs(temp);
			
			if(zone->Flowsol->w!=NULL)
			{
				temp=(*(zone->Flowsol->w+index) - *(zone->PrevFlowsol->w+index));
				errorw = errorw + temp*temp;
				
				if(fabs(temp)>maxw)
				maxw=fabs(temp);

			}
			index++;
		}
	}
	
	
	//printf("Max. change in u,v,w velocities are\t%e\t%e\t%e\n",maxu,maxv,maxw);
	
	rmsu = sqrt(erroru/totalcells);
	rmsv = sqrt(errorv/totalcells);
	rmsw = sqrt(errorw/totalcells);
	
	printf("U velocity residual= %e\nV velocity residual= %e\nW velocity residual= %e\n",rmsu,rmsv,rmsw);
	
	
	if(rmsu > rmsv)
	{
		if(rmsu > rmsw)
		 maxerror = rmsu;
		else
		 maxerror = rmsw; 
	}
	else
	{
		if(rmsv > rmsw)
		 maxerror = rmsv;
		else
		 maxerror = rmsw;
	}
	return(maxerror);
}


double CalculateErrorTemp(Base *base,int vartag)
{
	int numZones,cntZ,*ptrSize,NI,NJ,NK,index;
	int i,j,k,totalcells=0;
	double rmsT,maxerror;
	double errorT=0.0,delta_t;
	Zone *zone;
	double *domain,*Prevdomain;
	numZones= base->Number_of_zones;
	if(vartag==1||vartag==2||vartag==3)
		totalcells = base->TotalCellFluid;
	else
		totalcells = base->TotalCellDomain;
		
	if(vartag==600 || vartag==601)
		delta_t = base->BaseParameters->Turb->DualTimeStep;
	else
		delta_t = base->BaseParameters->TimeStep;
	
	//printf("delta_t from solver Seg=%lf\n",delta_t);
	for(cntZ=0;cntZ<numZones;cntZ++)
	{
		zone = (base->CgnsZone+cntZ);
		ptrSize = zone->Size;
		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);
	
		domain= MakeDomainArray(zone,vartag);
		Prevdomain= MakeDomainArray(zone,vartag*10);
		
		for(k=0;k<NK-1;k++)   /*Storing current solution in Previous solution*/
		for(j=0;j<NJ-1;j++)
		for(i=0;i<NI-1;i++)
		{
			index = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
			//errorT = errorT + ((*(domain+index) - *(Prevdomain+index))*(*(domain+index) - *(Prevdomain+index)))/(delta_t*delta_t);
			errorT = errorT + ((*(domain+index) - *(Prevdomain+index))*(*(domain+index) - *(Prevdomain+index)));
		}
	}
	rmsT = sqrt(errorT/totalcells);
	//rmsT = rmsT/delta_t;
	maxerror=rmsT;
	//printf("maxerror=%12.8e\n",maxerror);
	return(maxerror);
}

int make_w_zero(Base* base)
{
	int numZones,cntZ,*ptrSize,NI,NJ,NK,index;
	int i,j,k;
	
	Zone *zone;
	
	numZones= base->Number_of_zones;
	for(cntZ=0;cntZ<numZones;cntZ++)
	{
		zone = (base->CgnsZone+cntZ);
		ptrSize = zone->Size;
		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);
		
		for(k=0;k<NK-1;k++)
		for(j=0;j<NJ-1;j++)
		for(i=0;i<NI-1;i++)
		{
			index = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
			*(zone->Flowsol->w+index)=0.0;
		}
	}
	return(0);
}


#include "APKCalcPressureGrad.h"
/**This functions is used to calculate the initial velocities from the potential solution.  The potential equation is discretized using finite volume formulation and solved over the domain.*/
int find_velocities(Base *base)
{
	int numZones,cntZ,*ptrSize,NI,NJ,NK,index;
	int i,j,k,ind[3];
	double *slabPHI[6],*PHI,*slabVolume[6],*ptrVolume,*SurfaceArea;
	Zone *zone;
	
	numZones= base->Number_of_zones;
	for(cntZ=0;cntZ<numZones;cntZ++)
	{
		zone = (base->CgnsZone+cntZ);
		ptrSize = zone->Size;
		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);
		
		MakeSlabPtrArray(zone,slabPHI,990);
		MakeSlabPtrArray(zone,slabVolume,5);
		PHI = zone->Flowsol->Flow_Potential;
		ptrVolume= zone->ZGeom->Volume;
		SurfaceArea=zone->ZGeom->SurfaceArea;
		for(k=0;k<NK-1;k++)
		for(j=0;j<NJ-1;j++)
		for(i=0;i<NI-1;i++)
		{
			index = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
			ind[0]=k;ind[1]=j;ind[2]=i;
			*(zone->Flowsol->u+index)=1.0/(*(ptrVolume+index))*CalcPressureGrad(PHI,slabPHI,ptrVolume,slabVolume,SurfaceArea,ind,ptrSize,1);
			*(zone->Flowsol->v+index)=1.0/(*(ptrVolume+index))*CalcPressureGrad(PHI,slabPHI,ptrVolume,slabVolume,SurfaceArea,ind,ptrSize,2);
			if(zone->Flowsol->w!=NULL)
			*(zone->Flowsol->w+index)=1.0/(*(ptrVolume+index))*CalcPressureGrad(PHI,slabPHI,ptrVolume,slabVolume,SurfaceArea,ind,ptrSize,3);
		}
	}
	return 0;
}



