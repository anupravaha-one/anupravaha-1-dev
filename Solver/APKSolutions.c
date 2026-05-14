#include "APKSolutions.h"
/**
The function Solutions() does the following part of the Segregated algoritm for 3D domains:
1.This function calls the massvelocity() to calculate massvelocities.
2.Calculates volumeflux from the calculated value of massvelocities.
3.Finds out Pressure using the function PressureIter().
4.Repeats the above steps till Volume Flux converges.
5.Computes velocities by using the function CalcVariabler().
6.Returns to massvelocity calculation till the total residual of all the velocity components are less than their
	individual convergence criteria.
Note that the step 6 is performed in the case of when Implicit approach is selected.
The operations stated above are performed according to the type of solver scheme(Imlicit, SemiImplicit, PseudoCompressibility)
selected. Hence this function checks for the tag of the scheme selected and performs the operations stated above accordingly.
Scheme tags:
0-> Implicit scheme
1-> SemiImplicit Scheme.
21->Pseudo Compressibility(One T one Tau step)
22->Pseudo Compressibility(Enforcing Continuity at each timestep)
*/
int Solutions(Base *base,int counter_I)
{
	double residMassFlux;
	int scheme,nZones;
	scheme=base->SolverScheme ;
	nZones=base-> Number_of_zones;
	if(scheme==0)
	{
		do
		{
			savevolumeflux(base);
			
			Conv_Coeff_Calculation(base,0,1);
			CDcoeffCalculation(base,1);
			
			CalcSource_opti(base,100);
			massVelocity(base, 1);/*supply 1 for u, 2 for v and  3 for w and no other value*/
			
			CalcSource_opti(base,200);
			massVelocity(base, 2);/*supply 1 for u, 2 for v and  3 for w and no other value*/
			
			CalcSource_opti(base,300);
			massVelocity(base, 3);/*supply 1 for u, 2 for v and  3 for w and no other value*/
			
			volumefluxstar(base);/*Calculates new Massflux from mass velocities*/
			
			if(base->BiCgflag==0)
				PressureIter(base);//Calculate new pressure from new Massflux
			else if(base->BiCgflag==1)
				Conj_PressureIter(base);
			else if(base->BiCgflag==2)
				ModiPressureIter(base);
	
			CalcNewVolumeFlux(base);
			
			residMassFlux=FindVolumeFluxDifference(base);

			printf("residMassFlux=%12.8e\n",residMassFlux);
		}
		while(residMassFlux >base->BaseParameters->MassFluxConv);
		
		Conv_Coeff_Calculation(base,0,1);
		CDcoeffCalculation(base,1);
		
		Pressure_Grad_Calculation(base);
		
		CalcSource_opti(base,1);
		CalcVariable(base,1);/*calculate for flowsol->u*/
		
		CalcSource_opti(base,2);
		CalcVariable(base,2);/*calculate for flowsol->v*/
		
		CalcSource_opti(base,3);
		CalcVariable(base,3);/*calculate for flowsol->w*/
	}
	else if(scheme==1)
	{
		if(counter_I>1)
		{
			savevolumeflux(base);
			
			Conv_Coeff_Calculation(base,0,1);
			CDcoeffCalculation(base,1);
			
			CalcSource_opti(base,100);
			massVelocity(base, 1);/*supply 1 for u, 2 for v and  3 for w and no other value*/
			
			CalcSource_opti(base,200);
			massVelocity(base, 2);/*supply 1 for u, 2 for v and  3 for w and no other value*/
			
			CalcSource_opti(base,300);
			massVelocity(base, 3);/*supply 1 for u, 2 for v and  3 for w and no other value*/
			
			volumefluxstar(base);/*Calculates new Massflux from mass velocities*/
			
			if(base->BiCgflag==0)
				PressureIter(base);//Calculate new pressure from new Massflux
			else if(base->BiCgflag==1)
				Conj_PressureIter(base);
			else if(base->BiCgflag==2)
				ModiPressureIter(base);
	
			CalcNewVolumeFlux(base);
			
			Conv_Coeff_Calculation(base,0,1);
			CDcoeffCalculation(base,1);
			
			Pressure_Grad_Calculation(base);//For Optimizing testing purpose..  Can be removed later.
			
			CalcSource_opti(base,1);
			CalcVariable(base,1);/*calculate for flowsol->u*/
			
			CalcSource_opti(base,2);
			CalcVariable(base,2);/*calculate for flowsol->v*/
			
			CalcSource_opti(base,3);
			CalcVariable(base,3);/*calculate for flowsol->w*/	
		}
		else
		{
			do
			{
				savevolumeflux(base);
				
				Conv_Coeff_Calculation(base,0,1);
				CDcoeffCalculation(base,1);
				
				CalcSource_opti(base,100);
				massVelocity(base, 1);/*supply 1 for u, 2 for v and  3 for w and no other value*/

				CalcSource_opti(base,200);
				massVelocity(base, 2);/*supply 1 for u, 2 for v and  3 for w and no other value*/
				
				CalcSource_opti(base,300);
				massVelocity(base, 3);/*supply 1 for u, 2 for v and  3 for w and no other value*/
				
				volumefluxstar(base);/*Calculates new Massflux from mass velocities*/
				
				if(base->BiCgflag==0)
					PressureIter(base);//Calculate new pressure from new Massflux
				else if(base->BiCgflag==1)
					Conj_PressureIter(base);
				else if(base->BiCgflag==2)
					ModiPressureIter(base);
	
				CalcNewVolumeFlux(base);
				
				residMassFlux=FindVolumeFluxDifference(base);
				printf("residMassFlux=%12.8e\n",residMassFlux);
			}
			while(residMassFlux > base->BaseParameters->MassFluxConv);
			
			Conv_Coeff_Calculation(base,0,1);
			CDcoeffCalculation(base,1);
			
			Pressure_Grad_Calculation(base);
			
			CalcSource_opti(base,1);
			CalcVariable(base,1);/*calculate for flowsol->u*/
			
			CalcSource_opti(base,2);
			CalcVariable(base,2);/*calculate for flowsol->v*/
			
			CalcSource_opti(base,3);
			CalcVariable(base,3);/*calculate for flowsol->w*/	
		}
	}
// 	else if (scheme==21)
// 	{
// 		if(counter_I>1)
// 		{
// 			savevolumeflux(base);
// 			
// 			CDcoeffCalculation(base,1);
// 			
// 			massVelocity(base, 1);/*supply 1 for u, 2 for v and  3 for w and no other value*/
// 			
// 			massVelocity(base, 2);/*supply 1 for u, 2 for v and  3 for w and no other value*/
// 			
// 			massVelocity(base, 3);/*supply 1 for u, 2 for v and  3 for w and no other value*/
// 			
// 			volumefluxstar(base);/*Calculates new Massflux from mass velocities*/
// 			
// 			PressureIterArtificial(base);/*Seems O .K.*//*Calculate new pressure from new Massflux*/
// 			
// 			CalcNewVolumeFlux(base);
// 			
// 			CDcoeffCalculation(base,1);
// 			
// 			CalcVariable(base,1);/*calculate for flowsol->u*/
// 			CalcVariable(base,2);/*calculate for flowsol->v*/
// 			CalcVariable(base,3);/*calculate for flowsol->w*/	
// 		}
// 		else
// 		{
// 			do
// 			{
// 				savevolumeflux(base);
// 				
// 				CDcoeffCalculation(base,1);
// 				
// 				massVelocity(base, 1);/*supply 1 for u, 2 for v and  3 for w and no other value*/
// 				
// 				massVelocity(base, 2);/*supply 1 for u, 2 for v and  3 for w and no other value*/
// 				
// 				massVelocity(base, 3);/*supply 1 for u, 2 for v and  3 for w and no other value*/
// 				
// 				volumefluxstar(base);/*Calculates new Massflux from mass velocities*/
// 				
// 				PressureIter(base);/*Calculate new pressure from new Massflux*/
// 				
// 				CalcNewVolumeFlux(base);
// 				
// 				residMassFlux=FindVolumeFluxDifference(base);
// 			}
// 			while(residMassFlux > base->BaseParameters->MassFluxConv);
// 			
// 			CDcoeffCalculation(base,1);
// 			
// 			CalcVariable(base,1);/*calculate for flowsol->u*/
// 			CalcVariable(base,2);/*calculate for flowsol->v*/
// 			CalcVariable(base,3);/*calculate for flowsol->w*/	
// 		}
// 	}
// 	else if (scheme==22)
// 	{
// 		if(counter_I>1)
// 		{
// 			savevolumeflux(base);
// 			
// 			CDcoeffCalculation(base,1);
// 			
// 			massVelocity(base, 1);/*supply 1 for u, 2 for v and  3 for w and no other value*/
// 			
// 			massVelocity(base, 2);/*supply 1 for u, 2 for v and  3 for w and no other value*/
// 			
// 			massVelocity(base, 3);/*supply 1 for u, 2 for v and  3 for w and no other value*/
// 
// 			volumefluxstar(base);/*Calculates new Massflux from mass velocities*/
// 
// 			do
// 			{
// 				for(cntZ=0;cntZ<nZones;cntZ++)
// 				{
// 					fillnewtoprev(((base->CgnsZone+cntZ)->Flowsol)->P,((base->CgnsZone+cntZ)->PrevFlowsol)->P,(base->CgnsZone+cntZ)->Size);
// 				}
// 				PressureIterArtificial(base);/*Seems O .K.*//*Calculate new pressure from new Massflux*/
// 				maxerrorP=CalculateErrorPressure(base);
// 				//printf("maxerrorP=%lf\n",maxerrorP);
// 			}
// 			while(maxerrorP>1e-6);//*Calculate new pressure from new Massflux*/
// 
// 			CalcNewVolumeFlux(base);
// 			
// 			CDcoeffCalculation(base,1);
// 			
// 			CalcVariable(base,1);/*calculate for flowsol->u*/
// 			CalcVariable(base,2);/*calculate for flowsol->v*/
// 			CalcVariable(base,3);/*calculate for flowsol->w*/	
// 		}
// 		else
// 		{
// 			do
// 			{
// 				savevolumeflux(base);
// 				
// 				CDcoeffCalculation(base,1);
// 				
// 				massVelocity(base, 1);/*supply 1 for u, 2 for v and  3 for w and no other value*/
// 				
// 				massVelocity(base, 2);/*supply 1 for u, 2 for v and  3 for w and no other value*/
// 				
// 				massVelocity(base, 3);/*supply 1 for u, 2 for v and  3 for w and no other value*/
// 				
// 				volumefluxstar(base);/*Calculates new Massflux from mass velocities*/
// 				
// 				PressureIter(base);/*Seems O .K.*//*Calculate new pressure from new Massflux*/
// 				
// 				CalcNewVolumeFlux(base);/*seems O. K.   */
// 				
// 				residMassFlux=FindVolumeFluxDifference(base);
// 			}
// 			while(residMassFlux > base->BaseParameters->MassFluxConv);
// 			
// 			CDcoeffCalculation(base,1);
// 			
// 			CalcVariable(base,1);/*calculate for flowsol->u*/
// 			CalcVariable(base,2);/*calculate for flowsol->v*/
// 			CalcVariable(base,3);/*calculate for flowsol->w*/	
// 		}
// 	}
	return 0;
}
