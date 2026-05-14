#include "APKCalcNewParticleVelocity.h"
 
int CalcNewParticleVelocity(Base *base)
{
	int i,j,k,NI,NJ,NK,Index,*ptrSize;
	double Rep,VelA[3],EfieldA[3],Efielda,ParticleVelocity[3];
	double ModRelVelocity,ParticleCharge,dt,RhoIon,MobIon;
	//double dragcoeff,chargetest;
	double Mp,dp,RhoPart,RelPermitivity,Rho,Viscoeff,Permitivity;
	double pi=3.142857;
	Zone *desiredZone;
	
	dt=((base->BaseParameters)->deltat);
//	printf("i m inside CalcNewParticleVelocity");
	Mp=((base->BaseParameters)->Mass_ofParticles);
	dp=((base->BaseParameters)->Diameter_ofParticles);
	RhoPart=((base->BaseParameters)->Density_ofParticles);
	RelPermitivity=((base->BaseParameters)->RelativePermitivity);

	VelA[0]=((base->BaseParameters)->ParticlePosition->VelA[0]);
	VelA[1]=((base->BaseParameters)->ParticlePosition->VelA[1]);
	VelA[2]=((base->BaseParameters)->ParticlePosition->VelA[2]);
	EfieldA[0]=((base->BaseParameters)->ParticlePosition->EfieldA[0]);
	EfieldA[1]=((base->BaseParameters)->ParticlePosition->EfieldA[1]);
	EfieldA[2]=((base->BaseParameters)->ParticlePosition->EfieldA[2]);
	ParticleVelocity[0]=((base->BaseParameters)->ParticlePosition->ParticleVelocity[0]);
	ParticleVelocity[1]=((base->BaseParameters)->ParticlePosition->ParticleVelocity[1]);
	ParticleVelocity[2]=((base->BaseParameters)->ParticlePosition->ParticleVelocity[2]);
	
	//ParticleCharge=((base->BaseParameters)->ParticlePosition->ParticleCharge);
	
	desiredZone=((base->BaseParameters)->ParticlePosition->zone);
	i=((base->BaseParameters)->ParticlePosition->i);
	j=((base->BaseParameters)->ParticlePosition->j);
	k=((base->BaseParameters)->ParticlePosition->k);
	
	ptrSize=((desiredZone)->Size);
	NI=*(ptrSize+0);
	NJ=*(ptrSize+1);
	NK=*(ptrSize+2);
	Index=(NJ-1)*(NI-1)*k+(NI-1)*j+i;
	
	Rho=*((desiredZone->Flowsol)->Rho);
	Viscoeff=*((desiredZone->Flowsol)->Viscoeff);
	Permitivity=8.854*pow(10,-12);
	
/***********************************************************************************/	
/*************************calculation of Drag Force********************************/
/***********************************************************************************/	
	
	ModRelVelocity=sqrt((VelA[0]-ParticleVelocity[0])*(VelA[0]-ParticleVelocity[0])+(VelA[1]-ParticleVelocity[1])*(VelA[1]-ParticleVelocity[1])+(VelA[2]-ParticleVelocity[2])*(VelA[2]-ParticleVelocity[2]));
	Rep=(Rho*ModRelVelocity*dp)/Viscoeff;
	
	if(Rep<1000)
	{
	
		/*dragcoeff=(24/Rep)*(1+0.6667*pow(Rep,0.67));
		((base->BaseParameters)->ParticlePosition->DragForce[0])=(3*Rho*dragcoeff)/(4*RhoPart*dp)*ModRelVelocity*(VelA[0]-ParticleVelocity[0]);
		((base->BaseParameters)->ParticlePosition->DragForce[1])=(3*Rho*dragcoeff)/(4*RhoPart*dp)*ModRelVelocity*(VelA[1]-ParticleVelocity[1]);
		((base->BaseParameters)->ParticlePosition->DragForce[2])=(3*Rho*dragcoeff)/(4*RhoPart*dp)*ModRelVelocity*(VelA[2]-ParticleVelocity[2]);*/
				
		((base->BaseParameters)->ParticlePosition->DragForce[0])=pow(10,-3)*(18*Viscoeff*(VelA[0]-ParticleVelocity[0])*(1+0.16667*pow(Rep,0.67)))/(RhoPart*dp*dp);
		((base->BaseParameters)->ParticlePosition->DragForce[1])=pow(10,-3)*(18*Viscoeff*(VelA[1]-ParticleVelocity[1])*(1+0.16667*pow(Rep,0.67)))/(RhoPart*dp*dp);
		((base->BaseParameters)->ParticlePosition->DragForce[2])=pow(10,-3)*(18*Viscoeff*(VelA[2]-ParticleVelocity[2])*(1+0.16667*pow(Rep,0.67)))/(RhoPart*dp*dp);
			
		
//		printf("\nDragforceX= %12.8lf\tDragforceY= %12.8lf\tDragforceZ= %12.8lf\t ",((base->BaseParameters)->ParticlePosition->DragForce[0]),((base->BaseParameters)->ParticlePosition->DragForce[1]),((base->BaseParameters)->ParticlePosition->DragForce[2]));
	}
	else if(Rep>1000)
	{
				
		((base->BaseParameters)->ParticlePosition->DragForce[0])=(0.33*Rho*ModRelVelocity*(VelA[0]-ParticleVelocity[0]))/(dp*dp*RhoPart);
		((base->BaseParameters)->ParticlePosition->DragForce[1])=(0.33*Rho*ModRelVelocity*(VelA[1]-ParticleVelocity[1]))/(dp*dp*RhoPart);
		((base->BaseParameters)->ParticlePosition->DragForce[2])=(0.33*Rho*ModRelVelocity*(VelA[2]-ParticleVelocity[2]))/(dp*dp*RhoPart);
//		printf("\nDragforceX1= %12.8lf\tDragforceY1= %12.8lf\tDragforceZ1= %12.8lf\t ",((base->BaseParameters)->ParticlePosition->DragForce[0]),((base->BaseParameters)->ParticlePosition->DragForce[1]),((base->BaseParameters)->ParticlePosition->DragForce[2]));
	}
	
/***********************************************************************************/	
/******************************calculation of Electric Force************************/
/***********************************************************************************/	
	
	MobIon=*((desiredZone->Flowsol)->MobIon);
        RhoIon=*((desiredZone->Flowsol)->RhoIon+Index);

	
	Efielda=sqrt((EfieldA[0]*EfieldA[0])+(EfieldA[1]*EfieldA[1])+(EfieldA[2]*EfieldA[2]));
	
	ParticleCharge=((base->BaseParameters)->ParticlePosition->ParticleCharge);
	((base->BaseParameters)->ParticlePosition->ElectricForce[0])=pow(10,6)*(ParticleCharge*EfieldA[0])/(pi*RhoPart*dp*dp*dp/6);
	((base->BaseParameters)->ParticlePosition->ElectricForce[1])=pow(10,6)*(ParticleCharge*EfieldA[1])/(pi*RhoPart*dp*dp*dp/6);
	((base->BaseParameters)->ParticlePosition->ElectricForce[2])=pow(10,6)*(ParticleCharge*EfieldA[2])/(pi*RhoPart*dp*dp*dp/6);
	
	
	((base->BaseParameters)->ParticlePosition->SaturationCharge)=(3*RelPermitivity*pi*Permitivity*dp*dp*Efielda)/(RelPermitivity+2);
	((base->BaseParameters)->ParticlePosition->ParticleCharge)=((base->BaseParameters)->ParticlePosition->SaturationCharge);
	
	
	
/***************************************************************************/	
/**********Calculation of U,V,W Particle Velocities****************/
/***************************************************************************/
	
	((base->BaseParameters)->ParticlePosition->ParticleVelocity[0])=((base->BaseParameters)->ParticlePosition->ParticleVelocity[0])+dt*(((base->BaseParameters)->ParticlePosition->DragForce[0])+((base->BaseParameters)->ParticlePosition->ElectricForce[0]));
	
	((base->BaseParameters)->ParticlePosition->ParticleVelocity[1])=((base->BaseParameters)->ParticlePosition->ParticleVelocity[1])+dt*(((base->BaseParameters)->ParticlePosition->DragForce[1])+((base->BaseParameters)->ParticlePosition->ElectricForce[1]));
	
	((base->BaseParameters)->ParticlePosition->ParticleVelocity[2])=((base->BaseParameters)->ParticlePosition->ParticleVelocity[2])+dt*(((base->BaseParameters)->ParticlePosition->DragForce[2])+((base->BaseParameters)->ParticlePosition->ElectricForce[2]));
	
	
	 printf("\nupnew= %12.8lf\tvpnew= %12.8lf\twpnew= %12.8lf\tsaturationCharge1 =%12.8lf\n ",((base->BaseParameters)->ParticlePosition->ParticleVelocity[0]),((base->BaseParameters)->ParticlePosition->ParticleVelocity[1]),((base->BaseParameters)->ParticlePosition->ParticleVelocity[2]),((base->BaseParameters)->ParticlePosition->SaturationCharge));
	
return 0;
}


