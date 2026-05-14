#include "APKParticleTracking.h"

int ParticleTracking(Root CgnsRoot)
{
	double xp,yp,zp,dt,Lx,Ly;
	int Particlenum,No_Time_Step,ParticlesOut,nv;
	//int dypart;
	//char *str1,*str2,*str3;
	FILE *ESPEff,*Particle;
	
	Particle=fopen("SOLUTIONS/ParticleSolution.dat","w");
	ESPEff=fopen("SOLUTIONS/ESPEfficency.dat","w");
	
	//str1=(char *)malloc(40*sizeof (char));
	//str2=(char *)malloc(40*sizeof (char));
	//str3=(char *)malloc(40*sizeof (char));
	
	nv=((CgnsRoot.CgnsBase)->BaseParameters->No_ofParticles);
	dt=((CgnsRoot.CgnsBase)->BaseParameters->deltat);
	Lx=((CgnsRoot.CgnsBase)->BaseParameters->ESPLength);
	Ly=((CgnsRoot.CgnsBase)->BaseParameters->WirePlateSpacing);
	
	ParticlesOut=0;
	//dypart=0.013409091;
	((CgnsRoot.CgnsBase)->BaseParameters->TotalWeightage)=0.0;
	((CgnsRoot.CgnsBase)->BaseParameters->Weightage)=0.0;
	for(Particlenum=0;Particlenum<nv;Particlenum++)
	{	
		//sprintf(str1,"particle");
		//sprintf(str2,"_%d.dat",Particlenum);
		//str3 = strcat (str1,str2);
		//Particle = fopen(str3,"w");
		
		fprintf(Particle,"\n*****************Particle number:%d******************\n",Particlenum);
		xp=*(((CgnsRoot.CgnsBase->BaseParameters)->ParticlePosition)->Xp+Particlenum);
		yp=*(((CgnsRoot.CgnsBase->BaseParameters)->ParticlePosition)->Yp+Particlenum);
		zp=*(((CgnsRoot.CgnsBase->BaseParameters)->ParticlePosition)->Zp+Particlenum);
		
		GlobalParticleTracking(CgnsRoot.CgnsBase,xp,yp,zp);

		InterpolateVelocity(CgnsRoot.CgnsBase,xp,yp,zp);

		((CgnsRoot.CgnsBase)->BaseParameters->ParticlePosition->Uactual)=((CgnsRoot.CgnsBase->BaseParameters)->ParticlePosition->VelA[0]);
		
		((CgnsRoot.CgnsBase)->BaseParameters->ParticlePosition->ParticleVelocity[0])=((CgnsRoot.CgnsBase->BaseParameters)->ParticlePosition->VelA[0]);
		((CgnsRoot.CgnsBase)->BaseParameters->ParticlePosition->ParticleVelocity[1])=((CgnsRoot.CgnsBase->BaseParameters)->ParticlePosition->VelA[1]);
		((CgnsRoot.CgnsBase)->BaseParameters->ParticlePosition->ParticleVelocity[2])=((CgnsRoot.CgnsBase->BaseParameters)->ParticlePosition->VelA[2]);
		((CgnsRoot.CgnsBase)->BaseParameters->ParticlePosition->ParticleCharge)=0.0;
		
		No_Time_Step=0;
		do
		{
			xp=*(((CgnsRoot.CgnsBase->BaseParameters)->ParticlePosition)->Xp+Particlenum);
			yp=*(((CgnsRoot.CgnsBase->BaseParameters)->ParticlePosition)->Yp+Particlenum);
			zp=*(((CgnsRoot.CgnsBase->BaseParameters)->ParticlePosition)->Zp+Particlenum);
//			printf("XYZ corordinates123 = %lf\t%lf\t%lf\n",xp,yp,zp);
			QuickParticleTracking(CgnsRoot.CgnsBase,xp,yp,zp);
			
			
			InterpolateVelocity(CgnsRoot.CgnsBase,xp,yp,zp);
				
			
			InterpolateEfield(CgnsRoot.CgnsBase,xp,yp,zp);
			
			
			CalcNewParticleVelocity(CgnsRoot.CgnsBase);
			
			/*******************************************************************************/
			/***************Calculation of New Particle Position****************************/
			/*******************************************************************************/
			
			/**(((CgnsRoot.CgnsBase->BaseParameters)->ParticlePosition)->Xp+Particlenum)=*(((CgnsRoot.CgnsBase->BaseParameters)->ParticlePosition)->Xp+Particlenum)+(dt*((CgnsRoot.CgnsBase->BaseParameters)->ParticlePosition->ParticleVelocity[0]));
			*(((CgnsRoot.CgnsBase->BaseParameters)->ParticlePosition)->Yp+Particlenum)=*(((CgnsRoot.CgnsBase->BaseParameters)->ParticlePosition)->Yp+Particlenum)+(dt*((CgnsRoot.CgnsBase->BaseParameters)->ParticlePosition->ParticleVelocity[1]));
			*(((CgnsRoot.CgnsBase->BaseParameters)->ParticlePosition)->Zp+Particlenum)=*(((CgnsRoot.CgnsBase->BaseParameters)->ParticlePosition)->Zp+Particlenum)+(dt*((CgnsRoot.CgnsBase->BaseParameters)->ParticlePosition->ParticleVelocity[2]));
			
			printf("\nXparticleNew=%12.8lf\tYparticleNew=%12.8lf\tZparticleNew=%12.8lf\n",*((CgnsRoot.CgnsBase->BaseParameters)->ParticlePosition->Xp+Particlenum),*((CgnsRoot.CgnsBase->BaseParameters)->ParticlePosition->Yp+Particlenum),*((CgnsRoot.CgnsBase->BaseParameters)->ParticlePosition->Zp+Particlenum));*/
			
			No_Time_Step=No_Time_Step+1;
			
			
			
			fprintf(Particle,"\n%d\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\n",No_Time_Step,*((CgnsRoot.CgnsBase->BaseParameters)->ParticlePosition->Xp+Particlenum),*((CgnsRoot.CgnsBase->BaseParameters)->ParticlePosition->Yp+Particlenum),*((CgnsRoot.CgnsBase->BaseParameters)->ParticlePosition->Zp+Particlenum),((CgnsRoot.CgnsBase->BaseParameters)->ParticlePosition->ParticleVelocity[0]),((CgnsRoot.CgnsBase->BaseParameters)->ParticlePosition->ParticleVelocity[1]),((CgnsRoot.CgnsBase->BaseParameters)->ParticlePosition->ParticleVelocity[2]));
			
			*(((CgnsRoot.CgnsBase->BaseParameters)->ParticlePosition)->Xp+Particlenum)=*(((CgnsRoot.CgnsBase->BaseParameters)->ParticlePosition)->Xp+Particlenum)+(dt*((CgnsRoot.CgnsBase->BaseParameters)->ParticlePosition->ParticleVelocity[0]));
			*(((CgnsRoot.CgnsBase->BaseParameters)->ParticlePosition)->Yp+Particlenum)=*(((CgnsRoot.CgnsBase->BaseParameters)->ParticlePosition)->Yp+Particlenum)+(dt*((CgnsRoot.CgnsBase->BaseParameters)->ParticlePosition->ParticleVelocity[1]));
			*(((CgnsRoot.CgnsBase->BaseParameters)->ParticlePosition)->Zp+Particlenum)=*(((CgnsRoot.CgnsBase->BaseParameters)->ParticlePosition)->Zp+Particlenum)+(dt*((CgnsRoot.CgnsBase->BaseParameters)->ParticlePosition->ParticleVelocity[2]));
			
			printf("\nXparticleNew=%12.8lf\tYparticleNew=%12.8lf\tZparticleNew=%12.8lf\n",*((CgnsRoot.CgnsBase->BaseParameters)->ParticlePosition->Xp+Particlenum),*((CgnsRoot.CgnsBase->BaseParameters)->ParticlePosition->Yp+Particlenum),*((CgnsRoot.CgnsBase->BaseParameters)->ParticlePosition->Zp+Particlenum));
			
			//printf("\nlength=%lf",Lx);
			
			
			
		}
		while(*((CgnsRoot.CgnsBase->BaseParameters)->ParticlePosition->Xp+Particlenum) < Lx && *((CgnsRoot.CgnsBase->BaseParameters)->ParticlePosition->Yp+Particlenum) < Ly && *((CgnsRoot.CgnsBase->BaseParameters)->ParticlePosition->Xp+Particlenum)>0.0 && *((CgnsRoot.CgnsBase->BaseParameters)->ParticlePosition->Yp+Particlenum)>0.0);
		
		
		((CgnsRoot.CgnsBase)->BaseParameters->TotalWeightage)=((CgnsRoot.CgnsBase)->BaseParameters->TotalWeightage)+(((CgnsRoot.CgnsBase)->BaseParameters->ParticlePosition->Uactual)*0.001460);
		
		if(*((CgnsRoot.CgnsBase->BaseParameters)->ParticlePosition->Yp+Particlenum) > Ly)
		{
			((CgnsRoot.CgnsBase)->BaseParameters->Weightage)=((CgnsRoot.CgnsBase)->BaseParameters->Weightage)+(((CgnsRoot.CgnsBase)->BaseParameters->ParticlePosition->Uactual)*0.001460);
			ParticlesOut=ParticlesOut+1;
		}
		else
		ParticlesOut=ParticlesOut+0;
	}
	
	((CgnsRoot.CgnsBase)->BaseParameters->ParticlesOut)=ParticlesOut;
	
//	printf("\nNumber of Particles Collected = %d\nCollection Efficiency of ESP=%lf\n",ParticlesOut,((CgnsRoot.CgnsBase)->BaseParameters->CollectionEfficiency));
	printf("\nNumber of Particles  = %d\n",((CgnsRoot.CgnsBase)->BaseParameters->No_ofParticles));
	printf("\nNumber of Particles collected  = %d\n",((CgnsRoot.CgnsBase)->BaseParameters->ParticlesOut));
	//printf("\ndyparticle=%lf",dypart);
	//((CgnsRoot.CgnsBase)->BaseParameters->CollectionEfficiency)=((CgnsRoot.CgnsBase)->BaseParameters->ParticlesOut)/((CgnsRoot.CgnsBase)->BaseParameters->No_ofParticles);
	
//	printf("\nUactaul  = %12.8lf\n",((CgnsRoot.CgnsBase)->BaseParameters->ParticlePosition->Uactual));
	printf("\nWeitage  = %12.8lf\n",((CgnsRoot.CgnsBase)->BaseParameters->Weightage));
	printf("\nTotalWeitage  = %12.8lf\n",((CgnsRoot.CgnsBase)->BaseParameters->TotalWeightage));
	
	((CgnsRoot.CgnsBase)->BaseParameters->CollectionEfficiency)=((CgnsRoot.CgnsBase)->BaseParameters->Weightage)/((CgnsRoot.CgnsBase)->BaseParameters->TotalWeightage);
	
	printf("\nCollection Efficiency  = %12.8lf\n",((CgnsRoot.CgnsBase)->BaseParameters->CollectionEfficiency));
	//CollectionEfficiency=((CgnsRoot.CgnsBase)->BaseParameters->ParticlesOut)/((CgnsRoot.CgnsBase)->BaseParameters->No_ofParticles);
	fprintf(ESPEff,"\nNumber of Particles Collected = %d\nCollection Efficiency of ESP=%12.8lf\n",((CgnsRoot.CgnsBase)->BaseParameters->ParticlesOut),((CgnsRoot.CgnsBase)->BaseParameters->CollectionEfficiency));
	
	//printf("\nlength of ESP=%12.8lf  wire to plate spacing=%12.8lf\n",Lx,Ly);
	//printf("Required123 i=%d j=%d k=%d zonename=%s \n",(((CgnsRoot.CgnsBase->BaseParameters)->ParticlePosition)->i),(((CgnsRoot.CgnsBase->BaseParameters)->ParticlePosition)->j),(((CgnsRoot.CgnsBase->BaseParameters)->ParticlePosition)->k),(((CgnsRoot.CgnsBase->BaseParameters)->ParticlePosition)->zone->zonename));
	//fclose(Particle);
	fclose(ESPEff);
	fclose(Particle);
	return 0;
}
