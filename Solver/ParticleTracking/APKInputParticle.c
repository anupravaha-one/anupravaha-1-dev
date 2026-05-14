#include "APKInputParticle.h"
 
int InputParticle(Root CgnsRoot)
{
        FILE *var2,*var3;
        int nv,Particlenum,flagpart=0;
	double xp,yp,zp,dt,Lx,Ly,RhoPart,dp,RelPermitivity;
	//double mp;
        var2=fopen("Settings/coefficients.dat","a");	
        var3=fopen("Settings/ProgramReadInputValues.dat","a");

        do
	{
		printf("Particle Tracking is only implemented for ESP (Flows afftected by electric fields).  This case can also be extended to other variables easily\n\nRight now user can only use it for ESP problems");
		printf("\nEnter whether you want Particle Tracking[Enter->1] or No Particle Tracking[Enter->0]");
		scanf("%d",&flagpart);
	}while(!(flagpart == 1 || flagpart == 0));        
        fprintf(var2,"%d\n",flagpart);        
        fprintf(var3,"%d \tEnter whether you want Particle Tracking[Enter->1] or No Particle Tracking[Enter->0] \n\n",flagpart);
	CgnsRoot.CgnsBase->flagParticle_tracking=flagpart;
	
	if(CgnsRoot.CgnsBase->flagParticle_tracking==1)
	{
		printf("Enter the number of particles to be tracked:\n");
		scanf("%d",&nv);        
                fprintf(var2,"%d\n",nv);        
                fprintf(var3,"%d \tEnter the number of particles to be tracked: \n\n",nv);
		(CgnsRoot.CgnsBase)->BaseParameters->No_ofParticles = nv;
	//	printf("no of particles = %d\n",(CgnsRoot.CgnsBase->BaseParameters->No_ofParticles));
		(CgnsRoot.CgnsBase)->BaseParameters->ParticlePosition = (Position *)malloc(sizeof(Position));
		(CgnsRoot.CgnsBase)->BaseParameters->ParticlePosition->Xp=(double *)malloc(nv*sizeof(double));
		(CgnsRoot.CgnsBase)->BaseParameters->ParticlePosition->Yp=(double *)malloc(nv*sizeof(double));
		(CgnsRoot.CgnsBase)->BaseParameters->ParticlePosition->Zp=(double *)malloc(nv*sizeof(double));
		for(Particlenum=0;Particlenum<nv;Particlenum++)
		{
			printf("\nEnter initial X,Y,Z position of the Particles \n");
			scanf("%lf %lf %lf",&xp,&yp,&zp);        
                        fprintf(var2,"%lf %lf %lf\n",xp,yp,zp);        
                        fprintf(var3,"%lf %lf %lf \tEnter initial X,Y,Z position of the Particles \n\n",xp,yp,zp);
			*(((CgnsRoot.CgnsBase->BaseParameters)->ParticlePosition)->Xp+Particlenum)=xp;
			*(((CgnsRoot.CgnsBase->BaseParameters)->ParticlePosition)->Yp+Particlenum)=yp;
			*(((CgnsRoot.CgnsBase->BaseParameters)->ParticlePosition)->Zp+Particlenum)=zp;
			
		}
		
		/*for(Particlenum=0;Particlenum<nv;Particlenum++)
		{
			printf("%lf\t%lf\t%lf\n",*(((CgnsRoot.CgnsBase->BaseParameters)->ParticlePosition)->Xp+Particlenum),*(((CgnsRoot.CgnsBase->BaseParameters)->ParticlePosition)->Yp+Particlenum),*(((CgnsRoot.CgnsBase->BaseParameters)->ParticlePosition)->Zp+Particlenum));
			
		}*/
	
		printf("\nEnter the time step for particle tracking \n");
		scanf("%lf",&dt);        
                fprintf(var2,"%lf\n",dt);        
                fprintf(var3,"%lf \tEnter the time step for particle tracking \n\n",dt);
		((CgnsRoot.CgnsBase)->BaseParameters->deltat)=dt;
		printf("\nEnter the length of an Electrostatic Precipitator \n");
		scanf("%lf",&Lx);        
                fprintf(var2,"%lf\n",Lx);        
                fprintf(var3,"%lf \tEnter the length of an Electrostatic Precipitator \n\n",Lx);
		((CgnsRoot.CgnsBase)->BaseParameters->ESPLength)=Lx;
		printf("\nEnter the spacing between wire and collector plate of an Electrostatic Precipitator \n");
		scanf("%lf",&Ly);        
                fprintf(var2,"%lf\n",Ly);        
                fprintf(var3,"%lf \tEnter the spacing between wire and collector plate of an Electrostatic Precipitator \n\n",Ly);
		((CgnsRoot.CgnsBase)->BaseParameters->WirePlateSpacing)=Ly;
		
		//printf("length of ESP= %lf",((CgnsRoot.CgnsBase)->BaseParameters->ESPLength));
		
		printf("enter the diameter of the particle\n");
		scanf("%lf",&dp);        
                fprintf(var2,"%lf\n",dp);        
                fprintf(var3,"%lf \tenter the diameter of the particle \n\n",dp);
		((CgnsRoot.CgnsBase)->BaseParameters->Diameter_ofParticles)=dp;
		//printf("\nEnter the mass of the particle \n");
		//scanf("%lf",&mp);        
                //fprintf(var2,"%lf\n",mp);        
                //fprintf(var3,"%lf \tEnter the mass of the particle \n\n",mp);
		//((CgnsRoot.CgnsBase)->BaseParameters->Mass_ofParticles)=mp;
		printf("\nEnter the density of the particle \n");
		scanf("%lf",&RhoPart);        
                fprintf(var2,"%lf\n",RhoPart);        
                fprintf(var3,"%lf \tEnter the density of the particle \n\n",RhoPart);
		((CgnsRoot.CgnsBase)->BaseParameters->Density_ofParticles)=RhoPart;
		printf("\nEnter relative permitivity \n");
		scanf("%lf",&RelPermitivity);        
                fprintf(var2,"%lf\n",RelPermitivity);        
                fprintf(var3,"%lf \tEnter relative permitivity \n\n",RelPermitivity);
		((CgnsRoot.CgnsBase)->BaseParameters->RelativePermitivity)=RelPermitivity;
		//printf("\nreached");
	}
	fclose(var2);
        fclose(var3);
        return (0);
}
