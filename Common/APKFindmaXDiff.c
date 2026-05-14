#include "APKFindmaXDiff.h"
/**
This function, FindmaxDiff() is used to find maximum value of (coefficient of viscosity)/(deltaX)^2,(coefficient of viscosity)/(deltaY)^2, 
(coefficient of viscosity)/(deltaZ)^2. These are used to find deltaT(time) by using CFL condition in the case of variable viscosity.
Since the user can specify the zone in any orientation, computing with separate deltaX, deltaY & deltaZ can lead to problem.
So, combining all three deltaX,deltaY, deltaZ into single term Del
*/

int FindmaXDiff(Zone *zone,int *size,double *maxVarNUX,double *maxVarNUY,double *maxVarNUZ,int flagdim)
{
	int Ni,Nj,Nk,i,j,k,flag;
	double comparX,CordX,comparZ,comparY,CordY,CordZ;
	double *x,*y,*z,*varNU,Nu;
        double Delx,Dely,Delz,Del;
	
	Ni=*(size+0);
	Nj=*(size+1);
	Nk=*(size+2);
        
	x=zone->CgnsGrid->X;
	y=zone->CgnsGrid->Y;
	z=zone->CgnsGrid->Z;
        
	varNU=zone->Flowsol->Viscoeff;
	flag=zone->Flowsol->flagViscoeff;
	
	//printf("FLAG=%d",flag);	
	if(flagdim==0)
	{
		*maxVarNUX=0.0;
		*maxVarNUY=0.0;
		*maxVarNUZ=0.0;
	
		for(k=0;k<Nk-1;k++)
		for(j=0;j<Nj-1;j++)
		for(i=0;i<Ni-1;i++)
		{
                        Delx = (*(x+(Ni*Nj*k)+(Ni*j)+i+1)) - (*(x+(Ni*Nj*k)+(Ni*j)+i)) ;
                        Dely = (*(y+(Ni*Nj*k)+(Ni*j)+i+1)) - (*(y+(Ni*Nj*k)+(Ni*j)+i)) ;
                        Delz = (*(z+(Ni*Nj*k)+(Ni*j)+i+1)) - (*(z+(Ni*Nj*k)+(Ni*j)+i)) ;
                        Del  = sqrt(Delx*Delx+Dely*Dely+Delz*Delz);
                   
			Nu=*(varNU+((Nj-1)*(Ni-1)*k+(Ni-1)*j+i)*flag);
			CordX=(*(x+(Nj)*(Ni)*k+(Ni)*j+i+1)-*(x+(Nj)*(Ni)*k+(Ni)*j+i));

			comparX=Nu/(Del*Del);

			if(*maxVarNUX<comparX)
			*maxVarNUX=comparX;

			comparY=Nu/(Del*Del);

			if(*maxVarNUY<comparY)
			*maxVarNUY=comparY;

			comparZ=Nu/(Del*Del);

			if(*maxVarNUZ<comparZ)
			*maxVarNUZ=comparZ;

		}
	}
	else
	{
		for(k=0;k<Nk-1;k++)
		for(j=0;j<Nj-1;j++)
		for(i=0;i<Ni-1;i++)
		{
                    Delx = (*(x+(Ni*Nj*k)+(Ni*j)+i+1)) - (*(x+(Ni*Nj*k)+(Ni*j)+i)) ;
                    Dely = (*(y+(Ni*Nj*k)+(Ni*j)+i+1)) - (*(y+(Ni*Nj*k)+(Ni*j)+i)) ;
                    Del  = sqrt(Delx*Delx+Dely*Dely);
                   
                    *maxVarNUX= *(varNU+0)/Del;
                    *maxVarNUY= *(varNU+0)/Del;

                    Nu=*(varNU+((Nj-1)*(Ni-1)*k+(Ni-1)*j+i)*flag);
                    comparX=Nu/(Del*Del);
                        
                    if(*maxVarNUX<comparX)
                    *maxVarNUX=comparX;

                    comparY=Nu/(Del*Del);
                        
                    if(*maxVarNUY<comparY)
                    *maxVarNUY=comparY;
                    
                        
                    *maxVarNUZ=0.0;
                    //printf("I am here in diff\n");
                    //printf("%lf \t %lf \t %lf\n", *maxVarNUX,*maxVarNUY,*maxVarNUZ);
		}
	}
	return(0);
}

