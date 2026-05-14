#include "APKFindmaXConv.h"
/**
This function, FindmaxConv() is used to find maximum value of U/deltaX, V/deltaY and W/deltaZ in the domain. These values
are used in the determination of deltaT(time) by using CFL condition.
*/

int FindmaXConv(Zone *zone,int *size,double *maxVarX,double *maxVarY,double *maxVarZ,int flag)
{
	int Ni,Nj,Nk,i,j,k;
	double u,v,w,comparX,comparZ,comparY;
	double *x,*y,*z,*varU,*varV,*varW;
        double Delx,Dely,Delz,Del;
	
	
	Ni=*(size+0);
	Nj=*(size+1);
	Nk=*(size+2);
	
	x=zone->CgnsGrid->X;
	y=zone->CgnsGrid->Y;
	z=zone->CgnsGrid->Z;
	   
	varU=zone->Flowsol->u;
	varV=zone->Flowsol->v;
	varW=zone->Flowsol->w;
	//volume=zone->ZGeom->Volume;
	
	//maxVar=(sqrt((*varU)*(*varU)+(*varV)*(*varV)+(*varW)*(*varW)))/(pow(*volume,0.33333));
	
	
	if(flag==0)
	{
		*maxVarX=0.0;
		*maxVarY=0.0;
		*maxVarZ=0.0;
                
		for(k=0;k<Nk-1;k++)
		for(j=0;j<Nj-1;j++)
		for(i=0;i<Ni-1;i++)
		{
                        Delx = (*(x+(Ni*Nj*k)+(Ni*j)+i+1)) - (*(x+(Ni*Nj*k)+(Ni*j)+i)) ;
                        Dely = (*(y+(Ni*Nj*k)+(Ni*j)+i+1)) - (*(y+(Ni*Nj*k)+(Ni*j)+i)) ;
                        Delz = (*(z+(Ni*Nj*k)+(Ni*j)+i+1)) - (*(z+(Ni*Nj*k)+(Ni*j)+i)) ;
                        Del  = sqrt(Delx*Delx+Dely*Dely+Delz*Delz);
                        
			u=(*(varU+(Nj-1)*(Ni-1)*k+(Ni-1)*j+i));

			comparX=fabs(u/Del);

			if(*maxVarX<comparX)
			*maxVarX=comparX;

			v=fabs(*(varV+(Nj-1)*(Ni-1)*k+(Ni-1)*j+i));

			comparY=fabs(v/Del);

			if(*maxVarY<comparY)
			*maxVarY=comparY;

			w=fabs(*(varW+(Nj-1)*(Ni-1)*k+(Ni-1)*j+i));

			comparZ=fabs(w/Del);

			if(*maxVarZ<comparZ)
			*maxVarZ=comparZ;

		}
	}
	else
	{
		*maxVarX=0.0;
		*maxVarY=0.0;
                
		for(k=0;k<Nk-1;k++)
		for(j=0;j<Nj-1;j++)
		for(i=0;i<Ni-1;i++)
		{
                        Delx = (*(x+(Ni*Nj*k)+(Ni*j)+i+1)) - (*(x+(Ni*Nj*k)+(Ni*j)+i)) ;
                        Dely = (*(y+(Ni*Nj*k)+(Ni*j)+i+1)) - (*(y+(Ni*Nj*k)+(Ni*j)+i)) ;
                        Del  = sqrt(Delx*Delx+Dely*Dely);
                        
			u=(*(varU+(Nj-1)*(Ni-1)*k+(Ni-1)*j+i));

			comparX=fabs(u/Del);

			if(*maxVarX<comparX)
			*maxVarX=comparX;

			v=fabs(*(varV+(Nj-1)*(Ni-1)*k+(Ni-1)*j+i));

			comparY=fabs(v/Del);

			if(*maxVarY<comparY)
			*maxVarY=comparY;

			*maxVarZ=0.00;
			//printf("I am here in conv\n");
		}
	}
	return(0);
}

