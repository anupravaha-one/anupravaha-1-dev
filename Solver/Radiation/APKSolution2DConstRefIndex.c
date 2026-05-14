#include "APKSolution2DConstRefIndex.h"

/** The function Solution2DConstRefIndex() is the main control file for calculating the coefficient of Gauss Seidel formulation. This function 
passes the control to the functions RadiationSolutionInt2D(), IntExtremeFace2D() and IntExtremeCell2D() depending on the cell location within the zone. If the cells are completely
interior cells then RadiationSolutionInt2D() is called. If a cell has only one boundary or an interface among any of its four 
faces(since 2D), then IntExtremeFace2D() is called. If a cell has more than one boundary or an interface among any of its
four faces(since 2D), then IntExtremeCell2D() is called (This condition occurs for completely corner cells). At the end it calculate all the residuals
depends on the different conditions and calculate the new Intensity.*/

int Solution2DConstRefIndex(Base *base, Zone *zone, double *NewInt, double *OldInt, double alpha)
	{
	 int i,j,k,l,m,li,mi,ni,nj,nk,nt,np,cellindex,face,zflag,flagcv,zflagT,choice,index[5];
	 
	 double *slabInt[6],*slabvolume[6],Esource,*Vsource,PI,Ssource,*Temperature;
	 double *Intface, *Intdiag, Irra, Diag[6],diagonal, Relax=1.0,volP;
	 double *volume, *Kp, *Sp, *omega, *refdex, sumF, AX, RS,residue;
	 
	 PI=4.0*atan(1.0);
	  
	 ni=*(zone->Size+0);
	 nj=*(zone->Size+1);
	 nk=*(zone->Size+2);
	 nt=base->BaseParameters->ntheta;
	 np=base->BaseParameters->nphi;
	 
	 zflagT=zone->flagZoneTemp;
	 zflag=zone->flagRadHeatGeneration;
	 flagcv=zone->flagRadHeatConVar;
	 choice=base->BaseParameters->Scatteringchoice;
	 
	 Kp=zone->Flowsol->kappa;
	 Sp=zone->Flowsol->sigma;
	 omega=zone->ZGeom->Omega;
	 refdex=zone->Flowsol->RefIndex;
	 Vsource=zone->Flowsol->RadEnergysource;
	 volume=zone->ZGeom->Volume;
	 if(base->BaseParameters->flagRadiation==1)
	 Temperature=zone->Flowsol->RadGasTemp;
	 if(base->BaseParameters->flagflowRadiation==1)
	 {
	 Temperature=zone->PrevFlowsol->T;
	 zflagT=1;
	 }
	 //printf("Temp=%e\n",*Temperature);
	 
	 
	 Intface=(double *)malloc((ni-1)*(nj-1)*(nk-1)*6*nt*np*sizeof(double));
	 Intdiag=(double *)malloc((ni-1)*(nj-1)*(nk-1)*nt*np*sizeof(double));
	 
	 MakeSlabPtrArray(zone,slabInt,16); // Will Return variable (Intensity) slab pointer
	 MakeSlabPtrArray(zone,slabvolume,5); // Will Return volume slab pointer
	
	 
	 for(l=0;l<nt;l++)
	 for(m=0;m<np;m++)
		for(k=0;k<(nk-1);k++)
		for(j=0;j<(nj-1);j++)
		for(i=0;i<(ni-1);i++)
		  {
		  cellindex=(ni-1)*(nj-1)*k+(ni-1)*j+i;
		  index[0]=i,index[1]=j,index[2]=k,index[3]=l,index[4]=m;
		  volP=*(volume+cellindex);
		  
		  if((i>0&&i<ni-2)&&(j>0&&j<nj-2)&&(k>0&&k<nk-2))
		  	{
			 RadiationSolutionInt2D(base, zone, Intface, Diag, NewInt, OldInt, index, alpha);
			 
			 diagonal=0.0;
			 for(face=0;face<6;face++)
			 	diagonal=diagonal+Diag[face];
			 *(Intdiag+cellindex*nt*np+np*l+m)=diagonal;
			}
		
		else
		{ 
			diagonal=0.0;
		 for(face=0;face<6;face++)
		 	{
		switch(face)
			{
		case 0:
			if(i==0)
			Diag[face]=IntExtremeFace2D(base, zone, Intface, NewInt, OldInt, slabInt, slabvolume, index, alpha, face);
			else
			Diag[face]=IntExtremeCell2D(base, zone, Intface, NewInt, OldInt, index, alpha, face);
		break;
		
		case 1:
			if(i==ni-2)
			Diag[face]=IntExtremeFace2D(base, zone, Intface, NewInt, OldInt, slabInt, slabvolume, index, alpha, face);
			else
			Diag[face]=IntExtremeCell2D(base, zone, Intface, NewInt, OldInt, index, alpha, face);
		break;
		
		case 2:
			if(j==0)
			Diag[face]=IntExtremeFace2D(base, zone, Intface, NewInt, OldInt, slabInt, slabvolume, index, alpha, face);
			else
			Diag[face]=IntExtremeCell2D(base, zone, Intface, NewInt, OldInt, index, alpha, face);
		break;
		
		case 3:
			if(j==nj-2)
			Diag[face]=IntExtremeFace2D(base, zone, Intface, NewInt, OldInt, slabInt, slabvolume, index, alpha, face);
			else
			Diag[face]=IntExtremeCell2D(base, zone, Intface, NewInt, OldInt, index, alpha, face);
		break;
		
		case 4:
		
			if(k==0)
			Diag[face]=IntExtremeFace2D(base, zone, Intface, NewInt, OldInt, slabInt, slabvolume, index, alpha, face);
			else
			Diag[face]=IntExtremeCell2D(base, zone, Intface, NewInt, OldInt, index, alpha, face);
		break;
		
		case 5:
			if(k==nk-2)
			Diag[face]=IntExtremeFace2D(base, zone, Intface, NewInt, OldInt, slabInt, slabvolume, index, alpha, face);
			else
			Diag[face]=IntExtremeCell2D(base, zone, Intface, NewInt, OldInt, index, alpha, face);
		break;
			}
		diagonal=diagonal+Diag[face];
		
			}
		*(Intdiag+cellindex*nt*np+np*l+m)=diagonal;
			 
		}
		
		if(((base->BaseParameters->flagAbsorption)==1)&&((base->BaseParameters->flagScattering)==1))
	 	{
		 	if((base->BaseParameters->flagRadSource)==1)
		 		{
		 		Irra=0.0;
		  		for(li=0;li<nt;li++)
	 	 		for(mi=0;mi<np;mi++)
		  		Irra=Irra+*(OldInt+((ni-1)*(nj-1)*k+(ni-1)*j+i)*nt*np+np*li+mi)*(*(omega+np*li+mi));
				     
				Esource=pow(*(refdex),2)*((*Kp)*Irra+zflag*(*(Vsource+flagcv*cellindex)))*volP*(*(omega+np*l+m))/(4.0*PI);
				}
			if((base->BaseParameters->flagTemp)==1)
				{
		 		Esource=5.67e-8*pow(*(refdex),2)*(*Kp)*pow(zflagT*(*(Temperature+flagcv*cellindex)),4)*volP*(*(omega+np*l+m))/PI;
		 		//printf("vol=%e\tomega=%e\n",volP,*(omega+np*l+m));
				}
		
				Ssource=ChoiceScatteringfunction(choice,cellindex,OldInt,l,m,nt,np);
		
			sumF=0.0;
		 	for(face=0;face<6;face++)
		 		sumF=sumF+ *(Intface+cellindex*nt*np*6+nt*np*face+np*l+m);
			
				AX=sumF+(*Kp+*Sp)*(*(NewInt+cellindex*nt*np+np*l+m))*volP*(*(omega+np*l+m));
				RS=Esource+(*Sp)*Ssource*volP/(4.0*PI);
			
				residue=RS-AX;
			
				*(NewInt+cellindex*nt*np+np*l+m)=*(NewInt+cellindex*nt*np+np*l+m) + Relax*residue/(*(Intdiag+cellindex*nt*np+np*l+m)+(*Kp+*Sp)*volP*(*(omega+np*l+m)));
				
				if(*(NewInt+cellindex*nt*np+np*l+m)<0)
				*(NewInt+cellindex*nt*np+np*l+m)=0.0;
		}
		else
		{   
			if((base->BaseParameters->flagAbsorption)==1)
		 	{
		 		if((base->BaseParameters->flagRadSource)==1)
		 		{
		 			Irra=0.0;
		  			for(li=0;li<nt;li++)
	 	  			for(mi=0;mi<np;mi++)
		  			Irra=Irra+*(OldInt+((ni-1)*(nj-1)*k+(ni-1)*j+i)*nt*np+np*li+mi)*(*(omega+np*li+mi));
				     
					Esource=pow(*(refdex),2)*((*Kp)*Irra+zflag*(*(Vsource+flagcv*cellindex)))*volP*(*(omega+np*l+m))/(4.0*PI);
				}
				if((base->BaseParameters->flagTemp)==1)
				{
		 		Esource=5.67e-8*pow(*(refdex),2)*(*Kp)*pow(zflagT*(*(Temperature+flagcv*cellindex)),4)*volP*(*(omega+np*l+m))/PI;
				//if(*(Temperature+flagcv*cellindex)>0)
				//printf("%e\t",*(Temperature+flagcv*cellindex));
		 		//printf("vol=%e\tomega=%e\n",volP,*(omega+np*l+m));
				}
				sumF=0.0;
		 		for(face=0;face<6;face++)
		 			sumF=sumF+ *(Intface+cellindex*nt*np*6+nt*np*face+np*l+m);
			
					AX=sumF+(*Kp)*(*(NewInt+cellindex*nt*np+np*l+m))*volP*(*(omega+np*l+m));
					RS=Esource;
			
					residue=RS-AX;
			
				*(NewInt+cellindex*nt*np+np*l+m)=*(NewInt+cellindex*nt*np+np*l+m) + Relax*residue/(*(Intdiag+cellindex*nt*np+np*l+m)+(*Kp)*volP*(*(omega+np*l+m)));
				
				if(*(NewInt+cellindex*nt*np+np*l+m)<0)
				*(NewInt+cellindex*nt*np+np*l+m)=0.0;
				
			}
					
			if((base->BaseParameters->flagScattering)==1)
				{
					Ssource=ChoiceScatteringfunction(choice,cellindex,OldInt,l,m,nt,np);
		
					sumF=0.0;
		 		for(face=0;face<6;face++)
		 			sumF=sumF+ *(Intface+cellindex*nt*np*6+nt*np*face+np*l+m);
			
					AX=sumF+(*Sp)*(*(NewInt+cellindex*nt*np+np*l+m))*volP*(*(omega+np*l+m));
					RS=(*Sp)*Ssource*volP/(4.0*PI);
			
					residue=RS-AX;
			
				*(NewInt+cellindex*nt*np+np*l+m)=*(NewInt+cellindex*nt*np+np*l+m) + Relax*residue/(*(Intdiag+cellindex*nt*np+np*l+m)+(*Sp)*volP*(*(omega+np*l+m)));
				
				if(*(NewInt+cellindex*nt*np+np*l+m)<0)
				*(NewInt+cellindex*nt*np+np*l+m)=0.0;
				
				}
			
		}
		}
		free(Intface);
		free(Intdiag);
	 
	 return 0;
	}
