#include "APKInterpolateEfield.h"

int InterpolateEfield(Base *base,double XP,double YP,double ZP)
{
	int i,j,k,ip,jp,kp,Particleindex;
	
	int NI,NJ,NK,indexB,indexT,indexE,indexW,indexS,indexN,indexP,*ptrSize,Index;
	double *coordX,*coordY,*coordZ,centroid[3],*SurfArea;
	double *ptrPoint1X,*ptrPoint1Y,*ptrPoint1Z;
	double *domainVolume,*domainEfieldX,*domainEfieldY,*domainEfieldZ;
	double *slabVolume[6],*slabEfieldX[6],*slabEfieldY[6],*slabEfieldZ[6];
	double defxdx,defxdy,defxdz,defydx,defydy,defydz,defzdx,defzdy,defzdz;
	Zone *desiredZone;
	double efxw,efyw,efzw,efxe,efye,efze,efxn,efyn,efzn,efxs,efys,efzs,efxt,efyt,efzt,efxb,efyb,efzb;
	double efx_P,efx_E,efx_W,efx_S,efx_N,efx_T,efx_B,efy_P,efy_E,efy_W,efy_S,efy_N,efy_T,efy_B,efz_P,efz_E,efz_W,efz_S,efz_N,efz_T,efz_B;
	double vol_P,vol_E,vol_W,vol_N,vol_S,vol_T,vol_B;
	FILE *ActEfield;
	ActEfield=fopen("SOLUTIONS/ActualEfield.dat","w");
	
//	printf("i m inside interpolate Efield\n");
	
	desiredZone=((base->BaseParameters)->ParticlePosition->zone);
	ip=((base->BaseParameters)->ParticlePosition->i);
	jp=((base->BaseParameters)->ParticlePosition->j);
	kp=((base->BaseParameters)->ParticlePosition->k);
	SurfArea=desiredZone->ZGeom->SurfaceArea;
	
	NI=*((desiredZone)->Size+0);
	NJ=*((desiredZone)->Size+1);
	NK=*((desiredZone)->Size+2);
	
//	printf("ip=%d\tjp=%d\tkp=%d\tzonename=%s\n",ip,jp,kp,desiredZone->zonename);
//	printf("EfieldXCentroid=%12.8lf\tEfieldYCentroid=%12.8lf\tEfieldZCentroid=%12.8lf\n",*((desiredZone->Flowsol)->EfieldX+((NJ-1)*(NI-1)*kp+(NI-1)*jp+ip)),*((desiredZone->Flowsol)->EfieldY+((NJ-1)*(NI-1)*kp+(NI-1)*jp+ip)),*((desiredZone->Flowsol)->EfieldZ+((NJ-1)*(NI-1)*kp+(NI-1)*jp+ip)));
	
	domainEfieldX=MakeDomainArray(desiredZone,11);
	domainEfieldY=MakeDomainArray(desiredZone,12);
	domainEfieldZ=MakeDomainArray(desiredZone,13);
	domainVolume=MakeDomainArray(desiredZone,5);

	MakeSlabPtrArray(desiredZone,slabEfieldX,11);
	MakeSlabPtrArray(desiredZone,slabEfieldY,12);
	MakeSlabPtrArray(desiredZone,slabEfieldZ,13);
	MakeSlabPtrArray(desiredZone,slabVolume,5);
		
	indexP=(NJ-1)*(NI-1)*kp+(NI-1)*jp+ip;
	indexE=(NJ-1)*(NI-1)*kp+(NI-1)*jp+(ip+1);
	indexW=(NJ-1)*(NI-1)*kp+(NI-1)*jp+(ip-1);
	indexN=(NJ-1)*(NI-1)*kp+(NI-1)*(jp+1)+ip;
	indexS=(NJ-1)*(NI-1)*kp+(NI-1)*(jp-1)+ip;
	indexT=(NJ-1)*(NI-1)*(kp+1)+(NI-1)*jp+ip;
	indexB=(NJ-1)*(NI-1)*(kp-1)+(NI-1)*jp+ip;
			
	efx_P=*(domainEfieldX+indexP);
	efy_P=*(domainEfieldY+indexP);
	efz_P=*(domainEfieldZ+indexP);
	vol_P=*(domainVolume+indexP);
	
	if(ip==(NI-2))
	{
		efx_E=*(slabEfieldX[0]+(NJ-1)*kp+jp);
		efy_E=*(slabEfieldY[0]+(NJ-1)*kp+jp);
		efz_E=*(slabEfieldZ[0]+(NJ-1)*kp+jp); 
		vol_E=*(slabVolume[0]+(NJ-1)*kp+jp);
	}
	else
	{
		efx_E=*(domainEfieldX+indexE);
		efy_E=*(domainEfieldY+indexE);
		efz_E=*(domainEfieldZ+indexE);
		vol_E=*(domainVolume+indexE);
	}

	if(ip==0)
	{
		efx_W =*(slabEfieldX[1]+(NJ-1)*kp+jp);
		efy_W =*(slabEfieldY[1]+(NJ-1)*kp+jp);
		efz_W =*(slabEfieldZ[1]+(NJ-1)*kp+jp);
		vol_W =*(slabVolume[1]+(NJ-1)*kp+jp);
	}
	else
	{
		efx_W=*(domainEfieldX+indexW);
		efy_W=*(domainEfieldY+indexW);
		efz_W=*(domainEfieldZ+indexW);
		vol_W=*(domainVolume+indexW);
	}

	if(jp==(NJ-2))
	{
		efx_N=*(slabEfieldX[2]+(NI-1)*kp+ip);
		efy_N=*(slabEfieldY[2]+(NI-1)*kp+ip);
		efz_N=*(slabEfieldZ[2]+(NI-1)*kp+ip);
		vol_N=*(slabVolume[2]+(NI-1)*kp+ip);
	}
	else
	{
		efx_N=*(domainEfieldX+indexN);
		efy_N=*(domainEfieldY+indexN);
		efz_N=*(domainEfieldZ+indexN);
		vol_N=*(domainVolume+indexN);
	}

	if(jp==0)
	{
		efx_S=*(slabEfieldX[3]+(NI-1)*kp+ip);
		efy_S=*(slabEfieldY[3]+(NI-1)*kp+ip);
		efz_S=*(slabEfieldZ[3]+(NI-1)*kp+ip);
		vol_S=*(slabVolume[3]+(NI-1)*kp+ip);
	}
	else
	{
		efx_S=*(domainEfieldX+indexS);
		efy_S=*(domainEfieldY+indexS);
		efz_S=*(domainEfieldZ+indexS);
		vol_S=*(domainVolume+indexS);
	}

	if(kp==(NK-2))
	{
		efx_T=*(slabEfieldX[4]+(NI-1)*jp+ip);
		efy_T=*(slabEfieldY[4]+(NI-1)*jp+ip);
		efz_T=*(slabEfieldZ[4]+(NI-1)*jp+ip);
		vol_T=*(slabVolume[4]+(NI-1)*jp+ip);
	}
	else
	{
		efx_T=*(domainEfieldX+indexT);
		efy_T=*(domainEfieldY+indexT);
		efz_T=*(domainEfieldZ+indexT);
		vol_T=*(domainVolume+indexT);
	}

	if(kp==0)
	{
		efx_B=*(slabEfieldX[5]+(NI-1)*jp+ip);
		efy_B=*(slabEfieldY[5]+(NI-1)*jp+ip);
		efz_B=*(slabEfieldZ[5]+(NI-1)*jp+ip);
		vol_B=*(slabVolume[5]+(NI-1)*jp+ip);
	}
	else
	{
		efx_B=*(domainEfieldX+indexB);
		efy_B=*(domainEfieldY+indexB);
		efz_B=*(domainEfieldZ+indexB);
		vol_B=*(domainVolume+indexB);
	}

	
	efxe=((efx_E*vol_P+efx_P*vol_E)/(vol_P+vol_E));
	efxw=((efx_W*vol_P+efx_P*vol_W)/(vol_P+vol_W));
	efxs=((efx_S*vol_P+efx_P*vol_S)/(vol_P+vol_S));
	efxn=((efx_N*vol_P+efx_P*vol_N)/(vol_P+vol_N));
	efxt=((efx_T*vol_P+efx_P*vol_T)/(vol_P+vol_T));
	efxb=((efx_B*vol_P+efx_P*vol_B)/(vol_P+vol_B));
	
	efye=((efy_E*vol_P+efy_P*vol_E)/(vol_P+vol_E));
	efyw=((efy_W*vol_P+efy_P*vol_W)/(vol_P+vol_W));
	efys=((efy_S*vol_P+efy_P*vol_S)/(vol_P+vol_S));
	efyn=((efy_N*vol_P+efy_P*vol_N)/(vol_P+vol_N));
	efyt=((efy_T*vol_P+efy_P*vol_T)/(vol_P+vol_T));
	efyb=((efy_B*vol_P+efy_P*vol_B)/(vol_P+vol_B));
	
	efze=((efz_E*vol_P+efz_P*vol_E)/(vol_P+vol_E));
	efzw=((efz_W*vol_P+efz_P*vol_W)/(vol_P+vol_W));
	efzs=((efz_S*vol_P+efz_P*vol_S)/(vol_P+vol_S));
	efzn=((efz_N*vol_P+efz_P*vol_N)/(vol_P+vol_N));
	efzt=((efz_T*vol_P+efz_P*vol_T)/(vol_P+vol_T));
	efzb=((efz_B*vol_P+efz_P*vol_B)/(vol_P+vol_B));
	
//	printf("efxe efxw efxn efxs =\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\n",efxe,efxw,efxn,efxs,efxt,efxb);
//	printf("efye efyw efys efyn =\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\n",efye,efyw,efys,efyn,efyt,efyb);
//	printf("efze efzw efzs efzn =\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\n",efze,efzw,efzs,efzn,efzt,efzb);
			
	
	Particleindex=((NJ-1)*(NI-1)*kp+(NI-1)*jp+ip)*6*3;
	
	defxdx=(efxw*(*(SurfArea+Particleindex+3*0+0))+efxe*(*(SurfArea+Particleindex+3*1+0))+efxs*(*(SurfArea+Particleindex+3*2+0))+efxn*(*(SurfArea+Particleindex+3*3+0))+efxb*(*(SurfArea+Particleindex+3*4+0))+efxt*(*(SurfArea+Particleindex+3*5+0)))/(vol_P);
		
	defxdy=(efxw*(*(SurfArea+Particleindex+3*0+1))+efxe*(*(SurfArea+Particleindex+3*1+1))+efxs*(*(SurfArea+Particleindex+3*2+1))+efxn*(*(SurfArea+Particleindex+3*3+1))+efxb*(*(SurfArea+Particleindex+3*4+1))+efxt*(*(SurfArea+Particleindex+3*5+1)))/(vol_P);
		
	//defxdz=(efxw*(*(SurfArea+Particleindex+3*0+2))+efxe*(*(SurfArea+Particleindex+3*1+2))+efxs*(*(SurfArea+Particleindex+3*2+2))+efxn*(*(SurfArea+Particleindex+3*3+2))+efxb*(*(SurfArea+Particleindex+3*4+2))+efxt*(*(SurfArea+Particleindex+3*5+2)))/(vol_P);
		
	defydx=(efyw*(*(SurfArea+Particleindex+3*0+0))+efye*(*(SurfArea+Particleindex+3*1+0))+efys*(*(SurfArea+Particleindex+3*2+0))+efyn*(*(SurfArea+Particleindex+3*3+0))+efyb*(*(SurfArea+Particleindex+3*4+0))+efyt*(*(SurfArea+Particleindex+3*5+0)))/(vol_P);
		
	defydy=(efyw*(*(SurfArea+Particleindex+3*0+1))+efye*(*(SurfArea+Particleindex+3*1+1))+efys*(*(SurfArea+Particleindex+3*2+1))+efyn*(*(SurfArea+Particleindex+3*3+1))+efyb*(*(SurfArea+Particleindex+3*4+1))+efyt*(*(SurfArea+Particleindex+3*5+1)))/(vol_P);
		
	//defydz=(efyw*(*(SurfArea+Particleindex+3*0+2))+efye*(*(SurfArea+Particleindex+3*1+2))+efys*(*(SurfArea+Particleindex+3*2+2))+efyn*(*(SurfArea+Particleindex+3*3+2))+efyb*(*(SurfArea+Particleindex+3*4+2))+efyt*(*(SurfArea+Particleindex+3*5+2)))/(vol_P);
		
	//defzdx=(efzw*(*(SurfArea+Particleindex+3*0+0))+efze*(*(SurfArea+Particleindex+3*1+0))+efzs*(*(SurfArea+Particleindex+3*2+0))+efzn*(*(SurfArea+Particleindex+3*3+0))+efzb*(*(SurfArea+Particleindex+3*4+0))+efzt*(*(SurfArea+Particleindex+3*5+0)))/(vol_P);
		
	//defzdy=(efzw*(*(SurfArea+Particleindex+3*0+1))+efze*(*(SurfArea+Particleindex+3*1+1))+efzs*(*(SurfArea+Particleindex+3*2+1))+efzn*(*(SurfArea+Particleindex+3*3+1))+efzb*(*(SurfArea+Particleindex+3*4+1))+efzt*(*(SurfArea+Particleindex+3*5+1)))/(vol_P);
		
	//defzdz=(efzw*(*(SurfArea+Particleindex+3*0+2))+efze*(*(SurfArea+Particleindex+3*1+2))+efzs*(*(SurfArea+Particleindex+3*2+2))+efzn*(*(SurfArea+Particleindex+3*3+2))+efzb*(*(SurfArea+Particleindex+3*4+2))+efzt*(*(SurfArea+Particleindex+3*5+2)))/(vol_P);
		
	defxdz=defydz=defzdx=defzdy=defzdz=0.0;
	
//	printf("defxdx,defxdy,defxdz,defydx,defydy,defydz,defzdx,defzdy,defzdz   =\t%12.9lf\t%12.8lf\t%12.8lf\%12.8lf\%12.8lf\%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\n",defxdx,defxdy,defxdz,defydx,defydy,defydz,defzdx,defzdy,defzdz);			
	i=ip;j=jp;k=kp;

	ptrSize=((desiredZone)->Size);
	NI=*(ptrSize+0);
	NJ=*(ptrSize+1);
	NK=*(ptrSize+2);
			
	coordX=((desiredZone)->CgnsGrid)->X;
	coordY=((desiredZone)->CgnsGrid)->Y; 
	coordZ=((desiredZone)->CgnsGrid)->Z;
			
	ptrPoint1X=coordX+NI*NJ*k+NI*j+i;
	ptrPoint1Y=coordY+NI*NJ*k+NI*j+i;
	ptrPoint1Z=coordZ+NI*NJ*k+NI*j+i; 
	FindCellCentre(ptrSize,ptrPoint1X,ptrPoint1Y,ptrPoint1Z,centroid);
	Index=(NJ-1)*(NI-1)*k+(NI-1)*j+i;
	
	//3D Case
	((base->BaseParameters)->ParticlePosition->EfieldA[0])=*((desiredZone->Flowsol)->EfieldX+Index)+defxdx*(XP-centroid[0])+defxdy*(YP-centroid[1])+defxdz*(ZP-centroid[2]);			
	((base->BaseParameters)->ParticlePosition->EfieldA[1])=*((desiredZone->Flowsol)->EfieldY+Index)+defydx*(XP-centroid[0])+defydy*(YP-centroid[1])+defydz*(ZP-centroid[2]);			
	((base->BaseParameters)->ParticlePosition->EfieldA[2])=*((desiredZone->Flowsol)->EfieldZ+Index)+defzdx*(XP-centroid[0])+defzdy*(YP-centroid[1])+defzdz*(ZP-centroid[2]);
	
//	printf("\nEfXactual=%12.8lf",((base->BaseParameters)->ParticlePosition->EfieldA[0]));
//	printf("\tEfYactual=%12.8lf\t",((base->BaseParameters)->ParticlePosition->EfieldA[1]));
//	printf("\tEfZactual=%12.8lf\n",((base->BaseParameters)->ParticlePosition->EfieldA[2]));
			
	fprintf(ActEfield,"\n%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\n",((base->BaseParameters)->ParticlePosition->EfieldA[0]),((base->BaseParameters)->ParticlePosition->EfieldA[1]),((base->BaseParameters)->ParticlePosition->EfieldA[2]),*((desiredZone->Flowsol)->EfieldX+Index),*((desiredZone->Flowsol)->EfieldY+Index),*((desiredZone->Flowsol)->EfieldZ+Index));
			
	fclose(ActEfield);
	
return 0;	
} 
