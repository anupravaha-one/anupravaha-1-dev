#include "APKInterpolateVelocity.h"

int InterpolateVelocity (Base *base,double XP,double YP,double ZP)
{
	int i,j,k,ip,jp,kp,Index,Particleindex;
//	double up,vp,wp,ua,va,wa;
	
	int NI,NJ,NK,indexB,indexT,indexE,indexW,indexS,indexN,indexP,*ptrSize;
	double *coordX,*coordY,*coordZ,centroid[3],*SurfArea;
	double *ptrPoint1X,*ptrPoint1Y,*ptrPoint1Z;
	double *domainVolume,*domainU,*domainV,*domainW;
	double *slabVolume[6],*slabU[6],*slabV[6],*slabW[6];
	double dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz;
	Zone  *desiredZone;
	
	double uw,vw,ww,ue,ve,we,un,vn,wn,us,vs,ws,ut,vt,wt,ub,vb,wb;
	double u_P,u_E,u_W,u_S,u_N,u_T,u_B,v_P,v_E,v_W,v_S,v_N,v_T,v_B,w_P,w_E,w_W,w_S,w_N,w_T,w_B;
	double vol_P,vol_E,vol_W,vol_N,vol_S,vol_T,vol_B;
	FILE *ActVel;
	ActVel=fopen("SOLUTIONS/AcualVelocity.dat","w");
	
//	printf("i m inside interpolate velocities\n");
	
	desiredZone=((base->BaseParameters)->ParticlePosition->zone);
	ip=((base->BaseParameters)->ParticlePosition->i);
	jp=((base->BaseParameters)->ParticlePosition->j);
	kp=((base->BaseParameters)->ParticlePosition->k);
	SurfArea=desiredZone->ZGeom->SurfaceArea;
	
	NI=*((desiredZone)->Size+0);
	NJ=*((desiredZone)->Size+1);
	NK=*((desiredZone)->Size+2);
	
	domainU=MakeDomainArray(desiredZone,1);
	domainV=MakeDomainArray(desiredZone,2);
	domainW=MakeDomainArray(desiredZone,3);
	domainVolume=MakeDomainArray(desiredZone,5);

	MakeSlabPtrArray(desiredZone,slabU,1);
	MakeSlabPtrArray(desiredZone,slabV,2);
	MakeSlabPtrArray(desiredZone,slabW,3);
	MakeSlabPtrArray(desiredZone,slabVolume,5);
		
	indexP=(NJ-1)*(NI-1)*kp+(NI-1)*jp+ip;
	indexE=(NJ-1)*(NI-1)*kp+(NI-1)*jp+(ip+1);
	indexW=(NJ-1)*(NI-1)*kp+(NI-1)*jp+(ip-1);
	indexN=(NJ-1)*(NI-1)*kp+(NI-1)*(jp+1)+ip;
	indexS=(NJ-1)*(NI-1)*kp+(NI-1)*(jp-1)+ip;
	indexT=(NJ-1)*(NI-1)*(kp+1)+(NI-1)*jp+ip;
	indexB=(NJ-1)*(NI-1)*(kp-1)+(NI-1)*jp+ip;
			
		
	u_P=*(domainU+indexP);
	v_P=*(domainV+indexP);
	w_P=*(domainW+indexP);
	vol_P=*(domainVolume+indexP);
		
	if(ip==(NI-2))
	{
		u_E=*(slabU[0]+(NJ-1)*kp+jp);
		v_E=*(slabV[0]+(NJ-1)*kp+jp);
		w_E=*(slabW[0]+(NJ-1)*kp+jp);
		vol_E=*(slabVolume[0]+(NJ-1)*kp+jp);
	}
	else
	{
		u_E=*(domainU+indexE);
		v_E=*(domainV+indexE);
		w_E=*(domainW+indexE);
		vol_E=*(domainVolume+indexE);
	 }

	if(ip==0)
	{
		u_W =*(slabU[1]+(NJ-1)*kp+jp);
		v_W =*(slabV[1]+(NJ-1)*kp+jp);
		w_W =*(slabW[1]+(NJ-1)*kp+jp);
		vol_W =*(slabVolume[1]+(NJ-1)*kp+jp);
	}
	else
	{
		u_W=*(domainU+indexW);
		v_W=*(domainV+indexW);
		w_W=*(domainW+indexW);
		vol_W=*(domainVolume+indexW);
	}
	if(jp==(NJ-2))
	{
		u_N=*(slabU[2]+(NI-1)*kp+ip);
		v_N=*(slabV[2]+(NI-1)*kp+ip);
		w_N=*(slabW[2]+(NI-1)*kp+ip);
		vol_N=*(slabVolume[2]+(NI-1)*kp+ip);
	}
	else
	{
		u_N=*(domainU+indexN);
		v_N=*(domainV+indexN);
		w_N=*(domainW+indexN);
		vol_N=*(domainVolume+indexN);
	}
	if(jp==0)
	{
		u_S=*(slabU[3]+(NI-1)*kp+ip);
		v_S=*(slabV[3]+(NI-1)*kp+ip);
		w_S=*(slabW[3]+(NI-1)*kp+ip);
		vol_S=*(slabVolume[3]+(NI-1)*kp+ip);
	}
	else
	{
		u_S=*(domainU+indexS);
		v_S=*(domainV+indexS);
		w_S=*(domainW+indexS);
		vol_S=*(domainVolume+indexS);
	}

	if(kp==(NK-2))
	{
		u_T=*(slabU[4]+(NI-1)*jp+ip);
		v_T=*(slabV[4]+(NI-1)*jp+ip);
		w_T=*(slabW[4]+(NI-1)*jp+ip);
		vol_T=*(slabVolume[4]+(NI-1)*jp+ip);
	}
	else
	{
		u_T=*(domainU+indexT);
		v_T=*(domainV+indexT);
		w_T=*(domainW+indexT);
		vol_T=*(domainVolume+indexT);
	}
	if(kp==0)
	{
		u_B=*(slabU[5]+(NI-1)*jp+ip);
		v_B=*(slabV[5]+(NI-1)*jp+ip);
		w_B=*(slabW[5]+(NI-1)*jp+ip);
		vol_B=*(slabVolume[5]+(NI-1)*jp+ip);
	}
	else
	{
		u_B=*(domainU+indexB);
		v_B=*(domainV+indexB);
		w_B=*(domainW+indexB);
		vol_B=*(domainVolume+indexB);
		
	}
	
//	printf("u_p =%12.8lf\tu_e=%12.8lf\tu_w=%12.8lf\tu_n=%12.8lf\tu_s=%12.8lf\tu_t=%12.8lf\tu_b=%12.8lf\n",u_P,u_E,u_W,u_N,u_S,u_T,u_B);
//	printf("v_p =%12.8lf\tv_e=%12.8lf\tv_w=%12.8lf\tv_n=%12.8lf\tv_s=%12.8lf\tv_t=%12.8lf\tv_b=%12.8lf\n",v_P,v_E,v_W,v_N,v_S,v_T,v_B);
//	printf("w_p =%12.8lf\tw_e=%12.8lf\tw_w=%12.8lf\tw_n=%12.8lf\tw_s=%12.8lf\tw_t=%12.8lf\tw_b=%12.8lf\n",w_P,w_E,w_W,w_N,w_S,w_T,w_B);
	ue=((u_E*vol_P+u_P*vol_E)/(vol_P+vol_E));
	uw=((u_W*vol_P+u_P*vol_W)/(vol_P+vol_W));
	us=((u_S*vol_P+u_P*vol_S)/(vol_P+vol_S));
	un=((u_N*vol_P+u_P*vol_N)/(vol_P+vol_N));
	ut=((u_T*vol_P+u_P*vol_T)/(vol_P+vol_T));
	ub=((u_B*vol_P+u_P*vol_B)/(vol_P+vol_B));
		
	ve=((v_E*vol_P+v_P*vol_E)/(vol_P+vol_E));
	vw=((v_W*vol_P+v_P*vol_W)/(vol_P+vol_W));
	vs=((v_S*vol_P+v_P*vol_S)/(vol_P+vol_S));
	vn=((v_N*vol_P+v_P*vol_N)/(vol_P+vol_N));
	vt=((v_T*vol_P+v_P*vol_T)/(vol_P+vol_T));
	vb=((v_B*vol_P+v_P*vol_B)/(vol_P+vol_B));
			
	we=((w_E*vol_P+w_P*vol_E)/(vol_P+vol_E));
	ww=((w_W*vol_P+w_P*vol_W)/(vol_P+vol_W));
	ws=((w_S*vol_P+w_P*vol_S)/(vol_P+vol_S));
	wn=((w_N*vol_P+w_P*vol_N)/(vol_P+vol_N));
	wt=((w_T*vol_P+w_P*vol_T)/(vol_P+vol_T));
	wb=((w_B*vol_P+w_P*vol_B)/(vol_P+vol_B));
		
//	printf("ue uw un us ut ub=\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\n",ue,uw,us,un,ut,ub);
//	printf("ve,vw,vs,vn,vt,vb=\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\n",ve,vw,vs,vn,vt,vb);
//	printf("we,ww,ws,wn,wt,wb=\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\n",we,ww,ws,wn,wt,wb);


	Particleindex=((NJ-1)*(NI-1)*kp+(NI-1)*jp+ip)*6*3;

	//3D case
	dudx=(uw*(*(SurfArea+Particleindex+3*0+0))+ue*(*(SurfArea+Particleindex+3*1+0))+us*(*(SurfArea+Particleindex+3*2+0))+un*(*(SurfArea+Particleindex+3*3+0))+ub*(*(SurfArea+Particleindex+3*4+0))+ut*(*(SurfArea+Particleindex+3*5+0)))/(vol_P);
			
	dudy=(uw*(*(SurfArea+Particleindex+3*0+1))+ue*(*(SurfArea+Particleindex+3*1+1))+us*(*(SurfArea+Particleindex+3*2+1))+un*(*(SurfArea+Particleindex+3*3+1))+ub*(*(SurfArea+Particleindex+3*4+1))+ut*(*(SurfArea+Particleindex+3*5+1)))/(vol_P);
			
	//dudz=(uw*(*(SurfArea+Particleindex+3*0+2))+ue*(*(SurfArea+Particleindex+3*1+2))+us*(*(SurfArea+Particleindex+3*2+2))+un*(*(SurfArea+Particleindex+3*3+2))+ub*(*(SurfArea+Particleindex+3*4+2))+ut*(*(SurfArea+Particleindex+3*5+2)))/(vol_P);
			
	dvdx=(vw*(*(SurfArea+Particleindex+3*0+0))+ve*(*(SurfArea+Particleindex+3*1+0))+vs*(*(SurfArea+Particleindex+3*2+0))+vn*(*(SurfArea+Particleindex+3*3+0))+vb*(*(SurfArea+Particleindex+3*4+0))+vt*(*(SurfArea+Particleindex+3*5+0)))/(vol_P);
			
	dvdy=(vw*(*(SurfArea+Particleindex+3*0+1))+ve*(*(SurfArea+Particleindex+3*1+1))+vs*(*(SurfArea+Particleindex+3*2+1))+vn*(*(SurfArea+Particleindex+3*3+1))+vb*(*(SurfArea+Particleindex+3*4+1))+vt*(*(SurfArea+Particleindex+3*5+1)))/(vol_P);
			
	//dvdz=(vw*(*(SurfArea+Particleindex+3*0+2))+ve*(*(SurfArea+Particleindex+3*1+2))+vs*(*(SurfArea+Particleindex+3*2+2))+vn*(*(SurfArea+Particleindex+3*3+2))+vb*(*(SurfArea+Particleindex+3*4+2))+vt*(*(SurfArea+Particleindex+3*5+2)))/(vol_P);
			
	//dwdx=(ww*(*(SurfArea+Particleindex+3*0+0))+we*(*(SurfArea+Particleindex+3*1+0))+ws*(*(SurfArea+Particleindex+3*2+0))+wn*(*(SurfArea+Particleindex+3*3+0))+wb*(*(SurfArea+Particleindex+3*4+0))+wt*(*(SurfArea+Particleindex+3*5+0)))/(vol_P);
			
	//dwdy=(ww*(*(SurfArea+Particleindex+3*0+1))+we*(*(SurfArea+Particleindex+3*1+1))+ws*(*(SurfArea+Particleindex+3*2+1))+wn*(*(SurfArea+Particleindex+3*3+1))+wb*(*(SurfArea+Particleindex+3*4+1))+wt*(*(SurfArea+Particleindex+3*5+1)))/(vol_P);
			
	//dwdz=(ww*(*(SurfArea+Particleindex+3*0+2))+we*(*(SurfArea+Particleindex+3*1+2))+ws*(*(SurfArea+Particleindex+3*2+2))+wn*(*(SurfArea+Particleindex+3*3+2))+wb*(*(SurfArea+Particleindex+3*4+2))+wt*(*(SurfArea+Particleindex+3*5+2)))/(vol_P);
	
	dudz=dvdz=dwdx=dwdy=dwdz=0.0;	
		
//	printf("dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz =\t%12.9lf\t%12.8lf\t%12.8lf\%12.8lf\%12.8lf\%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\n",dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz);			
	i=ip;j=jp;k=kp;
	ptrSize=(desiredZone)->Size;
	NI=*(ptrSize+0);
	NJ=*(ptrSize+1);
	NK=*(ptrSize+2);
	Index=((NJ-1)*(NI-1)*k+(NI-1)*j+i);
	coordX=((desiredZone)->CgnsGrid)->X;
	coordY=((desiredZone)->CgnsGrid)->Y; 
	coordZ=((desiredZone)->CgnsGrid)->Z;
	ptrPoint1X=coordX+NI*NJ*k+NI*j+i;
	ptrPoint1Y=coordY+NI*NJ*k+NI*j+i;
	ptrPoint1Z=coordZ+NI*NJ*k+NI*j+i;  
	
	FindCellCentre(ptrSize,ptrPoint1X,ptrPoint1Y,ptrPoint1Z,centroid);
		
//	printf("\nup  vp  wp  =\t%12.8lf\t%12.8lf\t%12.8lf\n",*((desiredZone->Flowsol)->u+Index),*((desiredZone->Flowsol)->v+Index),*((desiredZone->Flowsol)->w+Index));
		 
	 //3D Case
	((base->BaseParameters)->ParticlePosition->VelA[0])=*((desiredZone->Flowsol)->u+Index)+dudx*(XP-centroid[0])+dudy*(YP-centroid[1])+dudz*(ZP-centroid[2]);
		
	((base->BaseParameters)->ParticlePosition->VelA[1])=*((desiredZone->Flowsol)->v+Index)+dvdx*(XP-centroid[0])+dvdy*(YP-centroid[1])+dvdz*(ZP-centroid[2]);
		
	((base->BaseParameters)->ParticlePosition->VelA[2])=*((desiredZone->Flowsol)->w+Index)+dwdx*(XP-centroid[0])+dwdy*(YP-centroid[1])+dwdz*(ZP-centroid[2]);
	
		
//	printf("\nuactual=%12.8lf",((base->BaseParameters)->ParticlePosition->VelA[0]));
//	printf("\tvactual=%12.8lf",((base->BaseParameters)->ParticlePosition->VelA[1]));
//	printf("\twactual=%12.8lf\n",((base->BaseParameters)->ParticlePosition->VelA[2]));
		
	fprintf(ActVel,"\n%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf\t%12.8lf",((base->BaseParameters)->ParticlePosition->VelA[0]),((base->BaseParameters)->ParticlePosition->VelA[1]),((base->BaseParameters)->ParticlePosition->VelA[2]),*((desiredZone->Flowsol)->u+Index),*((desiredZone->Flowsol)->v+Index),*((desiredZone->Flowsol)->w+Index));
	
	fclose(ActVel);
	return 0;
}
