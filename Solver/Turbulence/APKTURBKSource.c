#include "APKTURBKSource.h"
/** This function is used to caluculate the source term in the RANS momentum equations which contains "k" term*/

double KSource(Base *base,Zone *zone,int index[3],int variableTag)
{
	int i,j,k,NI,NJ,NK,P,E,W,N,S,T,B,*ptrsize;
	double *rho,*slabrho[6],*K,*slabK[6],*volume,*slabVolume[6],*DispVolfrac,*slabDispVolfrac[6];
	double RE,RW,RN,RS,RT,RB,RP,KE,KW,KN,KS,KT,KB,KP,VE,VW,VN,VS,VT,VB,VP;
	double volP,volW,volE,volS,volN,volB,volT;
	double Re,Rw,Rs,Rn,Rt,Rb,Ke,Kw,Ks,Kn,Kt,Kb,Ve,Vw,Vs,Vn,Vt,Vb;
	double *sArea;
	double *Surf_indexP,grad[3],grad1[3];
	double variable;
	
	ptrsize = zone->Size;
	
	rho=MakeDomainArray(zone,710);
	MakeSlabPtrArray(zone,slabrho,710);

	K=MakeDomainArray(zone,600);
	MakeSlabPtrArray(zone,slabK,600);
	
	volume=MakeDomainArray(zone,5);
	MakeSlabPtrArray(zone,slabVolume,5);
	
	DispVolfrac=MakeDomainArray(zone,35);
	MakeSlabPtrArray(zone,slabDispVolfrac,35);
	
	sArea=zone->ZGeom->SurfaceArea;
	
	NI=*(ptrsize+0);
	NJ=*(ptrsize+1);
	NK=*(ptrsize+2);
	
	k=index[0];
	j=index[1];
	i=index[2];
	
	P=(NI-1)*(NJ-1)*k+(NI-1)*j+i; 
	E=(NI-1)*(NJ-1)*k+(NI-1)*j+(i+1);
	W=(NI-1)*(NJ-1)*k+(NI-1)*j+(i-1);
	N=(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i;
	S=(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i;
	T=(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i;
	B=(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i;
	
	RP=*(rho+P);
	KP=*(K+P);
	volP=*(volume+P);
	VP=*(DispVolfrac+P);
		
	if(i==(NI-2))
	{
		RE=*(slabrho[0]+(NJ-1)*k+j);
		KE=*(slabK[0]+(NJ-1)*k+j);
		volE=*(slabVolume[0]+(NJ-1)*k+j);
		VE=*(slabDispVolfrac[0]+(NJ-1)*k+j);
	}
	else
	{
		RE=*(rho+E);
		KE=*(K+E);
		volE=*(volume+E);
		VE=*(DispVolfrac+E);
	}
	if(i==0)
	{
		RW=*(slabrho[1]+(NJ-1)*k+j);
		KW=*(slabK[1]+(NJ-1)*k+j);
		volW =*(slabVolume[1]+(NJ-1)*k+j);
		VW=*(slabDispVolfrac[1]+(NJ-1)*k+j);
	}
	else
	{
		RW=*(rho+W);
		KW=*(K+W);
		volW=*(volume+W);
		VW=*(DispVolfrac+W);
	}
	if(j==(NJ-2))
	{
		RN=*(slabrho[2]+(NI-1)*k+i);
		KN=*(slabK[2]+(NI-1)*k+i);
		volN=*(slabVolume[2]+(NI-1)*k+i);
		VN=*(slabDispVolfrac[2]+(NJ-1)*k+j);
	}
	else
	{
		RN=*(rho+N);
		KN=*(K+N);
		volN=*(volume+N);
		VN=*(DispVolfrac+N);
	}
	if(j==0)
	{
		RS=*(slabrho[3]+(NI-1)*k+i);
		KS=*(slabK[3]+(NI-1)*k+i);
		volS=*(slabVolume[3]+(NI-1)*k+i);
		VS=*(slabDispVolfrac[3]+(NJ-1)*k+j);
	}
	else
	{
		RS=*(rho+S);
		KS=*(K+S);
		volS=*(volume+S);
		VS=*(DispVolfrac+S);
	}
	if(k==0)
	{
		RB=*(slabrho[5]+(NI-1)*j+i);
		KB=*(slabK[5]+(NI-1)*j+i);
		volB=*(slabVolume[5]+(NI-1)*j+i);
		VB=*(slabDispVolfrac[5]+(NJ-1)*k+j);
		
	}
	else
	{
		RB=*(rho+B);
		KB=*(K+B);
		volB=*(volume+B);
		VB=*(DispVolfrac+B);
	}
	if(k==(NK-2))
	{
		RT=*(slabrho[4]+(NI-1)*j+i);
		KT=*(slabK[4]+(NI-1)*j+i);
		volT=*(slabVolume[4]+(NI-1)*j+i);
		VT=*(slabDispVolfrac[4]+(NJ-1)*k+j);
	}
	else
	{
		RT=*(rho+T);
		KT=*(K+T);
		volT=*(volume+T);
		VT=*(DispVolfrac+T);
		
	}
	//printf("\nVP=%lf,VE=%lf,VW=%lf,VN=%lf,VS=%lf,VT=%lf,VB=%lf--[Zone %s]",VP,VE,VW,VN,VS,VT,VB,zone->zonename);

	Re=((RE*volP+RP*volE)/(volP+volE));
	Rw=((RW*volP+RP*volW)/(volP+volW));
	Rs=((RS*volP+RP*volS)/(volP+volS));
	Rn=((RN*volP+RP*volN)/(volP+volN));
	Rt=((RT*volP+RP*volT)/(volP+volT));
	Rb=((RB*volP+RP*volB)/(volP+volB));

	Ke=((KE*volP+KP*volE)/(volP+volE));
	Kw=((KW*volP+KP*volW)/(volP+volW));
	Ks=((KS*volP+KP*volS)/(volP+volS));
	Kn=((KN*volP+KP*volN)/(volP+volN));
	Kt=((KT*volP+KP*volT)/(volP+volT));
	Kb=((KB*volP+KP*volB)/(volP+volB));

	Ve=((VE*volP+VP*volE)/(volP+volE));
	Vw=((VW*volP+VP*volW)/(volP+volW));
	Vs=((VS*volP+VP*volS)/(volP+volS));
	Vn=((VN*volP+VP*volN)/(volP+volN));
	Vt=((VT*volP+VP*volT)/(volP+volT));
	Vb=((VB*volP+VP*volB)/(volP+volB));
	
	//printf("\nRe=%lf,Rw=%lf,Rn=%lf,Rs=%lf,Rt=%lf,Rb=%lf--[Zone %s]",Ve,Vw,Vn,Vs,Vt,Vb,zone->zonename);


	Surf_indexP = (sArea+P*6*3);
	
	//printf("\ne=%lf,w=%lf,n=%lf,s=%lf,t=%lf,b=%lf",*(Surf_indexP),*(Surf_indexP+3),*(Surf_indexP+6),*(Surf_indexP+9),*(Surf_indexP+12),*(Surf_indexP+15));
	//West direction, x-component+East direction, x-component+South direction, x-component+North direction,x-component+Bottom direction, x-component+Top direction, x-component
	grad[0]=Rw*Kw*(*(Surf_indexP))+Re*Ke*(*(Surf_indexP+3))+Rs*Ks*(*(Surf_indexP+6))+Rn*Kn*(*(Surf_indexP+9))+
			Rb*Kb*(*(Surf_indexP+12))+Rt*Kt*(*(Surf_indexP+15));
	
	
	grad[0]=(-2.0/3)*grad[0];
	
	//West direction, y-component//East direction, y-component//South direction, y-component//North direction, y-component//Bottom direction, y-component//Top direction, y-component
	grad[1]=Rw*Kw*(*(Surf_indexP+1))+Re*Ke*(*(Surf_indexP+4))+Rs*Ks*(*(Surf_indexP+7))+Rn*Kn*(*(Surf_indexP+10))+
			Rb*Kb*(*(Surf_indexP+13))+Rt*Kt*(*(Surf_indexP+16));
	
	grad[1]=(-2.0/3)*grad[1];
	
	if(base->BaseParameters->flag2d==0)
	{
		//West direction, z-component+East direction, z-component+South direction, z-component+North direction, z-component+Bottom direction, z-component+Top direction, z-component
		grad[2]=Rw*Kw*(*(Surf_indexP+2))+Re*Ke*(*(Surf_indexP+5))+Rs*Ks*(*(Surf_indexP+8))+Rn*Kn*(*(Surf_indexP+11))+
				Rb*Kb*(*(Surf_indexP+14))+Rt*Kt*(*(Surf_indexP+17));
		grad[2] = (-2.0/3) *grad[2];
	}

	//West direction, x-component+East direction, x-component+South direction, x-component+North direction,x-component+Bottom direction, x-component+Top direction, x-component
	grad1[0]=Rw*Kw*Vw*(*(Surf_indexP))+Re*Ke*Ve*(*(Surf_indexP+3))+Rs*Ks*Vs*(*(Surf_indexP+6))+Rn*Kn*Vn*(*(Surf_indexP+9))+
			Rb*Kb*Vb*(*(Surf_indexP+12))+Rt*Kt*Vt*(*(Surf_indexP+15));

	grad1[0]=(-2.0/3)*grad1[0];
	
	//West direction, y-component//East direction, y-component//South direction, y-component//North direction, y-component//Bottom direction, y-component//Top direction, y-component
	grad1[1]=Rw*Kw*Vw*(*(Surf_indexP+1))+Re*Ke*Ve*(*(Surf_indexP+4))+Rs*Ks*Vs*(*(Surf_indexP+7))+Rn*Kn*Vn*(*(Surf_indexP+10))+
			Rb*Kb*Vb*(*(Surf_indexP+13))+Rt*Kt*Vt*(*(Surf_indexP+16));
	
	grad1[1]=(-2.0/3)*grad1[1];
	
	if(base->BaseParameters->flag2d==0)
	{
		//West direction, z-component+East direction, z-component+South direction, z-component+North direction, z-component+Bottom direction, z-component+Top direction, z-component
		grad1[2]=Rw*Kw*Vw*(*(Surf_indexP+2))+Re*Ke*Ve*(*(Surf_indexP+5))+Rs*Ks*Vs*(*(Surf_indexP+8))+Rn*Kn*Vn*(*(Surf_indexP+11))+Rb*Kb*Vb*(*(Surf_indexP+14))+Rt*Kt*Vt*(*(Surf_indexP+17));
		grad1[2] = (-2.0/3) *grad1[2];
	}

	if(variableTag==3100||variableTag==31)
		variable=grad[0];
	else if(variableTag==3200||variableTag==32)
		variable=grad[1];
	else if(variableTag==3300||variableTag==33)
		variable=grad[2];
	else if(variableTag==44)
		variable=grad1[0];
	else if(variableTag==45)
		variable=grad1[1];
	else if (variableTag==46)
		variable=grad1[2];
	//printf("\nVariable=%lf",variable);

	return(variable);
}