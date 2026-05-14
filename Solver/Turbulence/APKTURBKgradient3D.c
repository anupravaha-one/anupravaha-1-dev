#include "APKTURBKgradient3D.h"

double Kgradient3D(Base *base,Zone *zone,double *variable,double *slabVariable[6],int INDEX[3])
{
	int i,j,k,NI,NJ,NK,P,E,W,N,S,T,B,*ptrsize;
	double KE,KW,KN,KS,KT,KB,KP;
	double DKx,DKy,DKz,phiK,pi[3],*deltainv;
	
	ptrsize = zone->Size;
	deltainv = zone->ZGeom->DeltaInverse; // Delta inverse matrix
	
	NI=*(ptrsize+0);
	NJ=*(ptrsize+1);
	NK=*(ptrsize+2);
	k=INDEX[0];
	j=INDEX[1];
	i=INDEX[2];
	P=(NI-1)*(NJ-1)*k+(NI-1)*j+i; 
	E=(NI-1)*(NJ-1)*k+(NI-1)*j+(i+1);
	W=(NI-1)*(NJ-1)*k+(NI-1)*j+(i-1);
	N=(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i;
	S=(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i;
	T=(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i;
	B=(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i;
	
	KP=sqrt(*(variable+P));
	
	if(i==(NI-2))
	{
		KE=sqrt(*(slabVariable[0]+(NJ-1)*k+j));
	}
	else
	{
		KE=sqrt(*(variable+E));
	}
	if(i==0)
	{
		KW=sqrt(*(slabVariable[1]+(NJ-1)*k+j));
	}
	else
	{
		KW=sqrt(*(variable+W));
	}
	if(j==(NJ-2))
	{
		KN=sqrt(*(slabVariable[2]+(NI-1)*k+i));
	}
	else
	{
		KN=sqrt(*(variable+N));
	}
	if(j==0)
	{
		KS=sqrt(*(slabVariable[3]+(NI-1)*k+i));
	}
	else
	{
		KS=sqrt(*(variable+S));
	}
	if(k==0)
	{
		KB=sqrt(*(slabVariable[5]+(NI-1)*j+i));
	}
	else
	{
		KB=sqrt(*(variable+B));
	}
	if(k==(NK-2))
	{
		KT=sqrt(*(slabVariable[4]+(NI-1)*j+i));
	}
	else
	{
		KT=sqrt(*(variable+T));
	}

	pi[0]=KE-KW;
	pi[1]=KN-KS;
	pi[2]=KT-KB;
	
	DKx=*((deltainv+P*9)+0) * pi[0] + *((deltainv+P*9)+1) * pi[1] + *((deltainv+P*9)+2) * pi[2];
	DKy=*((deltainv+P*9)+3) * pi[0] + *((deltainv+P*9)+4) * pi[1] + *((deltainv+P*9)+5) * pi[2];
	DKz=*((deltainv+P*9)+6) * pi[0] + *((deltainv+P*9)+7) * pi[1] + *((deltainv+P*9)+8) * pi[2];
	
	phiK = DKx*DKx + DKy*DKy + DKz*DKz;
	return(phiK);
}
