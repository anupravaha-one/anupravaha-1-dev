#include "APKSmoothing.h"


double Smoothing(Base *base,Zone *zone,int gdi)
{
//	printf("Entered APKUpdatePressure2d.c\n");
	double residZone=0.0,*pressureDomain,*diffdia;
	double *LHS,*RHS,resid;
	int NI,NJ,NK,i,j,k,index[3],indexP,*ptrSize,flagrho;
	double *slabPressure[6],rho;
	double delta_t,*CoeffN[18],V[18],Vg,omega=1.0;
	FILE *Pomega;
	
	LHS=(double *)malloc(sizeof(double));
	RHS=(double *)malloc(sizeof(double));
	Pomega = fopen("Settings/PressureRelaxation.dat","r");
	fscanf(Pomega,"%lf", &omega);
	fclose(Pomega);
	
	diffdia = (zone->grid+gdi)->DiffusionCoeffPhi_P;
	for(i=0;i<18;i++)
	{
		CoeffN[i] = (zone->grid+gdi)->DiffusionCoeffPhi_nb[i];
	}
	delta_t = base->BaseParameters->TimeStep;
	pressureDomain = (zone->grid+gdi)->V;
	for(i=0;i<6;i++)
	{
		slabPressure[i] = (zone->grid+gdi)->SlabV[i];
	}
	ptrSize = (zone->grid+gdi)->GSize;
	
	NI=*(ptrSize+0);
	NJ=*(ptrSize+1);
	NK=*(ptrSize+2);
	
	flagrho = zone->Flowsol->flagRho;

//--------------------------------------------- Boundary cells -----------------------------------(1)--------

	k = 0; j = 0; i = 0;
	indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
	index[0]=k,index[1]=j,index[2]=i;
	
	rho = *(zone->Flowsol->Rho+indexP*flagrho);
	Vg = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
	V[0] = *(slabPressure[1]+(NJ-1)*k+j);
	V[1] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
	V[2] = *(slabPressure[3]+(NI-1)*k+i);
	V[3] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
	V[4] = *(slabPressure[5]+(NI-1)*j+i);
	V[5] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
	V[6] = *(slabPressure[1]+(NJ-1)*k+j+1);
	V[7] = 0.0;
	V[8] = *(slabPressure[1]+(NJ-1)*(k+1)+j);
	V[9] = 0.0;
	V[10] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1);
	V[11] = *(slabPressure[3]+(NI-1)*k+i+1);
	V[12] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1);
	V[13] = *(slabPressure[5]+(NI-1)*j+i+1);
	V[14] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i);
	V[15] = *(slabPressure[3]+(NI-1)*(k+1)+i);
	V[16] = *(slabPressure[5]+(NI-1)*(j+1)+i);
	V[17] = 0.0;
	
	MatrixMult(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
	
	resid = *((zone->grid+gdi)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS);

	*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/((delta_t/rho)*(*(diffdia+indexP)));
	
	residZone=residZone+resid*resid;
//--------------------------------------------------------------------------------------------------(2)---------
	k = 0; j = 0;
	
	for(i=1;i<NI-2;i++)
	{
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		Vg = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		V[0] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
		V[1] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
		V[2] = *(slabPressure[3]+(NI-1)*k+i);
		V[3] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
		V[4] = *(slabPressure[5]+(NI-1)*j+i);
		V[5] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
		V[6] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i-1);
		V[7] = *(slabPressure[3]+(NI-1)*k+i-1);
		V[8] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i-1);
		V[9] = *(slabPressure[5]+(NI-1)*j+i-1);
		V[10] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1);
		V[11] = *(slabPressure[3]+(NI-1)*k+i+1);
		V[12] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1);
		V[13] = *(slabPressure[5]+(NI-1)*j+i+1);
		V[14] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i);
		V[15] = *(slabPressure[3]+(NI-1)*(k+1)+i);
		V[16] = *(slabPressure[5]+(NI-1)*(j+1)+i);
		V[17] = 0.0;
		
		MatrixMult(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
		
		resid = *((zone->grid+gdi)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS);
	
		*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/((delta_t/rho)*(*(diffdia+indexP)));
		
		residZone=residZone+resid*resid;
	}
//------------------------------------------------------------------------------------------------------(3)----------
	k = 0; j = NJ-2;
	for(i=1;i<NI-2;i++)
	{
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		Vg = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		V[0] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
		V[1] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
		V[2] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
		V[3] = *(slabPressure[2]+(NI-1)*k+i);
		V[4] = *(slabPressure[5]+(NI-1)*j+i);
		V[5] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
		V[6] = *(slabPressure[2]+(NI-1)*k+i-1);
		V[7] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i-1);
		V[8] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i-1);
		V[9] = *(slabPressure[5]+(NI-1)*j+i-1);
		V[10] = *(slabPressure[2]+(NI-1)*k+i+1);
		V[11] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i+1);
		V[12] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1);
		V[13] = *(slabPressure[5]+(NI-1)*j+i+1);
		V[14] = *(slabPressure[2]+(NI-1)*(k+1)+i);
		V[15] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j-1)+i);
		V[16] = 0.0;
		V[17] = *(slabPressure[5]+(NI-1)*(j-1)+i);
		
		MatrixMult(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
		
		resid = *((zone->grid+gdi)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS);
	
		*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/((delta_t/rho)*(*(diffdia+indexP)));
		
		residZone=residZone+resid*resid;
	}
//----------------------------------------------------------------------------------------------------(4)-----------
	j = NJ-2;
	for(k=1;k<NK-2;k++)
	{
		for(i=1;i<NI-2;i++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			Vg = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			V[0] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			V[1] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			V[2] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			V[3] = *(slabPressure[2]+(NI-1)*k+i);
			V[4] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
			V[5] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
			V[6] = *(slabPressure[2]+(NI-1)*k+i-1);
			V[7] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i-1);
			V[8] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i-1);
			V[9] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i-1);
			V[10] = *(slabPressure[2]+(NI-1)*k+i+1);
			V[11] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i+1);
			V[12] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1);
			V[13] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i+1);
			V[14] = *(slabPressure[2]+(NI-1)*(k+1)+i);
			V[15] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j-1)+i);
			V[16] = *(slabPressure[2]+(NI-1)*(k-1)+i);
			V[17] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j-1)+i);
			
			MatrixMult(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
			
			resid = *((zone->grid+gdi)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS);
		
			*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/((delta_t/rho)*(*(diffdia+indexP)));
			
			residZone=residZone+resid*resid;
		}
	}
//------------------------------------------------------------------------------------------------------(5)----------
	k = 0; j = NJ-2; i = NI-2; 
	indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
	index[0]=k,index[1]=j,index[2]=i;
	
	rho = *(zone->Flowsol->Rho+indexP*flagrho);
	Vg = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
	V[0] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
	V[1] = *(slabPressure[0]+(NJ-1)*k+j);
	V[2] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
	V[3] = *(slabPressure[2]+(NI-1)*k+i);
	V[4] = *(slabPressure[5]+(NI-1)*j+i);
	V[5] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
	V[6] = *(slabPressure[2]+(NI-1)*k+i-1);
	V[7] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i-1);
	V[8] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i-1);
	V[9] = *(slabPressure[5]+(NI-1)*j+i-1);
	V[10] = 0.0;
	V[11] = *(slabPressure[0]+(NJ-1)*k+j-1);
	V[12] = *(slabPressure[0]+(NJ-1)*(k+1)+j);
	V[13] = 0.0;
	V[14] = *(slabPressure[2]+(NI-1)*(k+1)+i);
	V[15] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j-1)+i);
	V[16] = 0.0;
	V[17] = *(slabPressure[5]+(NI-1)*(j-1)+i);
	
	MatrixMult(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
	
	resid = *((zone->grid+gdi)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS);

	*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/((delta_t/rho)*(*(diffdia+indexP)));
	
	residZone=residZone+resid*resid;
//-----------------------------------------------------------------------------------------------------(6)----------
	j = NJ-2; i = NI-2;
	for(k=1;k<NK-2;k++)
	{
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		Vg = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		V[0] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
		V[1] = *(slabPressure[0]+(NJ-1)*k+j);
		V[2] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
		V[3] = *(slabPressure[2]+(NI-1)*k+i);
		V[4] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
		V[5] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
		V[6] = *(slabPressure[2]+(NI-1)*k+i-1);
		V[7] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i-1);
		V[8] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i-1);
		V[9] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i-1);
		V[10] = 0.0;
		V[11] = *(slabPressure[0]+(NJ-1)*k+j-1);
		V[12] = *(slabPressure[0]+(NJ-1)*(k+1)+j);
		V[13] = *(slabPressure[0]+(NJ-1)*(k-1)+j);
		V[14] = *(slabPressure[2]+(NI-1)*(k+1)+i);
		V[15] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j-1)+i);
		V[16] = *(slabPressure[2]+(NI-1)*(k-1)+i);
		V[17] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j-1)+i);
		
		MatrixMult(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
		
		resid = *((zone->grid+gdi)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS);
	
		*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/((delta_t/rho)*(*(diffdia+indexP)));
		
		residZone=residZone+resid*resid;
	}
//------------------------------------------------------------------------------------------------------(7)-----------
	k = NK-2; j = 0;
	for(i=1;i<NI-2;i++)
	{
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		Vg = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		V[0] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
		V[1] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
		V[2] = *(slabPressure[3]+(NI-1)*k+i);
		V[3] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
		V[4] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
		V[5] = *(slabPressure[4]+(NI-1)*j+i);
		V[6] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i-1);
		V[7] = *(slabPressure[3]+(NI-1)*k+i-1);
		V[8] = *(slabPressure[4]+(NI-1)*j+i-1);
		V[9] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i-1);
		V[10] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1);
		V[11] = *(slabPressure[3]+(NI-1)*k+i+1);
		V[12] = *(slabPressure[4]+(NI-1)*j+i+1);
		V[13] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i+1);
		V[14] = *(slabPressure[4]+(NI-1)*(j+1)+i);
		V[15] = 0.0;
		V[16] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j+1)+i);
		V[17] = *(slabPressure[3]+(NI-1)*(k-1)+i);
		
		MatrixMult(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
		
		resid = *((zone->grid+gdi)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS);
	
		*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/((delta_t/rho)*(*(diffdia+indexP)));
		
		residZone=residZone+resid*resid;
	}
//-----------------------------------------------------------------------------------------------------(8)--------------
	k = 0; j = 0; i = NI-2;
	indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
	index[0]=k,index[1]=j,index[2]=i;
	
	rho = *(zone->Flowsol->Rho+indexP*flagrho);
	Vg = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
	V[0] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
	V[1] = *(slabPressure[0]+(NJ-1)*k+j);
	V[2] = *(slabPressure[3]+(NI-1)*k+i);
	V[3] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
	V[4] = *(slabPressure[5]+(NI-1)*j+i);
	V[5] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
	V[6] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i-1);
	V[7] = *(slabPressure[3]+(NI-1)*k+i-1);
	V[8] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i-1);
	V[9] = *(slabPressure[5]+(NI-1)*j+i-1);
	V[10] = *(slabPressure[0]+(NJ-1)*k+j+1);
	V[11] = 0.0;
	V[12] = *(slabPressure[0]+(NJ-1)*(k+1)+j);
	V[13] = 0.0;
	V[14] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i);
	V[15] = *(slabPressure[3]+(NI-1)*(k+1)+i);
	V[16] = *(slabPressure[5]+(NI-1)*(j+1)+i);
	V[17] = 0.0;
	
	MatrixMult(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
	
	resid = *((zone->grid+gdi)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS);

	*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/((delta_t/rho)*(*(diffdia+indexP)));
	
	residZone=residZone+resid*resid;

//------------------------------------------------------------------------------------------------------(9)---------
	j = 0; i = 0;
	for(k=1;k<NK-2;k++)
	{
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		Vg = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		V[0] = *(slabPressure[1]+(NJ-1)*k+j);
		V[1] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
		V[2] = *(slabPressure[3]+(NI-1)*k+i);
		V[3] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
		V[4] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
		V[5] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
		V[6] = *(slabPressure[1]+(NJ-1)*k+j+1);
		V[7] = 0.0;
		V[8] = *(slabPressure[1]+(NJ-1)*(k+1)+j);
		V[9] = *(slabPressure[1]+(NJ-1)*(k-1)+j);
		V[10] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1);
		V[11] = *(slabPressure[3]+(NI-1)*k+i+1);
		V[12] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1);
		V[13] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i+1);
		V[14] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i);
		V[15] = *(slabPressure[3]+(NI-1)*(k+1)+i);
		V[16] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j+1)+i);
		V[17] = *(slabPressure[3]+(NI-1)*(k-1)+i);
		
		MatrixMult(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
		
		resid = *((zone->grid+gdi)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS);
	
		*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/((delta_t/rho)*(*(diffdia+indexP)));
		
		residZone=residZone+resid*resid;
	}
//--------------------------------------------------------------------------------------------------(10)------------
	j = 0; i = NI-2;
	for(k=1;k<NK-2;k++)
	{
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		Vg = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		V[0] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
		V[1] = *(slabPressure[0]+(NJ-1)*k+j);
		V[2] = *(slabPressure[3]+(NI-1)*k+i);
		V[3] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
		V[4] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
		V[5] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
		V[6] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i-1);
		V[7] = *(slabPressure[3]+(NI-1)*k+i-1);
		V[8] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i-1);
		V[9] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i-1);
		V[10] = *(slabPressure[0]+(NJ-1)*k+j+1);
		V[11] = 0.0;
		V[12] = *(slabPressure[0]+(NJ-1)*(k+1)+j);
		V[13] = *(slabPressure[0]+(NJ-1)*(k-1)+j);
		V[14] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i);
		V[15] = *(slabPressure[3]+(NI-1)*(k+1)+i);
		V[16] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j+1)+i);
		V[17] = *(slabPressure[3]+(NI-1)*(k-1)+i);
		
		MatrixMult(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
		
		resid = *((zone->grid+gdi)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS);
	
		*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/((delta_t/rho)*(*(diffdia+indexP)));
		
		residZone=residZone+resid*resid;
	}
//---------------------------------------------------------------------------------------------------(11)-----------
	j = NJ-2; i = 0;
	for(k=1;k<NK-2;k++)
	{
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		Vg = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		V[0] = *(slabPressure[1]+(NJ-1)*k+j);
		V[1] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
		V[2] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
		V[3] = *(slabPressure[2]+(NI-1)*k+i);
		V[4] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
		V[5] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
		V[6] = 0.0;
		V[7] = *(slabPressure[1]+(NJ-1)*k+j-1);
		V[8] = *(slabPressure[1]+(NJ-1)*(k+1)+j);
		V[9] = *(slabPressure[1]+(NJ-1)*(k-1)+j);
		V[10] = *(slabPressure[2]+(NI-1)*k+i+1);
		V[11] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i+1);
		V[12] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1);
		V[13] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i+1);
		V[14] = *(slabPressure[2]+(NI-1)*(k+1)+i);
		V[15] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j-1)+i);
		V[16] = *(slabPressure[2]+(NI-1)*(k-1)+i);
		V[17] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j-1)+i);
		
		MatrixMult(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
		
		resid = *((zone->grid+gdi)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS);
	
		*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/((delta_t/rho)*(*(diffdia+indexP)));
		
		residZone=residZone+resid*resid;
	}
//-----------------------------------------------------------------------------------------------------(12)---------
	j = 0; i = 0; k = NK-2;
	indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
	index[0]=k,index[1]=j,index[2]=i;
	
	rho = *(zone->Flowsol->Rho+indexP*flagrho);
	Vg = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
	V[0] = *(slabPressure[1]+(NJ-1)*k+j);
	V[1] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
	V[2] = *(slabPressure[3]+(NI-1)*k+i);
	V[3] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
	V[4] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
	V[5] = *(slabPressure[4]+(NI-1)*j+i);
	V[6] = *(slabPressure[1]+(NJ-1)*k+j+1);
	V[7] = 0.0;
	V[8] = 0.0;
	V[9] = *(slabPressure[1]+(NJ-1)*(k-1)+j);
	V[10] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1);
	V[11] = *(slabPressure[3]+(NI-1)*k+i+1);
	V[12] = *(slabPressure[4]+(NI-1)*j+i+1);
	V[13] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i+1);
	V[14] = *(slabPressure[4]+(NI-1)*(j+1)+i);
	V[15] = 0.0;
	V[16] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j+1)+i);
	V[17] = *(slabPressure[3]+(NI-1)*(k-1)+i);
	
	MatrixMult(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
	
	resid = *((zone->grid+gdi)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS);

	*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/((delta_t/rho)*(*(diffdia+indexP)));
	
	residZone=residZone+resid*resid;
//-----------------------------------------------------------------------------------------------------(13)---------
	k = 0; i = 0;
	for(j=1;j<NJ-2;j++)
	{
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		Vg = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		V[0] = *(slabPressure[1]+(NJ-1)*k+j);
		V[1] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
		V[2] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
		V[3] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
		V[4] = *(slabPressure[5]+(NI-1)*j+i);
		V[5] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
		V[6] = *(slabPressure[1]+(NJ-1)*k+j+1);
		V[7] = *(slabPressure[1]+(NJ-1)*k+j-1);
		V[8] = *(slabPressure[1]+(NJ-1)*(k+1)+j);
		V[9] = 0.0;
		V[10] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1);
		V[11] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i+1);
		V[12] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1);
		V[13] = *(slabPressure[5]+(NI-1)*j+i+1);
		V[14] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i);
		V[15] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j-1)+i);
		V[16] = *(slabPressure[5]+(NI-1)*(j+1)+i);
		V[17] = *(slabPressure[5]+(NI-1)*(j-1)+i);
		
		MatrixMult(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
		
		resid = *((zone->grid+gdi)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS);
	
		*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/((delta_t/rho)*(*(diffdia+indexP)));
		
		residZone=residZone+resid*resid;
	}
//-----------------------------------------------------------------------------------------------------(14)----------
	k = NK-2; i = 0;
	for(j=1;j<NJ-2;j++)
	{
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		Vg = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		V[0] = *(slabPressure[1]+(NJ-1)*k+j);
		V[1] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
		V[2] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
		V[3] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
		V[4] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
		V[5] = *(slabPressure[4]+(NI-1)*j+i);
		V[6] = *(slabPressure[1]+(NJ-1)*k+j+1);
		V[7] = *(slabPressure[1]+(NJ-1)*k+j-1);
		V[8] = 0.0;
		V[9] = *(slabPressure[1]+(NJ-1)*(k-1)+j);
		V[10] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1);
		V[11] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i+1);
		V[12] = *(slabPressure[4]+(NI-1)*j+i+1);
		V[13] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i+1);
		V[14] = *(slabPressure[4]+(NI-1)*(j+1)+i);
		V[15] = *(slabPressure[4]+(NI-1)*(j-1)+i);
		V[16] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j+1)+i);
		V[17] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j-1)+i);
		
		MatrixMult(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
		
		resid = *((zone->grid+gdi)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS);
	
		*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/((delta_t/rho)*(*(diffdia+indexP)));
		
		residZone=residZone+resid*resid;
	}
//---------------------------------------------------------------------------------------------------(15)------------
	k = 0; i = NI-2;
	for(j=1;j<NJ-2;j++)
	{
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		Vg = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		V[0] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
		V[1] = *(slabPressure[0]+(NJ-1)*k+j);
		V[2] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
		V[3] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
		V[4] = *(slabPressure[5]+(NI-1)*j+i);
		V[5] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
		V[6] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i-1);
		V[7] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i-1);
		V[8] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i-1);
		V[9] = *(slabPressure[5]+(NI-1)*j+i-1);
		V[10] = *(slabPressure[0]+(NJ-1)*k+j+1);
		V[11] = *(slabPressure[0]+(NJ-1)*k+j-1);
		V[12] = *(slabPressure[0]+(NJ-1)*(k+1)+j);
		V[13] = 0.0;
		V[14] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i);
		V[15] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j-1)+i);
		V[16] = *(slabPressure[5]+(NI-1)*(j+1)+i);
		V[17] = *(slabPressure[5]+(NI-1)*(j-1)+i);
		
		MatrixMult(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
		
		resid = *((zone->grid+gdi)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS);
	
		*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/((delta_t/rho)*(*(diffdia+indexP)));
		
		residZone=residZone+resid*resid;
	}

//-----------------------------------------------------------------------------------------------------(16)--------
	k = 0; i = 0; j = NJ-2;
	indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
	index[0]=k,index[1]=j,index[2]=i;
	
	rho = *(zone->Flowsol->Rho+indexP*flagrho);
	Vg = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
	V[0] = *(slabPressure[1]+(NJ-1)*k+j);
	V[1] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
	V[2] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
	V[3] = *(slabPressure[2]+(NI-1)*k+i);
	V[4] = *(slabPressure[5]+(NI-1)*j+i);
	V[5] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
	V[6] = 0.0;
	V[7] = *(slabPressure[1]+(NJ-1)*k+j-1);
	V[8] = *(slabPressure[1]+(NJ-1)*(k+1)+j);
	V[9] = 0.0;
	V[10] = *(slabPressure[2]+(NI-1)*k+i+1);
	V[11] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i+1);
	V[12] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1);
	V[13] = *(slabPressure[5]+(NI-1)*j+i+1);
	V[14] = *(slabPressure[2]+(NI-1)*(k+1)+i);
	V[15] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j-1)+i);
	V[16] = 0.0;
	V[17] = *(slabPressure[5]+(NI-1)*(j-1)+i);
	
	MatrixMult(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
	
	resid = *((zone->grid+gdi)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS);

	*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/((delta_t/rho)*(*(diffdia+indexP)));
	
	residZone=residZone+resid*resid;
//------------------------------------------------------------------------------------------------------(17)---------
	j = 0;
	for(k=1;k<NK-2;k++)
	{
		for(i=1;i<NI-2;i++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			Vg = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			V[0] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			V[1] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			V[2] = *(slabPressure[3]+(NI-1)*k+i);
			V[3] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			V[4] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
			V[5] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
			V[6] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i-1);
			V[7] = *(slabPressure[3]+(NI-1)*k+i-1);
			V[8] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i-1);
			V[9] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i-1);
			V[10] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1);
			V[11] = *(slabPressure[3]+(NI-1)*k+i+1);
			V[12] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1);
			V[13] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i+1);
			V[14] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i);
			V[15] = *(slabPressure[3]+(NI-1)*(k+1)+i);
			V[16] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j+1)+i);
			V[17] = *(slabPressure[3]+(NI-1)*(k-1)+i);
			
			MatrixMult(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
			
			resid = *((zone->grid+gdi)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS);
		
			*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/((delta_t/rho)*(*(diffdia+indexP)));
			
			residZone=residZone+resid*resid;
		}
	}
//-------------------------------------------------------------------------------------------------(18)-----------
	k = 0;
	for(j=1;j<NJ-2;j++)
	{
		for(i=1;i<NI-2;i++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			Vg = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			V[0] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			V[1] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			V[2] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			V[3] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			V[4] = *(slabPressure[5]+(NI-1)*j+i);
			V[5] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
			V[6] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i-1);
			V[7] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i-1);
			V[8] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i-1);
			V[9] = *(slabPressure[5]+(NI-1)*j+i-1);
			V[10] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1);
			V[11] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i+1);
			V[12] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1);
			V[13] = *(slabPressure[5]+(NI-1)*j+i+1);
			V[14] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i);
			V[15] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j-1)+i);
			V[16] = *(slabPressure[5]+(NI-1)*(j+1)+i);
			V[17] = *(slabPressure[5]+(NI-1)*(j-1)+i);
			
			MatrixMult(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
			
			resid = *((zone->grid+gdi)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS);
		
			*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/((delta_t/rho)*(*(diffdia+indexP)));
			
			residZone=residZone+resid*resid;
		}
	}
//--------------------------------------------------------------------------------------------------------(19)--------
	i=0;
	for(k=1;k<NK-2;k++)
	{
		for(j=1;j<NJ-2;j++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			Vg = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			V[0] = *(slabPressure[1]+(NJ-1)*k+j);
			V[1] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			V[2] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			V[3] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			V[4] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
			V[5] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
			V[6] = *(slabPressure[1]+(NJ-1)*k+j+1);
			V[7] = *(slabPressure[1]+(NJ-1)*k+j-1);
			V[8] = *(slabPressure[1]+(NJ-1)*(k+1)+j);
			V[9] = *(slabPressure[1]+(NJ-1)*(k-1)+j);
			V[10] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1);
			V[11] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i+1);
			V[12] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1);
			V[13] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i+1);
			V[14] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i);
			V[15] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j-1)+i);
			V[16] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j+1)+i);
			V[17] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j-1)+i);
			
			MatrixMult(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
			
			resid = *((zone->grid+gdi)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS);
		
			*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/((delta_t/rho)*(*(diffdia+indexP)));
			
			residZone=residZone+resid*resid;
		}
	}
//-----------------------------------------------------------------------------------------------------(20)---------
	k = NK-2; j = NJ-2; i = 0;
	indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
	index[0]=k,index[1]=j,index[2]=i;
	
	rho = *(zone->Flowsol->Rho+indexP*flagrho);
	Vg = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
	V[0] = *(slabPressure[1]+(NJ-1)*k+j);
	V[1] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
	V[2] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
	V[3] = *(slabPressure[2]+(NI-1)*k+i);
	V[4] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
	V[5] = *(slabPressure[4]+(NI-1)*j+i);
	V[6] = 0.0;
	V[7] = *(slabPressure[1]+(NJ-1)*k+j-1);
	V[8] = 0.0;
	V[9] = *(slabPressure[1]+(NJ-1)*(k-1)+j);
	V[10] = *(slabPressure[2]+(NI-1)*k+i+1);
	V[11] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i+1);
	V[12] = *(slabPressure[4]+(NI-1)*j+i+1);
	V[13] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i+1);
	V[14] = 0.0;
	V[15] = *(slabPressure[4]+(NI-1)*(j-1)+i);
	V[16] = *(slabPressure[2]+(NI-1)*(k-1)+i);
	V[17] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j-1)+i);
	
	MatrixMult(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
	
	resid = *((zone->grid+gdi)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS);

	*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/((delta_t/rho)*(*(diffdia+indexP)));
	
	residZone=residZone+resid*resid;
//-----------------------------------------------------------------------------------------------------(21)---------
	k = NK-2; j = NJ-2;
	for(i=1;i<NI-2;i++)
	{
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		Vg = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		V[0] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
		V[1] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
		V[2] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
		V[3] = *(slabPressure[2]+(NI-1)*k+i);
		V[4] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
		V[5] = *(slabPressure[4]+(NI-1)*j+i);
		V[6] = *(slabPressure[2]+(NI-1)*k+i-1);
		V[7] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i-1);
		V[8] = *(slabPressure[4]+(NI-1)*j+i-1);
		V[9] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i-1);
		V[10] = *(slabPressure[2]+(NI-1)*k+i+1);
		V[11] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i+1);
		V[12] = *(slabPressure[4]+(NI-1)*j+i+1);
		V[13] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i+1);
		V[14] = 0.0;
		V[15] = *(slabPressure[4]+(NI-1)*(j-1)+i);
		V[16] = *(slabPressure[2]+(NI-1)*(k-1)+i);
		V[17] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j-1)+i);
		
		MatrixMult(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
		
		resid = *((zone->grid+gdi)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS);
	
		*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/((delta_t/rho)*(*(diffdia+indexP)));
		
		residZone=residZone+resid*resid;
	}
//----------------------------------------------------------------------------------------------------(22)----------
	k = NK-2;
	for(j=1;j<NJ-2;j++)
	{
		for(i=1;i<NI-2;i++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			Vg = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			V[0] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			V[1] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			V[2] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			V[3] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			V[4] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
			V[5] = *(slabPressure[4]+(NI-1)*j+i);
			V[6] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i-1);
			V[7] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i-1);
			V[8] = *(slabPressure[4]+(NI-1)*j+i-1);
			V[9] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i-1);
			V[10] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1);
			V[11] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i+1);
			V[12] = *(slabPressure[4]+(NI-1)*j+i+1);
			V[13] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i+1);
			V[14] = *(slabPressure[4]+(NI-1)*(j+1)+i);
			V[15] = *(slabPressure[4]+(NI-1)*(j-1)+i);
			V[16] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j+1)+i);
			V[17] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j-1)+i);
			
			MatrixMult(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
			
			resid = *((zone->grid+gdi)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS);
		
			*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/((delta_t/rho)*(*(diffdia+indexP)));
			
			residZone=residZone+resid*resid;
		}
	}
//----------------------------------------------------------------------------------------------------(23)----------
	k = NK-2; j = 0; i = NI-2;
	indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
	index[0]=k,index[1]=j,index[2]=i;
	
	rho = *(zone->Flowsol->Rho+indexP*flagrho);
	Vg = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
	V[0] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
	V[1] = *(slabPressure[0]+(NJ-1)*k+j);
	V[2] = *(slabPressure[3]+(NI-1)*k+i);
	V[3] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
	V[4] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
	V[5] = *(slabPressure[4]+(NI-1)*j+i);
	V[6] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i-1);
	V[7] = *(slabPressure[3]+(NI-1)*k+i-1);
	V[8] = *(slabPressure[4]+(NI-1)*j+i-1);
	V[9] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i-1);
	V[10] = *(slabPressure[0]+(NJ-1)*k+j+1);
	V[11] = 0.0;
	V[12] = 0.0;
	V[13] = *(slabPressure[0]+(NJ-1)*(k-1)+j);
	V[14] = *(slabPressure[4]+(NI-1)*(j+1)+i);
	V[15] = 0.0;
	V[16] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j+1)+i);
	V[17] = *(slabPressure[3]+(NI-1)*(k-1)+i);
	
	MatrixMult(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
	
	resid = *((zone->grid+gdi)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS);

	*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/((delta_t/rho)*(*(diffdia+indexP)));
	
	residZone=residZone+resid*resid;
//-----------------------------------------------------------------------------------------------------(24)---------
	k = NK-2; i = NI-2;
	for(j=1;j<NJ-2;j++)
	{
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		Vg = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		V[0] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
		V[1] = *(slabPressure[0]+(NJ-1)*k+j);
		V[2] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
		V[3] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
		V[4] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
		V[5] = *(slabPressure[4]+(NI-1)*j+i);
		V[6] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i-1);
		V[7] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i-1);
		V[8] = *(slabPressure[4]+(NI-1)*j+i-1);
		V[9] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i-1);
		V[10] = *(slabPressure[0]+(NJ-1)*k+j+1);
		V[11] = *(slabPressure[0]+(NJ-1)*k+j-1);
		V[12] = 0.0;
		V[13] = *(slabPressure[0]+(NJ-1)*(k-1)+j);
		V[14] = *(slabPressure[4]+(NI-1)*(j+1)+i);
		V[15] = *(slabPressure[4]+(NI-1)*(j-1)+i);
		V[16] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j+1)+i);
		V[17] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j-1)+i);
		
		MatrixMult(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
		
		resid = *((zone->grid+gdi)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS);
	
		*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/((delta_t/rho)*(*(diffdia+indexP)));
		
		residZone=residZone+resid*resid;
	}
//-----------------------------------------------------------------------------------------------------(25)---------
	i = NI-2;
	for(k=1;k<NK-2;k++)
	{
		for(j=1;j<NJ-2;j++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			Vg = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			V[0] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			V[1] = *(slabPressure[0]+(NJ-1)*k+j);
			V[2] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			V[3] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			V[4] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
			V[5] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
			V[6] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i-1);
			V[7] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i-1);
			V[8] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i-1);
			V[9] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i-1);
			V[10] = *(slabPressure[0]+(NJ-1)*k+j+1);
			V[11] = *(slabPressure[0]+(NJ-1)*k+j-1);
			V[12] = *(slabPressure[0]+(NJ-1)*(k+1)+j);
			V[13] = *(slabPressure[0]+(NJ-1)*(k-1)+j);
			V[14] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i);
			V[15] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j-1)+i);
			V[16] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j+1)+i);
			V[17] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j-1)+i);
			
			MatrixMult(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
			
			resid = *((zone->grid+gdi)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS);
		
			*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/((delta_t/rho)*(*(diffdia+indexP)));
			
			residZone=residZone+resid*resid;
		}
	}
//---------------------------------------------------------------------------------------------------(26)-----------
	k = NK-2; j = NJ-2; i = NI-2;
	indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
	index[0]=k,index[1]=j,index[2]=i;
	
	rho = *(zone->Flowsol->Rho+indexP*flagrho);
	Vg = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
	V[0] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
	V[1] = *(slabPressure[0]+(NJ-1)*k+j);
	V[2] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
	V[3] = *(slabPressure[2]+(NI-1)*k+i);
	V[4] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
	V[5] = *(slabPressure[4]+(NI-1)*j+i);
	V[6] = *(slabPressure[2]+(NI-1)*k+i-1);
	V[7] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i-1);
	V[8] = *(slabPressure[4]+(NI-1)*j+i-1);
	V[9] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i-1);
	V[10] = 0.0;
	V[11] = *(slabPressure[0]+(NJ-1)*k+j-1);
	V[12] = 0.0;
	V[13] = *(slabPressure[0]+(NJ-1)*(k-1)+j);
	V[14] = 0.0;
	V[15] = *(slabPressure[4]+(NI-1)*(j-1)+i);
	V[16] = *(slabPressure[2]+(NI-1)*(k-1)+i);
	V[17] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j-1)+i);
	
	MatrixMult(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
	
	resid = *((zone->grid+gdi)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS);

	*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/((delta_t/rho)*(*(diffdia+indexP)));
	
	residZone=residZone+resid*resid;
//------------------------------------------- End of Boundary cells ----------------------------------------
	
//------------------------------------------- Interior cells -----------------------------------------(27)------
	for(k=1;k<NK-2;k++)
	for(j=1;j<NJ-2;j++)
	for(i=1;i<NI-2;i++)
	{
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
	
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		Vg = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		V[0] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
		V[1] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
		V[2] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
		V[3] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
		V[4] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
		V[5] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
		V[6] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i-1);
		V[7] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i-1);
		V[8] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i-1);
		V[9] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i-1);
		V[10] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1);
		V[11] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i+1);
		V[12] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1);
		V[13] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i+1);
		V[14] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i);
		V[15] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j-1)+i);
		V[16] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j+1)+i);
		V[17] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j-1)+i);
		
		MatrixMult(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
		
		resid = *((zone->grid+gdi)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS);
	
		*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/((delta_t/rho)*(*(diffdia+indexP)));
		
		residZone=residZone+resid*resid;
	}
//------------------------------------------- End of Interior cells ---------------------------------------------------
	residZone=sqrt(residZone/((NI-1)*(NJ-1)*(NK-1)));
	
	free(LHS);
	free(RHS);
	return(residZone);
}

