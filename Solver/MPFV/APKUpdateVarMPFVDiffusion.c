#include "APKUpdateVarMPFVDiffusion.h"

/**
This function evaluate diffusion term for individual momentum equation by using overall momentum value. This is implemented for 2D/3D case having  variable diffusion coefficient.
*/

double UpdateVarMPFVDiffusion(Base *base,Zone *zone,int vartag)
{
	int NI,NJ,NK,index[3],*ptrSize,cntDir,indexP;
	double AX,*CoeffN[18],*diffdia;
	int k,j,i,flagrho,nop;
	double Vg,V[18],*LHS,*RHS,*temp;
	double *prevIterVar,*variable,*prevTstepVar,*slabVariable[6],*ptrVolume;

	FILE *fd;
	fd = fopen("diff.dat","w");
	ptrSize= zone->Size;
	ptrVolume= zone->ZGeom->Volume;

	nop = base->nphase;
	LHS = malloc(sizeof(double));
	RHS = malloc(sizeof(double));
	
	

	if(vartag==41||vartag==44||vartag==47||vartag==50||vartag==53)
	temp = zone->Flowsol->diffu;
	else if(vartag==42||vartag==45||vartag==48||vartag==51||vartag==54)
	temp = zone->Flowsol->diffv;
	else if(vartag==43||vartag==46||vartag==49||vartag==52||vartag==55)
	temp = zone->Flowsol->diffw;
	
	
	
	switch(vartag)
	{
		case 41:
		case 44:
		case 47:
		case 50:
		case 53:
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->SlabVx[cntDir];
			}
			diffdia = (zone->grid)->DiffusionCoeffPhi_P;
			for(i=0;i<18;i++)
			{
				CoeffN[i] = (zone->grid)->DiffusionCoeffPhi_nb[i];
			}
			variable = zone->Flowsol->Momentumx;
			prevTstepVar = zone->PrevFlowsol->Momentumx;
		break;
		case 42:
		case 45:
		case 48:
		case 51:
		case 54:
			
			
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->SlabVy[cntDir];
			}
			diffdia = (zone->grid)->DiffusionCoeffPhi_P;
			for(i=0;i<18;i++)
			{
				CoeffN[i] = (zone->grid)->DiffusionCoeffPhi_nb[i];
				
			}
			
			variable = zone->Flowsol->Momentumy;
			
			prevTstepVar = zone->PrevFlowsol->Momentumy;
		break;
		case 43:
		case 46:
		case 49:
		case 52:
		case 55:
			for(cntDir=0;cntDir<6;cntDir++)
			{
				slabVariable[cntDir] = zone->Flowsol->SlabVz[cntDir];
			}
			diffdia = (zone->grid)->DiffusionCoeffPhi_P;
			for(i=0;i<18;i++)
			{
				CoeffN[i] = (zone->grid)->DiffusionCoeffPhi_nb[i];
			}
			variable = zone->Flowsol->Momentumz;
			prevTstepVar = zone->PrevFlowsol->Momentumz;
		break;
		
		default:
			printf("Wrong Tag supplied into UpdateVarConjugate2d.c...\n");
			exit(1);
		break;
	}
	
	NI=*(ptrSize+0);
	NJ=*(ptrSize+1);
	NK=*(ptrSize+2);
	
	prevIterVar = zone->Flowsol->prevIterVariable;
	for(k=0;k<NK-1;k++)
	for(j=0;j<NJ-1;j++)
	for(i=0;i<NI-1;i++)
	{
		indexP=(NJ-1)*(NI-1)*k+(NI-1)*j+i;
		*(prevIterVar+indexP) = *(variable+indexP);
	}
	
	
	flagrho = zone->Flowsol->flagRho;
	
	
	
	for(k=0;k<NK-1;k++)
	for(j=0;j<NJ-1;j++)
	for(i=0;i<NI-1;i++)
	{
	
		if(i!=0)
			V[0]=*(variable+(NJ-1)*(NI-1)*k+(NI-1)*j+(i-1));
		else
			V[0] =*(slabVariable[1]+(NJ-1)*k+j);
			
		if(i!=(NI-2))
			V[1]=*(variable+(NJ-1)*(NI-1)*k+(NI-1)*j+(i+1));
		else
			V[1]= *(slabVariable[0]+(NJ-1)*k+j);
	
		if(j!=0)
			V[2]=*(variable+(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i);
		else
			V[2]=*(slabVariable[3]+(NI-1)*k+i);
	
		if(j!=(NJ-2))
			V[3]=*(variable+(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+i);
		else
			V[3]=*(slabVariable[2]+(NI-1)*k+i);
	
		if(k!=0)
			V[4]=*(variable+(NJ-1)*(NI-1)*(k-1)+(NI-1)*j+i);
		else
			V[4]=*(slabVariable[5]+(NI-1)*j+i);
	
		if(k!=(NK-2))
			V[5]=*(variable+(NJ-1)*(NI-1)*(k+1)+(NI-1)*j+i);
		else
			V[5]=*(slabVariable[4]+(NI-1)*j+i);
		
			
		if(!(i==0 && j==(NJ-2)))
		{
			if(i==0)
			V[6]=*(slabVariable[1]+(NJ-1)*k+(j+1));//from western slab
			else
			{
				if(j==(NJ-2))
			V[6]=*(slabVariable[2]+(NI-1)*k+(i-1));//from northern slab
				else
				V[6]=*(variable+(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+(i-1));
			}
		}
		else
			V[6]=0.0;
		
		if(!(i==0 && j==0))
		{
			if(i==0)
			V[7]=*(slabVariable[1]+(NJ-1)*k+(j-1));//from west slab
			else
			{
				if(j==0)
				V[7]=*(slabVariable[3]+(NI-1)*k+(i-1));//from southern slab
				else
				V[7]=*(variable+(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+(i-1));
			}
		}
		else
			V[7]=0.0;
		
		if(!(i==0 && k==(NK-2)))
		{
			if(i==0)
			V[8]=*(slabVariable[1]+(NJ-1)*(k+1)+j);
			else
			{	
				if(k==(NK-2))
				V[8]=*(slabVariable[4]+(NI-1)*j+(i-1));
				else
				V[8]=*(variable+(NJ-1)*(NI-1)*(k+1)+(NI-1)*j+(i-1));
			}
			
		}
		else
			V[8]=0.0;
		
		if(!(i==0 && k==0))
		{
			if(i==0)
			V[9]=*(slabVariable[1]+(NJ-1)*(k-1)+j);
			else
			{
				if(k==0)
				V[9]=*(slabVariable[5]+(NI-1)*j+(i-1));
				else
				V[9]=*(variable+(NJ-1)*(NI-1)*(k-1)+(NI-1)*j+(i-1));
			}
			
		}
		else
			V[9]=0.0;
			
		if(!(i==(NI-2) && j==(NJ-2)))
		{
			if(i==(NI-2))
			V[10]=*(slabVariable[0]+(NJ-1)*k+(j+1));
			else
			{
				if(j==(NJ-2))
				V[10]=*(slabVariable[2]+(NI-1)*k+(i+1));
				else
				V[10]=*(variable+(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+(i+1));
			}
		}
		else
			V[10]=0.0;
		
		if(!(i==(NI-2) && j==0))
		{
			if(i==(NI-2))
			V[11]=*(slabVariable[0]+(NJ-1)*k+(j-1));
			else
			{	
				if(j==0)
				V[11]=*(slabVariable[3]+(NI-1)*k+(i+1));
				else
				V[11]=*(variable+(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+(i+1));
			}
		}
		else
			V[11]=0.0;
		
		if(!(i==(NI-2) && k==(NK-2)))
		{
			if(i==(NI-2))
			{
				V[12]=*(slabVariable[0]+(NJ-1)*(k+1)+j);
				//printf("Slabvairable = %e\n",V[12]);
			}
			else
			{
				if(k==(NK-2))
				V[12]=*(slabVariable[4]+(NI-1)*j+(i+1));
				else
				V[12]=*(variable+(NJ-1)*(NI-1)*(k+1)+(NI-1)*j+(i+1));
			}
		}
		else
			V[12]=0.0;	
		
		if(!(i==(NI-2) && k==0))
		{
			if(i==(NI-2))
			V[13]=*(slabVariable[0]+(NJ-1)*(k-1)+j);
			else
			{
				if(k==0)
				V[13]=*(slabVariable[5]+(NI-1)*j+(i+1));
				else
				V[13]=*(variable+(NJ-1)*(NI-1)*(k-1)+(NI-1)*j+(i+1));
			}
		}
		else
			V[13]=0.0;
		
		if(!(k==(NK-2) && j==(NJ-2)))
		{
			if(k==(NK-2))
			V[14]=*(slabVariable[4]+(NI-1)*(j+1)+i);
			else
			{
				if(j==(NJ-2))
				V[14]=*(slabVariable[2]+(NI-1)*(k+1)+i);
				else
				V[14]=*(variable+(NJ-1)*(NI-1)*(k+1)+(NI-1)*(j+1)+i);
			}
		}
		else
			V[14]=0.0;
		
		if(!(k==(NK-2) && j==0))
		{
			if(k==(NK-2))
			V[15]=*(slabVariable[4]+(NI-1)*(j-1)+i);
			else
			{
				if(j==0)
				V[15]=*(slabVariable[3]+(NI-1)*(k+1)+i);
				else
				V[15]=*(variable+(NJ-1)*(NI-1)*(k+1)+(NI-1)*(j-1)+i);
			}
		}
		else
			V[15]=0.0;
		
		if(!(k==0 && j==(NJ-2)))
		{
			if(k==0)
			V[16]=*(slabVariable[5]+(NI-1)*(j+1)+i);
			else
			{
				if(j==(NJ-2))
				V[16]=*(slabVariable[2]+(NI-1)*(k-1)+i);
				else
				V[16]=*(variable+(NJ-1)*(NI-1)*(k-1)+(NI-1)*(j+1)+i);
			}
		}
		else
			V[16]=0.0;
		
		if(!(k==0 && j==0))
		{
			if(k==0)
			V[17]=*(slabVariable[5]+(NI-1)*(j-1)+i);
			else
			{
				if(j==0)
				V[17]=*(slabVariable[3]+(NI-1)*(k-1)+i);
				else
				V[17]=*(variable+(NJ-1)*(NI-1)*(k-1)+(NI-1)*(j-1)+i);
			}
		}
		else
			V[17]=0.0;
			
		index[0]=k,index[1]=j,index[2]=i;
		
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		
		
		
		Vg = *(variable+indexP);
		
		
		MatrixMult(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
		

		AX = -(1.0)*(*LHS+*RHS);
		
		*(temp+indexP) = AX;
// 		if(vartag==44)
// 			printf("\nindexP=%d\tvalue=%e",indexP,*(temp+indexP));
		
		fprintf(fd,"index=%d\tvartag=%d\tdiffusion=%e\n",indexP,vartag,*(temp+indexP));
	
	
			
	}
// 	printf("a");
//  printf("vartag=%d\tdiffusion=%e\n",vartag,)
	fclose(fd);
	return(0);
}


