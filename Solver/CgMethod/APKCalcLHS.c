/**
This function performs the matrix and vector multiplication, used in both Optigs and Bicgstab to do the same.
Current function is for 2d cases
In case of BI-Cg it is used to find A*P and A*S while in the case of Optimized GS it is used to find A*delX.

If Var is 0-it finds A*delx and if it is 1 it calculates A*P and for var=2 it finds A*T.

 **/



#include "APKCalcLHS.h"


double CalcLHS(Base *base,Zone *zone,int Var,int variableTag)
{

	double Diag,Sum_Nb,*domainPointer[6],*Res_ini,sum=0.0;//,*slabDelx[6];
	
	double *delx,*Ax,*diffdia,*CoeffN[18],V[18],Vg,delta_t,*Rho,rho;
	
	int NI,NJ,NK,i,j,indexP,*ptrSize,flagrho,tag,nv;
	
	int J_Increment,TCell,ab,bc_type[6],*T_Indices[6];//,k;
	
	
	void *ptr;
	ZoneBoundary *bcPtr;
	Interface *facePtr;
	
	Zone *ZoneNeighbor;
	
	ptrSize = zone->Size ;
	
	NI=*(ptrSize+0);
	NJ=*(ptrSize+1);
	NK=*(ptrSize+2);
	
	TCell=(NI-1)*(NJ-1);
	
	if(Var==1)
	{
		delx = zone->Pi;
		
		Ax = zone->Vi;
		
		Res_ini = zone->Residue_init;
	}
	else if(Var==2)
	{
		delx = zone->Si;
		
		Ax = zone->Ti;
	}
	else if(Var==0)
	{
		delx=zone->Correction;
	
		Ax=zone->AP;
	}
	
	flagrho = zone->Flowsol->flagRho;
	
	Rho=zone->Flowsol->Rho;
	
	if(variableTag==4)
	{	
		delta_t = base->BaseParameters->TimeStep;
		nv=2;
	}
	else
	{
		delta_t=1.0;
		nv=find_nv(base,variableTag);
	}
	
	diffdia = (zone->grid)->DiffusionCoeff_P;
	
	for(i=0;i<18;i++)
	{
		CoeffN[i] = (zone->grid)->DiffusionCoeff_nb[i];
	}
	
	for(i=0;i<6;i++)
	{
		ptr = zone->ZONE_FACES[i];
		tag = zone->BC_IFACE_TAG[i];
		
		ab=zone->ZONE_NEIGHBOUR_NUM[i];
		
		ZoneNeighbor=base->CgnsZone+ab;
		
		if(tag==2)
		{
			facePtr=(Interface *)ptr;
			
			T_Indices[i]=facePtr->ZSlab->Target_Indices;
			
			bc_type[i]=9;
			
			if(Var==1)
			{
				domainPointer[i]=ZoneNeighbor->Pi;
			}
			else if(Var==2)
			{
				domainPointer[i]=ZoneNeighbor->Si;
			}
			else if(Var==0)
			{
				domainPointer[i]=ZoneNeighbor->Correction;
			}
			
			if((zone->flagSolidFluid == 0)&&(ZoneNeighbor->flagSolidFluid == 1))
			{
				bc_type[i]=0;
			}
			else if(zone->flagSolidFluid == 1)
			{
				bc_type[i]=3;
			}
		}
		else if(tag==1)
		{
			bcPtr=(ZoneBoundary *)ptr;
			
			bc_type[i]=((bcPtr)->Data_Set+nv)->Bc_Type;//nv=2 is for pressure
		}
// 		printf("Bc Type %d\n",bc_type[i]);
	}
	
	J_Increment = NI-1;
	
	
/// To find A*Delx
				
	indexP=0;

	for(j=0;j<NJ-1;j++)
	for(i=0;i<NI-1;i++)
	{
		if(i!=0)
			V[0]=*(delx+indexP-1);
		else
		{
			if((bc_type[1]==0)||(bc_type[1]==1))
				V[0] = *(delx+indexP);//0.0;
			else if((bc_type[1]==3)||(bc_type[1]==2))
				V[0] = 0.0;
			else if(bc_type[1]==9)
				V[0] = *(domainPointer[1]+*(T_Indices[1]+j));//*(slabDelx[1]+j);
		}
	
		if(i!=(NI-2))
			V[1]=*(delx+indexP+1);
		else
		{
			if((bc_type[0]==0)||(bc_type[0]==1))
				V[1] = *(delx+indexP);//0.0;
			else if((bc_type[0]==3)||(bc_type[0]==2))
				V[1] = 0.0;
			else if(bc_type[0]==9)
				V[1] = *(domainPointer[0]+*(T_Indices[0]+j));//*(slabDelx[0]+j);
		}
		
		if(j!=0)
			V[2]=*(delx+indexP-J_Increment);
		else
		{
			if((bc_type[3]==0)||(bc_type[3]==1))
				V[2] = *(delx+indexP);//0.0;
			else if((bc_type[3]==3)||(bc_type[3]==2))
				V[2] = 0.0;
			else if(bc_type[3]==9)
				V[2] = *(domainPointer[3]+*(T_Indices[3]+i));//*(slabDelx[3]+i);
		}
		
		if(j!=(NJ-2))
			V[3]=*(delx+indexP+J_Increment);
		else
		{
			if((bc_type[2]==0)||(bc_type[2]==1))
				V[3] = *(delx+indexP);//0.0;
			else if((bc_type[2]==3)||(bc_type[2]==2))
				V[3] = 0.0;
			else if(bc_type[2]==9)
				V[3] = *(domainPointer[2]+*(T_Indices[2]+i));//*(slabDelx[2]+i);
		}
	
		
		if(!(i==0 && j==(NJ-2)))
		{
			if(i==0)
			{
				if((bc_type[1]==0)||(bc_type[1]==1))
					V[6] = *(delx+indexP+J_Increment);
				else if((bc_type[1]==3)||(bc_type[1]==2))
					V[6] = 0.0;
				else if(bc_type[1]==9)
					V[6] = *(domainPointer[1]+*(T_Indices[1]+j+1));
			}	
			else
			{
				if(j==(NJ-2))
				{
					if((bc_type[2]==0)||(bc_type[2]==1))
						V[6] = *(delx+indexP-1);//0.0;
					else if((bc_type[2]==3)||(bc_type[2]==2))
						V[6] = 0.0;
					else if(bc_type[2]==9)
						V[6] = *(domainPointer[2]+*(T_Indices[2]+i-1));
				}	
				else
					V[6]=*(delx+indexP+J_Increment-1);
			}
		}
		else
			V[6]=0.0;
		
		
		
		if(!(i==0 && j==0))
		{
			if(i==0)
			{
				if((bc_type[1]==0)||(bc_type[1]==1))
					V[7] = *(delx+indexP-J_Increment);
				else if((bc_type[1]==3)||(bc_type[1]==2))
					V[7] = 0.0;
				else if(bc_type[1]==9)
					V[7] = *(domainPointer[1]+*(T_Indices[1]+j-1));
			}	
			else
			{
				if(j==0)
				{
					if((bc_type[3]==0)||(bc_type[3]==1))
						V[7] = *(delx+indexP-1);
					else if((bc_type[3]==3)||(bc_type[3]==2))
						V[7] = 0.0;
					else if(bc_type[3]==9)
						V[7] = *(domainPointer[3]+*(T_Indices[3]+i-1));
				}	
				else
					V[7]=*(delx+indexP-J_Increment-1);
			}
		}
		else
			V[7]=0.0;
		
		
		
		if(!(i==(NI-2) && j==(NJ-2)))
		{
			if(i==(NI-2))
			{
				if((bc_type[0]==0)||(bc_type[0]==1))
					V[10] = *(delx+indexP+J_Increment);
				else if((bc_type[0]==3)||(bc_type[0]==2))
					V[10] = 0.0;
				else if(bc_type[0]==9)
					V[10] = *(domainPointer[0]+*(T_Indices[0]+j+1));
			}
			else
			{
				if(j==(NJ-2))
				{
					if((bc_type[2]==0)||(bc_type[2]==1))
						V[10] = *(delx+indexP+1);
					else if((bc_type[2]==3)||(bc_type[2]==2))
						V[10] = 0.0;
					else if(bc_type[2]==9)
						V[10] = *(domainPointer[2]+*(T_Indices[2]+i+1));
				}
				else
					V[10]=*(delx+indexP+J_Increment+1);
			}
		}
		else
			V[10]=0.0;
		
		if(!(i==(NI-2) && j==0))
		{
			if(i==(NI-2))
			{
				if((bc_type[0]==0)||(bc_type[0]==1))
					V[11] = *(delx+indexP-J_Increment);
				else if((bc_type[0]==3)||(bc_type[0]==2))
					V[11] = 0.0;
				else if(bc_type[0]==9)
					V[11] = *(domainPointer[0]+*(T_Indices[0]+j-1));
			}
			else
			{	
				if(j==0)
				{
					if((bc_type[3]==0)||(bc_type[3]==1))
						V[11] = *(delx+indexP+1);
					else if((bc_type[3]==3)||(bc_type[3]==2))
						V[11] = 0.0;
					else if(bc_type[3]==9)
						V[11] = *(domainPointer[3]+*(T_Indices[3]+i+1));
				}	
				else
					V[11]=*(delx+indexP-J_Increment+1);
			}
		}
		else
			V[11]=0.0;
			
		rho = *(Rho+indexP*flagrho);
		
		Vg = *(delx+indexP);
	
		Diag = *(diffdia+indexP) * Vg;
		
		Sum_Nb = *(CoeffN[0]+indexP) *V[0] + *(CoeffN[1]+indexP) *V[1] + *(CoeffN[2]+indexP) *V[2] +*(CoeffN[3]+indexP) *V[3] + *(CoeffN[6]+indexP) *V[6] + *(CoeffN[7]+indexP) *V[7] + *(CoeffN[10]+indexP) *V[10]+*(CoeffN[11]+indexP) *V[11];
		
				
		*(Ax+indexP) = (delta_t/rho)*(Sum_Nb+Diag);
		
/*		
		//if(indexP==15)
		{	
			//printf("%d Ax %.18e %.18e %.18e\n",indexP,*(Ax+indexP),Sum_Nb,Diag);
				
			printf("%d Diffdia %.18e Vg %.18e\n",indexP,*(diffdia+indexP),Vg);
			for(ab=0;ab<4;ab++)
			{
				printf("Coeff %.18e Vg %.18e\n",*(CoeffN[ab]+indexP),V[ab]);
			}
			printf("6Coeff %.18e Vg %.18e\n",*(CoeffN[6]+indexP),V[6]);
			printf("7Coeff %.18e Vg %.18e\n",*(CoeffN[7]+indexP),V[7]);
			printf("10Coeff %.18e Vg %.18e\n",*(CoeffN[10]+indexP),V[10]);
			printf("11Coeff %.18e Vg %.18e\n",*(CoeffN[11]+indexP),V[11]);
			//exit(0);
		}
*/		
// 		printf("%d Ax %e\n",indexP,*(Ax+indexP));
		
// 		if(Var==1)
// 		{
// 			sum += *(Ax+indexP) * *(Res_ini+indexP);
// 		}
		
		
		indexP++;
	}
		
		
	return(sum);
}


