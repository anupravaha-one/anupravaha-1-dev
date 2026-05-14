/**
This function performs the matrix and vector multiplication, used in both Optigs and Bicgstab to do the same.
In case of BI-Cg it is used to find A*P and A*S while in the case of Optimized GS it is used to find A*delX.

If Var is 0-it finds A*delx and if it is 1 it calculates A*P and for var=2 it finds A*T.
Current function is for 3d cases

 **/

#include "APKCalcLHS3d.h"

double CalcLHS3d(Base *base,Zone *zone,int Var,int variableTag)
{

	double *delx,*Ax,*diffdia,*domainPointer[6];//,*Res_ini;//,*slabDelx[6];
	
	double Diag,Sum_Nb,*CoeffN[18],V[18],Vg,delta_t,*Rho,rho;
	
	int NI,NJ,NK,i,j,k,indexP,*ptrSize,flagrho,tag,nv;
	
	int J_Increment,K_Increment,ab,bc_type[6],*T_Indices[6];
	
	
	void *ptr;
	ZoneBoundary *bcPtr;
	Interface *facePtr;
	
	Zone *ZoneNeighbor;
	
	ptrSize = zone->Size ;
	
	NI=*(ptrSize+0);
	NJ=*(ptrSize+1);
	NK=*(ptrSize+2);
	
// 	TCell=(NI-1)*(NJ-1)*(NK-1);
	
	flagrho = zone->Flowsol->flagRho;
	
	Rho=zone->Flowsol->Rho;
	
	diffdia = (zone->grid)->DiffusionCoeff_P;
	
	if(variableTag==4)
		delta_t = base->BaseParameters->TimeStep;
	else
		delta_t=1.0;
	
	for(i=0;i<18;i++)
	{
		CoeffN[i] = (zone->grid)->DiffusionCoeff_nb[i];
	}
	
	if(Var==1)
	{
		delx = zone->Pi;
		
		Ax = zone->Vi;
		
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
	
	nv = find_nv(base,variableTag);
	
	for(i=0;i<6;i++)
	{
		//bc_type[i]=zone->bc_type[i];
		
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
			bc_type[i]=((bcPtr)->Data_Set+nv)->Bc_Type;//nv=2 for pressure
		}
	
	}
	
	K_Increment = (NI-1)*(NJ-1);
	J_Increment = (NI-1);
	//I_Increment = 1;
	
	
/// To find A*Delx
				
	indexP=0;

	for(k=0;k<NK-1;k++)
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
				V[0] = *(domainPointer[1]+*(T_Indices[1]+k*(NJ-1)+j));
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
				V[1] = *(domainPointer[0]+*(T_Indices[0]+k*(NJ-1)+j));
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
				V[2] = *(domainPointer[3]+*(T_Indices[3]+k*(NI-1)+i));
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
				V[3] = *(domainPointer[2]+*(T_Indices[2]+k*(NI-1)+i));
		}
		
		
		if(k!=0)
			V[4]=*(delx+indexP-K_Increment);
		else
		{
			if((bc_type[5]==0)||(bc_type[5]==1))
				V[4] = *(delx+indexP);//0.0;
			else if((bc_type[5]==3)||(bc_type[5]==2))
				V[4] = 0.0;
			else if(bc_type[5]==9)
				V[4] = *(domainPointer[5]+*(T_Indices[5]+j*(NI-1)+i));
		}	
	
		if(k!=(NK-2))
			V[5]=*(delx+indexP+K_Increment);
		else
		{
			if((bc_type[4]==0)||(bc_type[4]==1))
				V[5] = *(delx+indexP);//0.0;
			else if((bc_type[4]==3)||(bc_type[4]==2))
				V[5] = 0.0;
			else if(bc_type[4]==9)
				V[5] = *(domainPointer[4]+*(T_Indices[4]+j*(NI-1)+i));
		}
		
/*
		if(!(i==0 && j==(NJ-2)))
		{
			if(i==0)
			{
				if((bc_type[1]==0)||(bc_type[1]==1))
					V[6] = *(delx+indexP+J_Increment);
				else if((bc_type[1]==3)||(bc_type[1]==2))
					V[6] = 0.0;
				else if(bc_type[1]==9)
					V[6] = *(domainPointer[1]+*(T_Indices[1]+k*(NJ-1)+j+1));
			}	//V[6]=*(slabPressure[1]+(NJ-1)*k+(j+1));//from western slab
			else
			{
				if(j==(NJ-2))
				{
					if((bc_type[2]==0)||(bc_type[2]==1))
						V[6] = *(delx+indexP-1);//0.0;
					else if((bc_type[2]==3)||(bc_type[2]==2))
						V[6] = 0.0;
					else if(bc_type[2]==9)
						V[6] = *(domainPointer[2]+*(T_Indices[2]+k*(NI-1)+i-1));
				} //V[6]=*(slabPressure[2]+(NI-1)*k+(i-1));//from northern slab
				else
					V[6]=*(delx+indexP+J_Increment-I_Increment);
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
					V[7] = *(domainPointer[1]+*(T_Indices[1]+k*(NJ-1)+j-1));
			}//	V[7]=*(slabPressure[1]+(NJ-1)*k+(j-1));//from west slab
			else
			{
				if(j==0)
				{
					if((bc_type[3]==0)||(bc_type[3]==1))
						V[7] = *(delx+indexP-1);//0.0;
					else if((bc_type[3]==3)||(bc_type[3]==2))
						V[7] = 0.0;
					else if(bc_type[3]==9)
						V[7] = *(domainPointer[3]+*(T_Indices[3]+k*(NI-1)+i-1));
				}//	V[7]=*(slabPressure[3]+(NI-1)*k+(i-1));//from southern slab
				else
					V[7]=*(delx+indexP-J_Increment-I_Increment);
			}
		}
		else
			V[7]=0.0;
		
		if(!(i==0 && k==(NK-2)))
		{
			if(i==0)
			{
				if((bc_type[1]==0)||(bc_type[1]==1))
					V[8] = *(delx+indexP+K_Increment);//0.0;
				else if((bc_type[1]==3)||(bc_type[1]==2))
					V[8] = 0.0;
				else if(bc_type[1]==9)
					V[8] = *(domainPointer[1]+*(T_Indices[1]+(k+1)*(NJ-1)+j));
			} //	V[8]=*(slabPressure[1]+(NJ-1)*(k+1)+j);
			else
			{	
				if(k==(NK-2))
				{
					if((bc_type[4]==0)||(bc_type[4]==1))
						V[8] = *(delx+indexP-1);//0.0;
					else if((bc_type[4]==3)||(bc_type[4]==2))
						V[8] = 0.0;
					else if(bc_type[4]==9)
						V[8] = *(domainPointer[4]+*(T_Indices[4]+j*(NI-1)+i-1));
				} //	V[8]=*(slabPressure[4]+(NI-1)*j+(i-1));
				else
					V[8]=*(delx+indexP+K_Increment-I_Increment);
			}
			
		}
		else
			V[8]=0.0;
		
		if(!(i==0 && k==0))
		{
			if(i==0)
			{
				if((bc_type[1]==0)||(bc_type[1]==1))
					V[9] = *(delx+indexP-K_Increment);
				else if((bc_type[1]==3)||(bc_type[1]==2))
					V[9] = 0.0;
				else if(bc_type[1]==9)
					V[9] = *(domainPointer[1]+*(T_Indices[1]+(k-1)*(NJ-1)+j));
			}	//V[9]=*(slabPressure[1]+(NJ-1)*(k-1)+j);
			else
			{
				if(k==0)
				{
					if((bc_type[5]==0)||(bc_type[5]==1))
						V[9] = *(delx+indexP-1);//0.0;
					else if((bc_type[5]==3)||(bc_type[5]==2))
						V[9] = 0.0;
					else if(bc_type[5]==9)
						V[9] = *(domainPointer[5]+*(T_Indices[5]+j*(NI-1)+i-1));
				}	//V[9]=*(slabPressure[5]+(NI-1)*j+(i-1));
				else
					V[9]=*(delx+indexP-K_Increment-I_Increment);
			}
			
		}
		else
			V[9]=0.0;
			
		if(!(i==(NI-2) && j==(NJ-2)))
		{
			if(i==(NI-2))
			{
				if((bc_type[0]==0)||(bc_type[0]==1))
					V[10] = *(delx+indexP+J_Increment);//0.0;
				else if((bc_type[0]==3)||(bc_type[0]==2))
					V[10] = 0.0;
				else if(bc_type[0]==9)
					V[10] = *(domainPointer[0]+*(T_Indices[0]+k*(NJ-1)+j+1));
			}	//V[10]=*(slabPressure[0]+(NJ-1)*k+(j+1));
			else
			{
				if(j==(NJ-2))
				{
					if((bc_type[2]==0)||(bc_type[2]==1))
						V[10] = *(delx+indexP+1);//0.0;
					else if((bc_type[2]==3)||(bc_type[2]==2))
						V[10] = 0.0;
					else if(bc_type[2]==9)
						V[10] = *(domainPointer[2]+*(T_Indices[2]+k*(NI-1)+i+1));
				}//	V[10]=*(slabPressure[2]+(NI-1)*k+(i+1));
				else
					V[10]=*(delx+indexP+J_Increment+I_Increment);
			}
		}
		else
			V[10]=0.0;
		
		if(!(i==(NI-2) && j==0))
		{
			if(i==(NI-2))
			{
				if((bc_type[0]==0)||(bc_type[0]==1))
					V[11] = *(delx+indexP-J_Increment);//0.0;
				else if((bc_type[0]==3)||(bc_type[0]==2))
					V[11] = 0.0;
				else if(bc_type[0]==9)
					V[11] = *(domainPointer[0]+*(T_Indices[0]+k*(NJ-1)+j-1));
			}//	V[11]=*(slabPressure[0]+(NJ-1)*k+(j-1));
			else
			{	
				if(j==0)
				{
					if((bc_type[3]==0)||(bc_type[3]==1))
						V[11] = *(delx+indexP+1);//0.0;
					else if((bc_type[3]==3)||(bc_type[3]==2))
						V[11] = 0.0;
					else if(bc_type[3]==9)
						V[11] = *(domainPointer[3]+*(T_Indices[3]+k*(NI-1)+i+1));
				}//	V[11]=*(slabPressure[3]+(NI-1)*k+(i+1));
				else
					V[11]=*(delx+indexP-J_Increment+I_Increment);
			}
		}
		else
			V[11]=0.0;
		
		if(!(i==(NI-2) && k==(NK-2)))
		{
			if(i==(NI-2))
			{
				if((bc_type[0]==0)||(bc_type[0]==1))
					V[12] = *(delx+indexP+K_Increment);//0.0;
				else if((bc_type[0]==3)||(bc_type[0]==2))
					V[12] = 0.0;
				else if(bc_type[0]==9)
					V[12] = *(domainPointer[0]+*(T_Indices[0]+(k+1)*(NJ-1)+j));
			}//	V[12]=*(slabPressure[0]+(NJ-1)*(k+1)+j);
			else
			{
				if(k==(NK-2))
				{
					if((bc_type[4]==0)||(bc_type[4]==1))
						V[12] = *(delx+indexP+1);//0.0;
					else if((bc_type[4]==3)||(bc_type[4]==2))
						V[12] = 0.0;
					else if(bc_type[4]==9)
						V[12] = *(domainPointer[4]+*(T_Indices[4]+j*(NI-1)+i+1));
				} //	V[12]=*(slabPressure[4]+(NI-1)*j+(i+1));
				else
					V[12]=*(delx+indexP+K_Increment+I_Increment);
			}
		}
		else
			V[12]=0.0;	
		
		if(!(i==(NI-2) && k==0))
		{
			if(i==(NI-2))
			{
				if((bc_type[0]==0)||(bc_type[0]==1))
					V[13] = *(delx+indexP-K_Increment);//0.0;
				else if((bc_type[0]==3)||(bc_type[0]==2))
					V[13] = 0.0;
				else if(bc_type[0]==9)
					V[13] = *(domainPointer[0]+*(T_Indices[0]+(k-1)*(NJ-1)+j));
			}//	V[13]=*(slabPressure[0]+(NJ-1)*(k-1)+j);
			else
			{
				if(k==0)
				{
					if((bc_type[5]==0)||(bc_type[5]==1))
						V[13] = *(delx+indexP+1);//0.0;
					else if((bc_type[5]==3)||(bc_type[5]==2))
						V[13] = 0.0;
					else if(bc_type[5]==9)
						V[13] = *(domainPointer[5]+*(T_Indices[5]+j*(NI-1)+i+1));
				}//	V[13]=*(slabPressure[5]+(NI-1)*j+(i+1));
				else
					V[13]=*(delx+indexP-K_Increment+I_Increment);
			}
		}
		else
			V[13]=0.0;
		
		if(!(k==(NK-2) && j==(NJ-2)))
		{
			if(k==(NK-2))
			{
				if((bc_type[4]==0)||(bc_type[4]==1))
					V[14] = *(delx+indexP+J_Increment);//0.0;
				else if((bc_type[4]==3)||(bc_type[4]==2))
					V[14] = 0.0;
				else if(bc_type[4]==9)
					V[14] = *(domainPointer[4]+*(T_Indices[4]+(j+1)*(NI-1)+i));
			} //	V[14]=*(slabPressure[4]+(NI-1)*(j+1)+i);
			else
			{
				if(j==(NJ-2))
				{
					if((bc_type[2]==0)||(bc_type[2]==1))
						V[14] = *(delx+indexP+K_Increment);//0.0;
					else if((bc_type[2]==3)||(bc_type[2]==2))
						V[14] = 0.0;
					else if(bc_type[2]==9)
						V[14] = *(domainPointer[2]+*(T_Indices[2]+(k+1)*(NI-1)+i));
				}//	V[14]=*(slabPressure[2]+(NI-1)*(k+1)+i);
				else
					V[14]=*(delx+indexP+K_Increment+J_Increment);
			}
		}
		else
			V[14]=0.0;
		
		if(!(k==(NK-2) && j==0))
		{
			if(k==(NK-2))
			{
				if((bc_type[4]==0)||(bc_type[4]==1))
					V[15] = *(delx+indexP-J_Increment);//0.0;
				else if((bc_type[4]==3)||(bc_type[4]==2))
					V[15] = 0.0;
				else if(bc_type[4]==9)
					V[15] = *(domainPointer[4]+*(T_Indices[4]+(j-1)*(NI-1)+i));
			} //	V[15]=*(slabPressure[4]+(NI-1)*(j-1)+i);
			else
			{
				if(j==0)
				{
					if((bc_type[3]==0)||(bc_type[3]==1))
						V[15] = *(delx+indexP+K_Increment);//0.0;
					else if((bc_type[3]==3)||(bc_type[3]==2))
						V[15] = 0.0;
					else if(bc_type[3]==9)
						V[15] = *(domainPointer[3]+*(T_Indices[3]+(k+1)*(NI-1)+i));
				}//	V[15]=*(slabPressure[3]+(NI-1)*(k+1)+i);
				else
					V[15]=*(delx+indexP+K_Increment-J_Increment);
			}
		}
		else
			V[15]=0.0;
		
		if(!(k==0 && j==(NJ-2)))
		{
			if(k==0)
			{
				if((bc_type[5]==0)||(bc_type[5]==1))
					V[16] = *(delx+indexP+J_Increment);//0.0;
				else if((bc_type[5]==3)||(bc_type[5]==2))
					V[16] = 0.0;
				else if(bc_type[5]==9)
					V[16] = *(domainPointer[5]+*(T_Indices[5]+(j+1)*(NI-1)+i));
			}//	V[16]=*(slabPressure[5]+(NI-1)*(j+1)+i);
			else
			{
				if(j==(NJ-2))
				{
					if((bc_type[2]==0)||(bc_type[2]==1))
						V[16] = *(delx+indexP-K_Increment);//0.0;
					else if((bc_type[2]==3)||(bc_type[2]==2))
						V[16] = 0.0;
					else if(bc_type[2]==9)
						V[16] = *(domainPointer[2]+*(T_Indices[2]+(k-1)*(NI-1)+i));
				}//	V[16]=*(slabPressure[2]+(NI-1)*(k-1)+i);
				else
					V[16]=*(delx+indexP-K_Increment+J_Increment);
			}
		}
		else
			V[16]=0.0;
		
		if(!(k==0 && j==0))
		{
			if(k==0)
			{
				if((bc_type[5]==0)||(bc_type[5]==1))
					V[17] = *(delx+indexP-J_Increment);//0.0;
				else if((bc_type[5]==3)||(bc_type[5]==2))
					V[17] = 0.0;
				else if(bc_type[5]==9)
					V[17] = *(domainPointer[5]+*(T_Indices[5]+(j-1)*(NI-1)+i));
				}//	V[17]=*(slabPressure[5]+(NI-1)*(j-1)+i);
			else
			{
				if(j==0)
				{
					if((bc_type[3]==0)||(bc_type[3]==1))
						V[17] = *(delx+indexP-K_Increment);//0.0;
					else if((bc_type[3]==3)||(bc_type[3]==2))
						V[17] = 0.0;
					else if(bc_type[3]==9)
						V[17] = *(domainPointer[3]+*(T_Indices[3]+(k-1)*(NI-1)+i));
				}//	V[17]=*(slabPressure[3]+(NI-1)*(k-1)+i);
				else
					V[17]=*(delx+indexP-K_Increment-J_Increment);
			}
		}
		else
			V[17]=0.0;



*/		
		
		
		rho = *(Rho+indexP*flagrho);
		
		Vg = *(delx+indexP);
	
		Diag = *(diffdia+indexP) * Vg;
		
		Sum_Nb =  *(CoeffN[0 ]+indexP) *V[0 ] + *(CoeffN[1 ]+indexP) *V[1 ] + *(CoeffN[2 ]+indexP) *V[2 ]
			+ *(CoeffN[3 ]+indexP) *V[3 ] + *(CoeffN[4 ]+indexP) *V[4 ] + *(CoeffN[5 ]+indexP) *V[5 ];
//			+ *(CoeffN[6 ]+indexP) *V[6 ] + *(CoeffN[7 ]+indexP) *V[7 ] + *(CoeffN[8 ]+indexP) *V[8 ] 
//			+ *(CoeffN[9 ]+indexP) *V[9 ] + *(CoeffN[10]+indexP) *V[10] + *(CoeffN[11]+indexP) *V[11]
//			+ *(CoeffN[12]+indexP) *V[12] + *(CoeffN[13]+indexP) *V[13] + *(CoeffN[14]+indexP) *V[14]
//			+ *(CoeffN[15]+indexP) *V[15] + *(CoeffN[16]+indexP) *V[16] + *(CoeffN[17]+indexP)*V[17];
		
		
		*(Ax+indexP) = (delta_t/rho)*(Sum_Nb+Diag);
		
		indexP++;
	
	}
		
		
	return(0);
}


