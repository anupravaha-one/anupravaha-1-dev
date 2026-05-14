/**
This function calculates the correction to the current approximation using residue/diagonal term for all the cells
this is used in OPti Gs iterations.
APKCalcCorrection3d.c handles 3d cases only.
 **/

#include "APKCalcCorrection3d.h"

//UnHash Those lines where V[6] to V[17] will be calculated in case of Non-Orthogonal 3d Grid. Hashed to save unnecessary multiplications.


double CalcCorrection3d(Base *base,Zone *zone)
{

	double resid,Diag,Sum_Nb,*diffdia,*domainPointer[6];
	
	double *CoeffN[18],V[18],Vg,*Rho,rho;//,*slabDelx[6];
	
	double *residue,*delx,rms=0.0,delta_t;
	
	int NI,NJ,NK,i,j,k,indexP,*ptrSize,flagrho,tag;
	
	int J_Increment,K_Increment,TCell,ab,bc_type[6],*T_Indices[6];//,k;
	///use findnv ...
	
	void *ptr;
	ZoneBoundary *bcPtr;
	Interface *facePtr;
	
	Zone *ZoneNeighbor;
	
	ptrSize = zone->Size ;
	
	NI=*(ptrSize+0);
	NJ=*(ptrSize+1);
	NK=*(ptrSize+2);
	
	TCell=(NI-1)*(NJ-1)*(NK-1);
	
	delx=zone->Correction;
	
	residue=zone->Residue;
	
	for(i=0;i<6;i++)
	{
		ptr = zone->ZONE_FACES[i];
		tag = zone->BC_IFACE_TAG[i];
		
		//if(tag=='i')
		if(tag==2)
		{
			facePtr=(Interface *)ptr;
			
			T_Indices[i]=facePtr->ZSlab->Target_Indices;
			
			ab=zone->ZONE_NEIGHBOUR_NUM[i];
		
			ZoneNeighbor=base->CgnsZone+ab;
		
			domainPointer[i]=ZoneNeighbor->Correction;
		
			bc_type[i]=9;
		}
		else if(tag==1)
		{
			bcPtr=(ZoneBoundary *)ptr;
			bc_type[i]=((bcPtr)->Data_Set+2)->Bc_Type;//2 is for pressure
		}
	}
	
	diffdia = (zone->grid)->DiffusionCoeff_P;
	
	for(i=0;i<18;i++)
	{
		CoeffN[i] = (zone->grid)->DiffusionCoeff_nb[i];
	}
	
	flagrho = zone->Flowsol->flagRho;
	
	Rho = zone->Flowsol->Rho;
	
	delta_t = base->BaseParameters->TimeStep;
	
	J_Increment = NI-1;
	K_Increment = (NI-1)*(NJ-1);
	
	rms=0.0;
		
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
		
		
		//resid = *(residue+indexP)-(delta_t/rho)*(Sum_Nb+Diag);

		resid = *(residue+indexP);
		
		rms+= resid*resid;
		
		*(delx+indexP) = resid/((delta_t/rho)*(*(diffdia+indexP)));
		
		indexP++;
	}
	
// 	rms=sqrt(rms/TCell);
// 	printf("Rms zone %e \n",rms);
	return(rms);
}


