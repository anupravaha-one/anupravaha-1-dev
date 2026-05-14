/**
This function calculates the correction to the current approximation using residue/diagonal term for all the cells
this is used in OPti Gs iterations.
APKCalcCorrection.c handles 2d cases only.
 **/



#include "APKCalcCorrection.h"


double CalcCorrection(Base *base,Zone *zone,int vartag)
{

	double resid,Diag,Sum_Nb,*diffdia,*domainPointer[6];
	
	double *CoeffN[18],V[18],Vg,*Rho,rho;//,*slabDelx[6];
	
	double *residue,*delx,rms=0.0,delta_t;
	
	int NI,NJ,NK,i,j,indexP,*ptrSize,flagrho,tag,nv;
	
	int J_Increment,TCell,ab,bc_type[6],*T_Indices[6];//,k;
	///use findnv ...
	
	void *ptr;
	ZoneBoundary *bcPtr;
	Interface *facePtr;
	
	Zone *ZoneNeighbor;
	
	ptrSize = zone->Size ;
	
	NI=*(ptrSize+0);
	NJ=*(ptrSize+1);
	NK=*(ptrSize+2);
	
	TCell=(NI-1)*(NJ-1);
	
	delx=zone->Correction;
	
	residue=zone->Residue;
	
	if(vartag==4)
	{
		delta_t = base->BaseParameters->TimeStep;
		nv=2;
	}
	else
	{	
		delta_t=1.0;
		nv = find_nv(base,vartag);
	}
	
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
			bc_type[i]=((bcPtr)->Data_Set+nv)->Bc_Type;//2 is for pressure
		}
	}
	
	diffdia = (zone->grid)->DiffusionCoeff_P;
	
	for(i=0;i<18;i++)
	{
		CoeffN[i] = (zone->grid)->DiffusionCoeff_nb[i];
	}
	
	flagrho = zone->Flowsol->flagRho;
	
	Rho = zone->Flowsol->Rho;
	
	J_Increment = NI-1;
	
	rms=0.0;
		
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
		
		Sum_Nb = *(CoeffN[0]+indexP) *V[0] + *(CoeffN[1]+indexP) *V[1] + *(CoeffN[2]+indexP) *V[2] + *(CoeffN[3]+indexP) *V[3] + *(CoeffN[6]+indexP) *V[6] + *(CoeffN[7]+indexP) *V[7] + *(CoeffN[10]+indexP) *V[10]+*(CoeffN[11]+indexP) *V[11];
		
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


