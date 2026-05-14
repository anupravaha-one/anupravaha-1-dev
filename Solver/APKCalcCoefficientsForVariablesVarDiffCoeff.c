#include "APKCalcCoefficientsForVariablesVarDiffCoeff.h"

/**
This function calculates coeffients of diffusion term in the convection-diffusion equation for velocities, temperature
and scalars. This is used for the case of varying diffusion coefficients. It is used for 3D case. It is called only for
interior cells.
*/

double CalcCoefficientsForVariablesVarDiffCoeff(Zone *zone,int index[3],double vol_p,double u[18],int *size,double *beta,double *DiaN[18],int vartag)
{
	int NK,NJ,NI,k,j,i,comp,indexP,iFace;
//	double vol_e,vol_w,vol_s,vol_n,vol_t,vol_b,vol_ws,vol_wn,vol_wt,vol_wb,vol_es,vol_en,vol_et,vol_eb,vol_ts,vol_tn,vol_bs,vol_bn;
	double betaLocal[6][3],*diffusionCoeffdomain,*slabDiffusionCoeff[6];
	double coeffUp[12],diag,coeffWp[4],coeffEp[4],coeffNp[4],coeffSp[4],coeffTp[4],coeffBp[4];
	double coeffPWS,coeffPWN,coeffPWT,coeffPWB,coeffPEN,coeffPES,coeffPET,coeffPEB,coeffPTN,coeffPBN,coeffPTS,coeffPBS;
	double diffCoeff_p,diffCoeff_e,diffCoeff_w,diffCoeff_s,diffCoeff_n,diffCoeff_t,diffCoeff_b;//These are Cell center values
	double diffCoeff_Eface,diffCoeff_Wface,diffCoeff_Sface,diffCoeff_Nface,diffCoeff_Tface,diffCoeff_Bface;//These are face center values
	
	NI=*(size+0);
	NJ=*(size+1);
	NK=*(size+2);
	
	if((vartag>=31)&&(vartag<=55))
	{	
		diag=CalcCoefficientsForVariablesVarDiffCoeffMPFV(zone,index, vol_p, u,size,beta,DiaN,vartag);
	}
	else
	{
		diffusionCoeffdomain = MakeDiffusionCoeffDomain(zone,vartag);
		MakeDiffusionCoeffSlabPtrArray(zone,slabDiffusionCoeff,vartag);
		
		k=index[0],j=index[1],i=index[2];
		
		diffCoeff_p=*(diffusionCoeffdomain+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
		
		if(i!=0)
		{
			diffCoeff_w=*(diffusionCoeffdomain+(NJ-1)*(NI-1)*k+(NI-1)*j+(i-1));
		}
		else
		{
			diffCoeff_w =*(slabDiffusionCoeff[1]+(NJ-1)*k+j);
		}
		if(i!=(NI-2))
		{
			diffCoeff_e=*(diffusionCoeffdomain+(NJ-1)*(NI-1)*k+(NI-1)*j+(i+1));
		}
		else
		{
			diffCoeff_e=*(slabDiffusionCoeff[0]+(NJ-1)*k+j);
		}
		if(j!=0)
		{
			diffCoeff_s=*(diffusionCoeffdomain+(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i);
		}
		else
		{
			diffCoeff_s=*(slabDiffusionCoeff[3]+(NI-1)*k+i);
		}
		if(j!=(NJ-2))
		{
			diffCoeff_n=*(diffusionCoeffdomain+(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+i);
		}
		else
		{
			diffCoeff_n=*(slabDiffusionCoeff[2]+(NI-1)*k+i);
		}
		if(k!=0)
		{
			diffCoeff_b=*(diffusionCoeffdomain+(NJ-1)*(NI-1)*(k-1)+(NI-1)*j+i);
		}
		else
		{
			diffCoeff_b=*(slabDiffusionCoeff[5]+(NI-1)*j+i);
		}
		if(k!=(NK-2))
		{
			diffCoeff_t=*(diffusionCoeffdomain+(NJ-1)*(NI-1)*(k+1)+(NI-1)*j+i);
		}
		else
		{
			diffCoeff_t=*(slabDiffusionCoeff[4]+(NI-1)*j+i);
		}
		
	/*	vol_w = u[0];
		vol_e = u[1];
		vol_s = u[2];
		vol_n = u[3];
		vol_b = u[4];
		vol_t = u[5];
		vol_wn = u[6];
		vol_ws = u[7];
		vol_wt = u[8];
		vol_wb = u[9];
		vol_en = u[10];
		vol_es = u[11];
		vol_et = u[12];
		vol_eb = u[13];
		vol_tn = u[14];
		vol_ts = u[15];
		vol_bn = u[16];
		vol_bs = u[17];*/
		
		coeffPWN=vol_p/(u[6]+vol_p+u[3]+u[0]);
		coeffPWS=vol_p/(u[7]+vol_p+u[2]+u[0]);
		coeffPWT=vol_p/(u[8]+vol_p+u[5]+u[0]);
		coeffPWB=vol_p/(u[9]+vol_p+u[4]+u[0]);
		
		coeffPEN=vol_p/(u[10]+vol_p+u[3]+u[1]);
		coeffPES=vol_p/(u[11]+vol_p+u[2]+u[1]);
		coeffPET=vol_p/(u[12]+vol_p+u[5]+u[1]);
		coeffPEB=vol_p/(u[13]+vol_p+u[4]+u[1]);
		
		coeffPTN=vol_p/(u[14]+vol_p+u[3]+u[5]);
		coeffPTS=vol_p/(u[15]+vol_p+u[2]+u[5]);
		
		coeffPBN=vol_p/(u[16]+vol_p+u[3]+u[4]);
		coeffPBS=vol_p/(u[17]+vol_p+u[2]+u[4]);
		
		coeffUp[0]=u[7]/(u[7]+vol_p+u[2]+u[0])
		,coeffUp[1]=u[6]/(u[6]+vol_p+u[3]+u[0])
		,coeffUp[2]=u[9]/(u[9]+vol_p+u[4]+u[0])
		,coeffUp[3]=u[8]/(u[8]+vol_p+u[5]+u[0])
		,coeffUp[4]=u[11]/(u[11]+vol_p+u[2]+u[1])
		,coeffUp[5]=u[10]/(u[10]+vol_p+u[3]+u[1])
		,coeffUp[6]=u[13]/(u[13]+vol_p+u[4]+u[1])
		,coeffUp[7]=u[12]/(u[12]+vol_p+u[5]+u[1])
		,coeffUp[8]=u[17]/(u[17]+vol_p+u[4]+u[2])
		,coeffUp[9]=u[16]/(u[16]+vol_p+u[4]+u[3])
		,coeffUp[10]=u[15]/(u[15]+vol_p+u[2]+u[5])
		,coeffUp[11]=u[14]/(u[14]+vol_p+u[3]+u[5]);
		
		coeffWp[0]=u[2]/(u[7]+vol_p+u[2]+u[0])
		,coeffWp[1]=u[3]/(u[6]+vol_p+u[3]+u[0])
		,coeffWp[2]=u[4]/(u[9]+vol_p+u[4]+u[0])
		,coeffWp[3]=u[5]/(u[8]+vol_p+u[5]+u[0]);
		
		
		coeffEp[0]=u[2]/(u[11]+vol_p+u[2]+u[1])
		,coeffEp[1]=u[3]/(u[10]+vol_p+u[3]+u[1])
		,coeffEp[2]=u[4]/(u[13]+vol_p+u[4]+u[1])
		,coeffEp[3]=u[5]/(u[12]+vol_p+u[5]+u[1]); 
		
		coeffSp[1]=u[0]/(u[7]+vol_p+u[2]+u[0])
		,coeffSp[0]=u[1]/(u[11]+vol_p+u[2]+u[1])
		,coeffSp[3]=u[4]/(u[17]+vol_p+u[4]+u[2])
		,coeffSp[2]=u[5]/(u[15]+vol_p+u[2]+u[5]);
		
		
		coeffNp[1]=u[0]/(u[6]+vol_p+u[3]+u[0])
		,coeffNp[0]=u[1]/(u[10]+vol_p+u[3]+u[1])
		,coeffNp[3]=u[4]/(u[16]+vol_p+u[4]+u[3])
		,coeffNp[2]=u[5]/(u[14]+vol_p+u[3]+u[5]);
		
		coeffTp[1]=u[2]/(u[15]+vol_p+u[2]+u[5])
		,coeffTp[0]=u[3]/(u[14]+vol_p+u[3]+u[5])
		,coeffTp[3]=u[0]/(u[8]+vol_p+u[5]+u[0])
		,coeffTp[2]=u[1]/(u[12]+vol_p+u[5]+u[1]);
		
		
		coeffBp[1]=u[2]/(u[17]+vol_p+u[2]+u[4])
		,coeffBp[0]=u[3]/(u[16]+vol_p+u[3]+u[4])
		,coeffBp[3]=u[0]/(u[9]+vol_p+u[4]+u[0])
		,coeffBp[2]=u[1]/(u[13]+vol_p+u[4]+u[1]); 
		
	/*	if(((i==0)&&(j==0))||((i==0)&&(k==0))||((i==0)&&(j==(NJ-2)))||((i==0)&&(k==(NK-2)))||((i==(NI-2))&&(j==0))||((i==(NI-2))&&(k==0))||((i==(NI-2))&&(j==(NJ-2)))||((i==(NI-2))&&(k==(NK-2)))||((j==0)&&(k==0))||((j==0)&&(k==(NK-2)))||((j==(NJ-2))&&(k==0))||((j==(NJ-2))&&(k==(NK-2))))
		{
		
			coeffPWS=0
			,coeffPWN=0
			,coeffPWT=0
			,coeffPWB=0
			,coeffPEN=0
			,coeffPES=0
			,coeffPET=0
			,coeffPEB=0
			,coeffPTN=0
			,coeffPBN=0
			,coeffPTS=0
			,coeffPBS=0;
	
			coeffUp[0]=vol_w/(vol_w+vol_p)+vol_s/(vol_p+vol_s)-1.0
			,coeffUp[1]=vol_w/(vol_w+vol_p)+vol_n/(vol_p+vol_n)-1.0
			,coeffUp[2]=vol_w/(vol_w+vol_p)+vol_b/(vol_p+vol_b)-1.0
			,coeffUp[3]=vol_w/(vol_w+vol_p)+vol_t/(vol_p+vol_t)-1.0
			,coeffUp[4]=vol_e/(vol_p+vol_e)+vol_s/(vol_p+vol_s)-1.0
			,coeffUp[5]=vol_e/(vol_p+vol_e)+vol_n/(vol_p+vol_n)-1.0
			,coeffUp[6]=vol_e/(vol_p+vol_e)+vol_b/(vol_p+vol_b)-1.0
			,coeffUp[7]=vol_e/(vol_p+vol_e)+vol_t/(vol_p+vol_t)-1.0
			,coeffUp[8]=vol_b/(vol_b+vol_p)+vol_s/(vol_s+vol_p)-1.0
			,coeffUp[9]=vol_b/(vol_b+vol_p)+vol_n/(vol_n+vol_p)-1.0
			,coeffUp[10]=vol_t/(vol_t+vol_p)+vol_s/(vol_s+vol_p)-1.0
			,coeffUp[11]=vol_t/(vol_t+vol_p)+vol_n/(vol_n+vol_p)-1.0;
			
			
			coeffWp[0]=vol_p/(vol_p+vol_w)
			,coeffWp[1]=vol_p/(vol_p+vol_w)
			,coeffWp[2]=vol_p/(vol_p+vol_w)
			,coeffWp[3]=vol_p/(vol_p+vol_w);
		
		
			coeffEp[0]=vol_p/(vol_p+vol_e)
			,coeffEp[1]=vol_p/(vol_p+vol_e)
			,coeffEp[2]=vol_p/(vol_p+vol_e)
			,coeffEp[3]=vol_p/(vol_p+vol_e); 
		
			coeffSp[1]=vol_p/(vol_p+vol_s)
			,coeffSp[0]=vol_p/(vol_p+vol_s)
			,coeffSp[3]=vol_p/(vol_p+vol_s)
			,coeffSp[2]=vol_p/(vol_p+vol_s);
		
		
			coeffNp[1]=vol_p/(vol_p+vol_n)
			,coeffNp[0]=vol_p/(vol_p+vol_n)
			,coeffNp[3]=vol_p/(vol_p+vol_n)
			,coeffNp[2]=vol_p/(vol_p+vol_n);
			
			coeffTp[1]=vol_p/(vol_t+vol_p)
			,coeffTp[0]=vol_p/(vol_t+vol_p)
			,coeffTp[3]=vol_p/(vol_t+vol_p)
			,coeffTp[2]=vol_p/(vol_t+vol_p);
			
			
			coeffBp[1]=vol_p/(vol_b+vol_p)
			,coeffBp[0]=vol_p/(vol_b+vol_p)
			,coeffBp[3]=vol_p/(vol_b+vol_p)
			,coeffBp[2]=vol_p/(vol_b+vol_p); 
			
		}*/
		
		if((zone->Flowsol->Momentumx!=NULL)&&(vartag==4))
		{
			diffCoeff_Eface = 1/((diffCoeff_e*vol_p + diffCoeff_p*u[1])/(vol_p+u[1]));
			diffCoeff_Wface = 1/((diffCoeff_w*vol_p + diffCoeff_p*u[0])/(vol_p+u[0]));
			diffCoeff_Nface = 1/((diffCoeff_n*vol_p + diffCoeff_p*u[3])/(vol_p+u[3]));
			diffCoeff_Sface = 1/((diffCoeff_s*vol_p + diffCoeff_p*u[2])/(vol_p+u[2]));
			diffCoeff_Tface = 1/((diffCoeff_t*vol_p + diffCoeff_p*u[5])/(vol_p+u[5]));
			diffCoeff_Bface = 1/((diffCoeff_b*vol_p + diffCoeff_p*u[4])/(vol_p+u[4]));
		}
		else
		{
			diffCoeff_Eface = (diffCoeff_e*vol_p + diffCoeff_p*u[1])/(vol_p+u[1]);
			diffCoeff_Wface = (diffCoeff_w*vol_p + diffCoeff_p*u[0])/(vol_p+u[0]);
			diffCoeff_Nface = (diffCoeff_n*vol_p + diffCoeff_p*u[3])/(vol_p+u[3]);
			diffCoeff_Sface = (diffCoeff_s*vol_p + diffCoeff_p*u[2])/(vol_p+u[2]);
			diffCoeff_Tface = (diffCoeff_t*vol_p + diffCoeff_p*u[5])/(vol_p+u[5]);
			diffCoeff_Bface = (diffCoeff_b*vol_p + diffCoeff_p*u[4])/(vol_p+u[4]);
		}
		for(iFace=0;iFace<6;iFace++)
			for(comp=0;comp<3;comp++)
				betaLocal[iFace][comp]=*(beta+((NI-1)*(NJ-1)*k+(NI-1)*j+i)*6*3+3*iFace+comp);
				
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		
		diag=(diffCoeff_Wface*(betaLocal[0][0] +betaLocal[0][1]*(coeffUp[1]-coeffUp[0])+betaLocal[0][2]*(coeffUp[3]-coeffUp[2]))
			+diffCoeff_Eface*(-betaLocal[1][0] +betaLocal[1][1]*(coeffUp[5]-coeffUp[4])+betaLocal[1][2]*(coeffUp[7]-coeffUp[6]))
			+diffCoeff_Sface*(betaLocal[2][0]*(coeffUp[4]-coeffUp[0])+betaLocal[2][1] +betaLocal[2][2]*(coeffUp[10]-coeffUp[8]))
			+diffCoeff_Nface*(betaLocal[3][0]*(coeffUp[5]-coeffUp[1])-betaLocal[3][1] +betaLocal[3][2]*(coeffUp[11]-coeffUp[9]))
			+diffCoeff_Bface*(betaLocal[4][0]*(coeffUp[6]-coeffUp[2])+betaLocal[4][1]*(coeffUp[9]-coeffUp[8])+betaLocal[4][2])
			+diffCoeff_Tface*(betaLocal[5][0]*(coeffUp[7]-coeffUp[3])+betaLocal[5][1]*(coeffUp[11]-coeffUp[10])-betaLocal[5][2]));
			
		*(DiaN[0]+indexP)=-1.0*(diffCoeff_Wface*betaLocal[0][0]-diffCoeff_Wface*betaLocal[0][1]*(-coeffWp[0]+coeffWp[1])-diffCoeff_Wface*betaLocal[0][2]*(coeffWp[3]-coeffWp[2])
			+diffCoeff_Nface*betaLocal[3][0]*(coeffWp[1])+diffCoeff_Sface*betaLocal[2][0]*(coeffWp[0])+diffCoeff_Tface*betaLocal[5][0]*(coeffWp[3])+diffCoeff_Bface*betaLocal[4][0]*(coeffWp[2]));
			
		*(DiaN[1]+indexP)=-1.0*(-diffCoeff_Eface*betaLocal[1][0]-diffCoeff_Eface*betaLocal[1][1]*(-coeffEp[0]+coeffEp[1])-diffCoeff_Eface*betaLocal[1][2]*(coeffEp[3]-coeffEp[2])
			-diffCoeff_Nface*betaLocal[3][0]*(coeffEp[1])-diffCoeff_Sface*betaLocal[2][0]*(coeffEp[0])-diffCoeff_Tface*betaLocal[5][0]*(coeffEp[3])-diffCoeff_Bface*betaLocal[4][0]*(coeffEp[2]));
			
		*(DiaN[2]+indexP)=-1.0*(diffCoeff_Sface*betaLocal[2][1]+diffCoeff_Eface*betaLocal[1][1]*(coeffSp[0])+diffCoeff_Wface*betaLocal[0][1]*(coeffSp[1])-diffCoeff_Sface*betaLocal[2][0]*(coeffSp[0]-coeffSp[1])
			-diffCoeff_Sface*betaLocal[2][2]*(coeffSp[2]-coeffSp[3])+diffCoeff_Tface*betaLocal[5][1]*(coeffSp[2])+diffCoeff_Bface*betaLocal[4][1]*(coeffSp[3]));
			
		*(DiaN[3]+indexP)=-1.0*(-diffCoeff_Nface*betaLocal[3][1]-diffCoeff_Eface*betaLocal[1][1]*(coeffNp[0])-diffCoeff_Wface*betaLocal[0][1]*(coeffNp[1])-diffCoeff_Nface*betaLocal[3][0]*(coeffNp[0]-coeffNp[1])
			-diffCoeff_Nface*betaLocal[3][2]*(coeffNp[2]-coeffNp[3])-diffCoeff_Tface*betaLocal[5][1]*(coeffNp[2])-diffCoeff_Bface*betaLocal[4][1]*(coeffNp[3]));
			
		*(DiaN[4]+indexP)=-1.0*(diffCoeff_Bface*betaLocal[4][2]+diffCoeff_Wface*betaLocal[0][2]*(coeffBp[3])+diffCoeff_Eface*betaLocal[1][2]*(coeffBp[2])+diffCoeff_Nface*betaLocal[3][2]*(coeffBp[0])+diffCoeff_Sface*betaLocal[2][2]*(coeffBp[1])-diffCoeff_Bface*betaLocal[4][0]*(coeffBp[2]-coeffBp[3])-diffCoeff_Bface*betaLocal[4][1]*(coeffBp[0]-coeffBp[1]));
		
		*(DiaN[5]+indexP)=-1.0*(-diffCoeff_Tface*betaLocal[5][2]-diffCoeff_Wface*betaLocal[0][2]*(coeffTp[3])-diffCoeff_Eface*betaLocal[1][2]*(coeffTp[2])-diffCoeff_Nface*betaLocal[3][2]*(coeffTp[0])-diffCoeff_Sface*betaLocal[2][2]*(coeffTp[1])-diffCoeff_Tface*betaLocal[5][0]*(coeffTp[2]-coeffTp[3])-diffCoeff_Tface*betaLocal[5][1]*(coeffTp[0]-coeffTp[1]));
		
		*(DiaN[6]+indexP)=coeffPWN*(diffCoeff_Wface*betaLocal[0][1]-diffCoeff_Nface*betaLocal[3][0]);
		
		*(DiaN[7]+indexP)=-1.0*coeffPWS*(diffCoeff_Wface*betaLocal[0][1]+diffCoeff_Sface*betaLocal[2][0]);
		
		*(DiaN[8]+indexP)=coeffPWT*(diffCoeff_Wface*betaLocal[0][2]-diffCoeff_Tface*betaLocal[5][0]);
		
		*(DiaN[9]+indexP)=-1.0*coeffPWB*(diffCoeff_Wface*betaLocal[0][2]+diffCoeff_Bface*betaLocal[4][0]);
		
		*(DiaN[10]+indexP)=coeffPEN*(diffCoeff_Eface*betaLocal[1][1]+diffCoeff_Nface*betaLocal[3][0]);
		
		*(DiaN[11]+indexP)=coeffPES*(diffCoeff_Sface*betaLocal[2][0]-diffCoeff_Eface*betaLocal[1][1]);
		
		*(DiaN[12]+indexP)=coeffPET*(diffCoeff_Eface*betaLocal[1][2]+diffCoeff_Tface*betaLocal[5][0]);
		
		*(DiaN[13]+indexP)=-1.0*coeffPEB*(diffCoeff_Eface*betaLocal[1][2]-diffCoeff_Bface*betaLocal[4][0]);
		
		*(DiaN[14]+indexP)=coeffPTN*(diffCoeff_Nface*betaLocal[3][2]+diffCoeff_Tface*betaLocal[5][1]);
		
		*(DiaN[15]+indexP)=coeffPTS*(diffCoeff_Sface*betaLocal[2][2]-diffCoeff_Tface*betaLocal[5][1]);
		
		*(DiaN[16]+indexP)=-1.0*coeffPBN*(diffCoeff_Nface*betaLocal[3][2]-diffCoeff_Bface*betaLocal[4][1]);
		
		*(DiaN[17]+indexP)=-1.0*coeffPBS*(diffCoeff_Sface*betaLocal[2][2]+diffCoeff_Bface*betaLocal[4][1]);
	}
	
	return (diag);
}

