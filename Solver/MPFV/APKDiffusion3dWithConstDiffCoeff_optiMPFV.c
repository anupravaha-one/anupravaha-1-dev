#include "APKDiffusion3dWithConstDiffCoeff_optiMPFV.h"
/**
The function Diffusion3dWithConstDiffCoeff() is used to compute the diffusion flux of the form (Del.(Gamma Del(phi)))
for 3 dimensional domain, where the diffusion coefficient Gamma is a constant. This function gets all the neighbouring cell
pointers either from the domain or from the zone slabs, depending upon the cell index, and computes the diffusion flux
on all the six faces of the cell. Volumetric interpolation of the diffusion coefficient is not done for computing the
diffusion flux through a face. See Eqn(2.22) to Eqn(2.31) in Project Report No 1, Dec 2005, for further details.
*/
double Diffusion3dWithConstDiffCoeff_optiMPFV(int index[3],Zone *zone,double *domain,double *volumes,int *size,double *slabVariable[6],double *slabVolume[6],double *beta,double DiffFlux[6])
{
	
	int NK,NJ,NI,k,j,i,iFace,comp,INDEX;
	double u_p,u_e,u_w,u_s,u_n,u_t,u_b,vol_p,vol_e,vol_w,vol_s,vol_n,vol_t,vol_b,u_ws,u_wn,u_wt,u_wb,u_es,
	u_en,u_et,u_eb,vol_ws,vol_wn,vol_wt,vol_wb,vol_es,vol_en,vol_et,vol_eb,u_ts,u_tn,u_bs,u_bn,vol_ts,
	vol_tn,vol_bs,vol_bn;
	double u_ws_e,u_wn_e,u_wb_e,u_wt_e,u_es_e,u_en_e,u_eb_e,u_et_e,u_bs_e,u_bn_e,u_ts_e,u_tn_e;
	double betaLocal[6][3];
	double diagonal=1.0;
	int K_Increment,J_Increment,I_Increment,IMAX,JMAX,KMAX;
	double *slabRho[6],*domainRho,r_p,r_w,r_s,r_n,r_t,r_b,r_e,re,rw,rt,rb,rs,rn;
	NI=*(size+0);
	NJ=*(size+1);
	NK=*(size+2);
	
	
	domainRho=MakeDomainArray(zone,710);/*Density*/
	MakeSlabPtrArray(zone,slabRho,710);
	
	k=index[0],j=index[1],i=index[2];
	//diffCoeff=*(diffusionCoeffdomain);
	
	INDEX = (NJ-1)*(NI-1)*k+(NI-1)*j+i;
	u_p=*(domain+INDEX);
	vol_p=*(volumes+INDEX);
	r_p = *(domainRho+INDEX);
	
	K_Increment = (NI-1)*(NJ-1);
	J_Increment = (NI-1);
	I_Increment = 1;
	IMAX = NI-2;
	JMAX = NJ-2;
	KMAX = NK-2;
	
	
	if(i!=0)
	{
		u_w=*(domain+INDEX-1);
		vol_w=*(volumes+INDEX-1);
		r_w=*(domainRho+INDEX-1);
	}
	else
	{
		u_w =*(slabVariable[1]+(NJ-1)*k+j);
		vol_w =*(slabVolume[1]+(NJ-1)*k+j);
		r_w =*(slabRho[1]+(NJ-1)*k+j);
	}
	if(i!=(NI-2))
	{
		u_e=*(domain+INDEX+1);
		vol_e=*(volumes+INDEX+1);
		r_e=*(domainRho+INDEX+1);
	}
	else
	{
		u_e= *(slabVariable[0]+(NJ-1)*k+j);
		vol_e=*(slabVolume[0]+(NJ-1)*k+j);
		r_e= *(slabRho[0]+(NJ-1)*k+j);
	}
	if(j!=0)
	{
		u_s=*(domain+INDEX-J_Increment);
		vol_s=*(volumes+INDEX-J_Increment);
		r_s=*(domainRho+INDEX-J_Increment);
	}
	else
	{
		u_s=*(slabVariable[3]+(NI-1)*k+i);
		vol_s=*(slabVolume[3]+(NI-1)*k+i);
		r_s=*(slabRho[3]+(NI-1)*k+i);
	}
	if(j!=(NJ-2))
	{
		u_n=*(domain+INDEX+J_Increment);
		vol_n=*(volumes+INDEX+J_Increment);
		r_n=*(domainRho+INDEX+J_Increment);
	}
	else
	{
		u_n=*(slabVariable[2]+(NI-1)*k+i);
		vol_n=*(slabVolume[2]+(NI-1)*k+i);
		r_n=*(slabRho[2]+(NI-1)*k+i);
	}
	if(k!=0)
	{
		u_b=*(domain+INDEX-K_Increment);
		vol_b=*(volumes+INDEX-K_Increment);
		r_b=*(domainRho+INDEX-K_Increment);
	}
	else
	{
		u_b=*(slabVariable[5]+(NI-1)*j+i);
		vol_b=*(slabVolume[5]+(NI-1)*j+i);
		r_b=*(slabRho[5]+(NI-1)*j+i);
	}
	if(k!=(NK-2))
	{
		u_t=*(domain+INDEX+K_Increment);
		vol_t=*(volumes+INDEX+K_Increment);
		r_t=*(domainRho+INDEX+K_Increment);
	}
	else
	{
		u_t=*(slabVariable[4]+(NI-1)*j+i);
		vol_t=*(slabVolume[4]+(NI-1)*j+i);
		r_t=*(slabRho[4]+(NI-1)*j+i);
	}
	
	if(!(i==0 && j==0))
	{
		if(i==0)
		{
			u_ws=*(slabVariable[1]+(NJ-1)*k+(j-1));/*from west slab*/
			vol_ws=*(slabVolume[1]+(NJ-1)*k+(j-1));
		}
		else
		{
			if(j==0)
			{
				u_ws=*(slabVariable[3]+(NI-1)*k+(i-1));/*from southern slab*/
				vol_ws=*(slabVolume[3]+(NI-1)*k+(i-1));
			}
			else
			{
				u_ws=*(domain+INDEX-J_Increment-1);
				vol_ws=*(volumes+INDEX-J_Increment-1);
			}
		}
	}
		
	if(!(i==0 && j==(NJ-2)))
	{
		if(i==0)
		{
			u_wn=*(slabVariable[1]+(NJ-1)*k+(j+1));/*from western slab*/
			vol_wn=*(slabVolume[1]+(NJ-1)*k+(j+1));
		}
		else
		{
			if(j==(NJ-2))
			{
				u_wn=*(slabVariable[2]+(NI-1)*k+(i-1));/*from northern slab*/
				vol_wn=*(slabVolume[2]+(NI-1)*k+(i-1));
			}
			else
			{
				u_wn=*(domain+INDEX+J_Increment-1);
				vol_wn=*(volumes+INDEX+J_Increment-1);
			}
		}
	}
	if(!(i==0 && k==0))
	{
		if(i==0)
		{
			u_wb=*(slabVariable[1]+(NJ-1)*(k-1)+j);
			vol_wb=*(slabVolume[1]+(NJ-1)*(k-1)+j);
		}
		else
		{
			if(k==0)
			{
				u_wb=*(slabVariable[5]+(NI-1)*j+(i-1));
				vol_wb=*(slabVolume[5]+(NI-1)*j+(i-1));
			}
			else
			{
				u_wb=*(domain+INDEX-K_Increment-1);
				vol_wb=*(volumes+INDEX-K_Increment-1);
			}
		}
		
	}
	
	if(!(i==0 && k==(NK-2)))
	{
		if(i==0)
		{
			u_wt=*(slabVariable[1]+(NJ-1)*(k+1)+j);
			vol_wt=*(slabVolume[1]+(NJ-1)*(k+1)+j);
		}
		else
		{
			if(k==(NK-2))
			{
				u_wt=*(slabVariable[4]+(NI-1)*j+(i-1));
				vol_wt=*(slabVolume[4]+(NI-1)*j+(i-1));
			}
			else
			{
				u_wt=*(domain+INDEX+K_Increment-1);
				vol_wt=*(volumes+INDEX+K_Increment-1);
			}
		}
		
	}
	if(!(i==(NI-2) && j==0))
	{
		if(i==(NI-2))
		{
			u_es=*(slabVariable[0]+(NJ-1)*k+(j-1));
			vol_es=*(slabVolume[0]+(NJ-1)*k+(j-1));
		}
		else
		{
			if(j==0)
			{
				u_es=*(slabVariable[3]+(NI-1)*k+(i+1));
				vol_es=*(slabVolume[3]+(NI-1)*k+(i+1));
			}
			else
			{
				u_es=*(domain+INDEX-J_Increment+1);
				vol_es=*(volumes+INDEX-J_Increment+1);
			}
		}
		
	}
		
	if(!(i==(NI-2) && j==(NJ-2)))
	{
		if(i==(NI-2))
		{
			u_en=*(slabVariable[0]+(NJ-1)*k+(j+1));
			vol_en=*(slabVolume[0]+(NJ-1)*k+(j+1));
		}
		else
		{
			if(j==(NJ-2))
			{
				u_en=*(slabVariable[2]+(NI-1)*k+(i+1));
				vol_en=*(slabVolume[2]+(NI-1)*k+(i+1));
			}
			else
			{
				u_en=*(domain+INDEX+J_Increment+1);
				vol_en=*(volumes+INDEX+J_Increment+1);
			}
		}
	}
	
	if(!(i==(NI-2) && k==0))
	{
		if(i==(NI-2))
		{
			u_eb=*(slabVariable[0]+(NJ-1)*(k-1)+j);
			vol_eb=*(slabVolume[0]+(NJ-1)*(k-1)+j);
		}
		else
		{
			if(k==0)
			{
				u_eb=*(slabVariable[5]+(NI-1)*j+(i+1));
				vol_eb=*(slabVolume[5]+(NI-1)*j+(i+1));
			}
			else
			{
				u_eb=*(domain+INDEX-K_Increment+1);
				vol_eb=*(volumes+INDEX-K_Increment+1);
			}
		}
	}
	
	if(!(i==(NI-2) && k==(NK-2)))
	{
		if(i==(NI-2))
		{
			u_et=*(slabVariable[0]+(NJ-1)*(k+1)+j);
			vol_et= *(slabVolume[0]+(NJ-1)*(k+1)+j);
		}
		else
		{
			if(k==(NK-2))
			{
				u_et=*(slabVariable[4]+(NI-1)*j+(i+1));
				vol_et=*(slabVolume[4]+(NI-1)*j+(i+1));
			}
			else
			{
				u_et=*(domain+INDEX+K_Increment+1);
				vol_et=*(volumes+INDEX+K_Increment+1);
			}
		}
	}	
	if(!(k==0 && j==0))
	{
		if(k==0)
		{
			u_bs=*(slabVariable[5]+(NI-1)*(j-1)+i);
			vol_bs=*(slabVolume[5]+(NI-1)*(j-1)+i);
		}
		else
		{
			if(j==0)
			{
				u_bs=*(slabVariable[3]+(NI-1)*(k-1)+i);
				vol_bs=*(slabVolume[3]+(NI-1)*(k-1)+i);
			}
			else
			{
				u_bs=*(domain+INDEX-J_Increment-K_Increment);
				vol_bs=*(volumes+INDEX-J_Increment-K_Increment);
			}
		}
		
	}
	if(!(k==0 && j==(NJ-2)))
	{
		if(k==0 )
		{
			u_bn=*(slabVariable[5]+(NI-1)*(j+1)+i);
			vol_bn=*(slabVolume[5]+(NI-1)*(j+1)+i);
		}
		else
		{
			if(j==(NJ-2))
			{
				u_bn=*(slabVariable[2]+(NI-1)*(k-1)+i);
				vol_bn=*(slabVolume[2]+(NI-1)*(k-1)+i);
			}
			else
			{
				u_bn=*(domain+INDEX+J_Increment-K_Increment);
				vol_bn=*(volumes+INDEX+J_Increment-K_Increment);
			}
		}
		
	}
	if(!(k==(NK-2) && j==0))
	{
		if(k==(NK-2))
		{
			u_ts=*(slabVariable[4]+(NI-1)*(j-1)+i);
			vol_ts=*(slabVolume[4]+(NI-1)*(j-1)+i);
		}
		else
		{
			if(j==0)
			{
				u_ts=*(slabVariable[3]+(NI-1)*(k+1)+i);
				vol_ts=*(slabVolume[3]+(NI-1)*(k+1)+i);
			}
			else
			{
				u_ts=*(domain+INDEX-J_Increment+K_Increment);
				vol_ts=*(volumes+INDEX-J_Increment+K_Increment);
			}
		}
	}
	
	if(!(k==(NK-2) && j==(NJ-2)))
	{
		if(k==(NK-2))
		{
			u_tn=*(slabVariable[4]+(NI-1)*(j+1)+i);
			vol_tn=*(slabVolume[4]+(NI-1)*(j+1)+i);
		}
		else
		{
			if(j==(NJ-2))
			{
				u_tn=*(slabVariable[2]+(NI-1)*(k+1)+i);
				vol_tn=*(slabVolume[2]+(NI-1)*(k+1)+i);
			}
			else
			{
				u_tn=*(domain+INDEX+J_Increment+K_Increment);
				vol_tn=*(volumes+INDEX+J_Increment+K_Increment);
			}
		}
	}


	if(((i==0)&&(j==0))||((i==0)&&(k==0))||((i==0)&&(j==(NJ-2)))||((i==0)&&(k==(NK-2)))||((i==(NI-2))&&(j==0))||((i==(NI-2))&&(k==0))||((i==(NI-2))&&(j==(NJ-2)))||((i==(NI-2))&&(k==(NK-2)))||((j==0)&&(k==0))||((j==0)&&(k==(NK-2)))||((j==(NJ-2))&&(k==0))||((j==(NJ-2))&&(k==(NK-2))))
	{
		u_ws_e=(u_p*vol_w+u_w*vol_p)/(vol_p+vol_w)+(u_p*vol_s+u_s*vol_p)/(vol_s+vol_p)-u_p;
		u_wn_e=(u_p*vol_w+u_w*vol_p)/(vol_p+vol_w)+(u_p*vol_n+u_n*vol_p)/(vol_n+vol_p)-u_p;
		u_wb_e=(u_p*vol_w+u_w*vol_p)/(vol_p+vol_w)+(u_p*vol_b+u_b*vol_p)/(vol_b+vol_p)-u_p;
		u_wt_e=(u_p*vol_w+u_w*vol_p)/(vol_p+vol_w)+(u_p*vol_t+u_t*vol_p)/(vol_t+vol_p)-u_p;

		u_es_e=(u_p*vol_e+u_e*vol_p)/(vol_p+vol_e)+(u_p*vol_s+u_s*vol_p)/(vol_p+vol_s)-u_p;
		u_en_e=(u_p*vol_e+u_e*vol_p)/(vol_p+vol_e)+(u_p*vol_n+u_n*vol_p)/(vol_p+vol_n)-u_p;
		u_eb_e=(u_p*vol_e+u_e*vol_p)/(vol_p+vol_e)+(u_p*vol_b+u_b*vol_p)/(vol_p+vol_b)-u_p;
		u_et_e=(u_p*vol_e+u_e*vol_p)/(vol_p+vol_e)+(u_p*vol_t+u_t*vol_p)/(vol_p+vol_t)-u_p;

		u_bs_e=(u_p*vol_b+u_b*vol_p)/(vol_p+vol_b)+(u_p*vol_s+u_s*vol_p)/(vol_s+vol_p)-u_p;
		u_bn_e=(u_p*vol_b+u_b*vol_p)/(vol_p+vol_b)+(u_p*vol_n+u_n*vol_p)/(vol_n+vol_p)-u_p;

		u_ts_e=(u_p*vol_t+u_t*vol_p)/(vol_p+vol_t)+(u_p*vol_s+u_s*vol_p)/(vol_s+vol_p)-u_p;
		u_tn_e=(u_p*vol_t+u_t*vol_p)/(vol_p+vol_t)+(u_p*vol_n+u_n*vol_p)/(vol_n+vol_p)-u_p;

		
	}
	else
	{
		u_ws_e=(u_p*vol_ws+u_ws*vol_p+u_w*vol_s+u_s*vol_w)/(vol_ws+vol_p+vol_s+vol_w);
		u_wn_e=(u_p*vol_wn+u_wn*vol_p+u_w*vol_n+u_n*vol_w)/(vol_wn+vol_p+vol_n+vol_w);
		u_wb_e=(u_p*vol_wb+u_wb*vol_p+u_w*vol_b+u_b*vol_w)/(vol_wb+vol_p+vol_b+vol_w);
		u_wt_e=(u_p*vol_wt+u_wt*vol_p+u_w*vol_t+u_t*vol_w)/(vol_wt+vol_p+vol_t+vol_w);
			
		u_es_e=(u_p*vol_es+u_es*vol_p+u_e*vol_s+u_s*vol_e)/(vol_es+vol_p+vol_s+vol_e);
		u_en_e=(u_p*vol_en+u_en*vol_p+u_e*vol_n+u_n*vol_e)/(vol_en+vol_p+vol_n+vol_e);
		u_eb_e=(u_p*vol_eb+u_eb*vol_p+u_e*vol_b+u_b*vol_e)/(vol_eb+vol_p+vol_b+vol_e);
		u_et_e=(u_p*vol_et+u_et*vol_p+u_e*vol_t+u_t*vol_e)/(vol_et+vol_p+vol_t+vol_e);
			
		u_bs_e=(u_p*vol_bs+u_bs*vol_p+u_b*vol_s+u_s*vol_b)/(vol_bs+vol_p+vol_b+vol_s);
		u_bn_e=(u_p*vol_bn+u_bn*vol_p+u_b*vol_n+u_n*vol_b)/(vol_bn+vol_p+vol_b+vol_n);
			
		u_ts_e=(u_p*vol_ts+u_ts*vol_p+u_t*vol_s+u_s*vol_t)/(vol_ts+vol_p+vol_s+vol_t);
		u_tn_e=(u_p*vol_tn+u_tn*vol_p+u_t*vol_n+u_n*vol_t)/(vol_tn+vol_p+vol_n+vol_t);
	}
	
	re=1.0/(((r_e*vol_p+r_p*vol_e)/(vol_p+vol_e)));
	rw=1.0/(((r_w*vol_p+r_p*vol_w)/(vol_p+vol_w)));
	rs=1.0/(((r_s*vol_p+r_p*vol_s)/(vol_p+vol_s)));
	rn=1.0/(((r_n*vol_p+r_p*vol_n)/(vol_p+vol_n)));
	rt=1.0/(((r_t*vol_p+r_p*vol_t)/(vol_p+vol_t)));
	rb=1.0/(((r_b*vol_p+r_p*vol_b)/(vol_p+vol_b)));
	
	
	
		
	for(iFace=0;iFace<6;iFace++)
		for(comp=0;comp<3;comp++)
			betaLocal[iFace][comp]=*(beta+INDEX*18+3*iFace+comp);
	
	DiffFlux[0]/*West*/=	rw*((betaLocal[0][0]*(u_p-u_w)+betaLocal[0][1]*(u_wn_e-u_ws_e)+betaLocal[0][2]*(u_wt_e-u_wb_e)));
	
	DiffFlux[1]/*East*/=	re*((betaLocal[1][0]*(u_e-u_p)+betaLocal[1][1]*(u_en_e-u_es_e)+betaLocal[1][2]*(u_et_e-u_eb_e)));
	
	DiffFlux[2]/*South*/=	rs*((betaLocal[2][0]*(u_es_e-u_ws_e)+betaLocal[2][1]*(u_p-u_s)+betaLocal[2][2]*(u_ts_e-u_bs_e)));
	
	DiffFlux[3]/*North*/=	rn*((betaLocal[3][0]*(u_en_e-u_wn_e)+betaLocal[3][1]*(u_n-u_p)+betaLocal[3][2]*(u_tn_e-u_bn_e)));
	
	DiffFlux[4]/*Bottom*/=	rb*((betaLocal[4][0]*(u_eb_e-u_wb_e)+betaLocal[4][1]*(u_bn_e-u_bs_e)+betaLocal[4][2]*(u_p-u_b)));
	
	DiffFlux[5]/*Top*/=	rt*((betaLocal[5][0]*(u_et_e-u_wt_e)+betaLocal[5][1]*(u_tn_e-u_ts_e)+betaLocal[5][2]*(u_t-u_p)));
	
	
	/*printf("%12.8lf %12.8lf %12.8lf %12.8lf %12.8lf %12.8lf %12.8lf\n",DiffFlux[0],DiffFlux[1],DiffFlux[2],DiffFlux[3],DiffFlux[4],DiffFlux[5],diagonal);*/
	return (diagonal);
}

