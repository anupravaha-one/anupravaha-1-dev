#include "APKDiffusion2dWithConstDiffCoeff_opti.h"
/**
The function Diffusion2dWithConstDiffCoeff_opti() is used to compute the diffusion flux of the form (Del.(Gamma Del(phi))) for 2 dimensional domain, where the diffusion coefficient Gamma is a constant. This function gets all the neighbouring cell pointers either from the domain or from the zone slabs, depending upon the cell index, and computes the diffusion flux on all the four faces of the cell. Volumetric interpolation of the diffusion coefficient is not done for computing the diffusion flux through a face. See Eqn(2.22) to Eqn(2.31) in Project Report No 1, Dec 2005, for further details.

The above statements are copied from APKDiffusion2dWithConstDiffCoeff function.  This function in written in an efficient way by avoiding calculation of indices and etc.
*/
double Diffusion2dWithConstDiffCoeff_opti(int index[3],double *domain,double *volumes,int *size,double *slabVariable[6],double *slabVolume[6],double *diffusionCoeffdomain,double *beta,double DiffFlux[6],double DiagonalN[6])
{
	int NK,NJ,NI,k,j,i,iFace,comp,INDEX;
	double u_p,u_e,u_w,u_s,u_n,vol_p,vol_e,vol_w,vol_s,vol_n,u_ws,u_wn,u_es,u_en;
	double vol_ws,vol_wn,vol_es,vol_en;
	double u_ws_e,u_wn_e,u_es_e,u_en_e;
	double betaLocal[6][3],diffCoeff,diagonal;
	
	
	
	NI=*(size+0);
	NJ=*(size+1);
	NK=*(size+2);
	
	diagonal = 0.0;
	
	k=index[0],j=index[1],i=index[2];
	INDEX = (NJ-1)*(NI-1)*k+(NI-1)*j+i;
	diffCoeff = *diffusionCoeffdomain;
	
	//printf("\nDiffusionCoeff=%lf\n",diffCoeff);
	u_p=*(domain+INDEX);
	vol_p=*(volumes+INDEX);
		
	if(i!=0)
	{
		u_w=*(domain+INDEX-1);
		vol_w=*(volumes+INDEX-1);
	}
	else
	{
		u_w =*(slabVariable[1]+j);
		vol_w =*(slabVolume[1]+j);
	}
	if(i!=(NI-2))
	{
		u_e=*(domain+INDEX+1);
		vol_e=*(volumes+INDEX+1);
	}
	else
	{
		u_e= *(slabVariable[0]+j);
		vol_e=*(slabVolume[0]+j);
	}
	if(j!=0)
	{
		u_s=*(domain+INDEX-(NI-1));
		vol_s=*(volumes+INDEX-(NI-1));
	}
	else
	{
		u_s=*(slabVariable[3]+i);
		vol_s=*(slabVolume[3]+i);
	}
	if(j!=(NJ-2))
	{
		u_n=*(domain+INDEX+(NI-1));
		vol_n=*(volumes+INDEX+(NI-1));
	}
	else
	{
		u_n=*(slabVariable[2]+i);
		vol_n=*(slabVolume[2]+i);
	}
	
	if(!(i==0 && j==0))
	{
		if(i==0)
		{
			u_ws=*(slabVariable[1]+(j-1));/*from west slab*/
			vol_ws=*(slabVolume[1]+(j-1));
		}
		else
		{
			if(j==0)
			{
				u_ws=*(slabVariable[3]+(i-1));/*from southern slab*/
				vol_ws=*(slabVolume[3]+(i-1));
			}
			else
			{
				u_ws=*(domain+INDEX-(NI-1)-1);
				vol_ws=*(volumes+INDEX-(NI-1)-1);
			}
		}
	}
		
	if(!(i==0 && j==(NJ-2)))
	{
		if(i==0)
		{
			u_wn=*(slabVariable[1]+(j+1));/*from western slab*/
			vol_wn=*(slabVolume[1]+(j+1));
		}
		else
		{
			if(j==(NJ-2))
			{
				u_wn=*(slabVariable[2]+(i-1));/*from northern slab*/
				vol_wn=*(slabVolume[2]+(i-1));
			}
			else
			{
				u_wn=*(domain+INDEX+(NI-1)-1);
				vol_wn=*(volumes+INDEX+(NI-1)-1);
			}
		}
	}
	
	if(!(i==(NI-2) && j==0))
	{
		if(i==(NI-2))
		{
			u_es=*(slabVariable[0]+(j-1));
			vol_es=*(slabVolume[0]+(j-1));
		}
		else
		{
			if(j==0)
			{
				u_es=*(slabVariable[3]+(i+1));
				vol_es=*(slabVolume[3]+(i+1));
			}
			else
			{
				u_es=*(domain+INDEX-(NI-1)+1);
				vol_es=*(volumes+INDEX-(NI-1)+1);
			}
		}
	}
		
	if(!(i==(NI-2) && j==(NJ-2)))
	{
		if(i==(NI-2))
		{
			u_en=*(slabVariable[0]+(j+1));
			vol_en=*(slabVolume[0]+(j+1));
		}
		else
		{
			if(j==(NJ-2))
			{
				u_en=*(slabVariable[2]+(i+1));
				vol_en=*(slabVolume[2]+(i+1));
			}
			else
			{
				u_en=*(domain+INDEX+(NI-1)+1);
				vol_en=*(volumes+INDEX+(NI-1)+1);
			}
		}
	}
	
	
	u_ws_e=(u_p*vol_ws+u_ws*vol_p+u_w*vol_s+u_s*vol_w)/(vol_ws+vol_p+vol_s+vol_w);
	u_wn_e=(u_p*vol_wn+u_wn*vol_p+u_w*vol_n+u_n*vol_w)/(vol_wn+vol_p+vol_n+vol_w);
		
	u_es_e=(u_p*vol_es+u_es*vol_p+u_e*vol_s+u_s*vol_e)/(vol_es+vol_p+vol_s+vol_e);
	u_en_e=(u_p*vol_en+u_en*vol_p+u_e*vol_n+u_n*vol_e)/(vol_en+vol_p+vol_n+vol_e);
		
	//printf("%12.8lf,%12.8lf,%12.8lf,%12.8lf,%12.8lf,%12.8lf,%12.8lf,%12.8lf,%12.8lf,%12.8lf,%12.8lf,%12.8lf,\n",vol_ws,vol_wn,vol_wb,vol_wt,vol_es,vol_en,vol_eb,vol_et,vol_bs,vol_bn,vol_ts,vol_tn);
	//printf("%12.8lf,%12.8lf,%12.8lf,%12.8lf,%12.8lf,%12.8lf,%12.8lf,%12.8lf,%12.8lf,%12.8lf,%12.8lf,%12.8lf,\n",u_ws_e,u_wn_e,u_wb_e,u_wt_e,u_es_e,u_en_e,u_eb_e,u_et_e,u_bs_e,u_bn_e,u_ts_e,u_tn_e);
	
	if(((i==0)&&(j==0))||((i==0)&&(k==0))||((i==0)&&(j==(NJ-2)))||((i==0)&&(k==(NK-2)))||((i==(NI-2))&&(j==0))||((i==(NI-2))&&(k==0))||((i==(NI-2))&&(j==(NJ-2)))||((i==(NI-2))&&(k==(NK-2)))||((j==0)&&(k==0))||((j==0)&&(k==(NK-2)))||((j==(NJ-2))&&(k==0))||((j==(NJ-2))&&(k==(NK-2))))
	///This condition needs to be revisited as three point interpolation technique is not necessary for slabcells for 2-d domain.  Applying this technique to corner points will suffice our requirement
	//if(((i==0)&&(j==0))||((i==0)&&(j==(NJ-2)))||((i==(NI-2))&&(j==0))||((i==(NI-2))&&(j==(NJ-2))))
	{
		u_ws_e=(u_p*vol_w+u_w*vol_p)/(vol_p+vol_w)+(u_p*vol_s+u_s*vol_p)/(vol_s+vol_p)-u_p;
		u_wn_e=(u_p*vol_w+u_w*vol_p)/(vol_p+vol_w)+(u_p*vol_n+u_n*vol_p)/(vol_n+vol_p)-u_p;
		u_es_e=(u_p*vol_e+u_e*vol_p)/(vol_p+vol_e)+(u_p*vol_s+u_s*vol_p)/(vol_p+vol_s)-u_p;
		u_en_e=(u_p*vol_e+u_e*vol_p)/(vol_p+vol_e)+(u_p*vol_n+u_n*vol_p)/(vol_p+vol_n)-u_p;
		//printf("Entered \tINDEX=%d\n",INDEX);	
	}
	
		
	for(iFace=0;iFace<4;iFace++)
	for(comp=0;comp<2;comp++)
		betaLocal[iFace][comp]=*(beta+INDEX*18+3*iFace+comp);
	
	DiffFlux[0]/*West*/=	diffCoeff*(betaLocal[0][0]*(u_p-u_w)+betaLocal[0][1]*(u_wn_e-u_ws_e));//+betaLocal[0][2]*(u_wt_e-u_wb_e));
	
	DiffFlux[1]/*East*/=	diffCoeff*(betaLocal[1][0]*(u_e-u_p)+betaLocal[1][1]*(u_en_e-u_es_e));//+betaLocal[1][2]*(u_et_e-u_eb_e));
	
	DiffFlux[2]/*South*/=	diffCoeff*(betaLocal[2][0]*(u_es_e-u_ws_e)+betaLocal[2][1]*(u_p-u_s));//+betaLocal[2][2]*(u_ts_e-u_bs_e));
	
	DiffFlux[3]/*North*/=	diffCoeff*(betaLocal[3][0]*(u_en_e-u_wn_e)+betaLocal[3][1]*(u_n-u_p));//+betaLocal[3][2]*(u_tn_e-u_bn_e));
	
	DiffFlux[4]/*Bottom*/=	diffCoeff*0.0;//(betaLocal[4][0]*(u_eb_e-u_wb_e)+betaLocal[4][1]*(u_bn_e-u_bs_e)+betaLocal[4][2]*(u_p-u_b));
	
	DiffFlux[5]/*Top*/=	diffCoeff*0.0;//(betaLocal[5][0]*(u_et_e-u_wt_e)+betaLocal[5][1]*(u_tn_e-u_ts_e)+betaLocal[5][2]*(u_t-u_p));
	
	//printf("Diffflux[0]=%e\tDiffflux[1]=%e\tDiffflux[2]=%e\tDiffflux[3]=%e\tDiffflux[4]=%e\tDiffflux[5]=%e\n",DiffFlux[0],DiffFlux[1],DiffFlux[2],DiffFlux[3],DiffFlux[4],DiffFlux[5]);
		 
	//printf("%12.8lf %12.8lf %12.8lf %12.8lf %12.8lf %12.8lf %12.8lf\n",DiffFlux[0],DiffFlux[1],DiffFlux[2],DiffFlux[3],DiffFlux[4],DiffFlux[5],diagonal);
	return (diagonal);
}

