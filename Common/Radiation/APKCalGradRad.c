#include "APKCalGradRad.h"

int CalGradRad(int index[3],double *domain,double *volumes,int *size,double *slabVariable[6],double *slabVolume[6],double *beta,double DiffFlux[6], double diagonal[6])
{
	
	int NK,NJ,NI,k,j,i,iFace,comp;
	double u_p,u_e,u_w,u_s,u_n,u_t,u_b,vol_p,vol_e,vol_w,vol_s,vol_n,vol_t,vol_b,u_ws,u_wn,u_wt,u_wb,u_es,
	u_en,u_et,u_eb,vol_ws,vol_wn,vol_wt,vol_wb,vol_es,vol_en,vol_et,vol_eb,u_ts,u_tn,u_bs,u_bn,vol_ts,
	vol_tn,vol_bs,vol_bn;
	double u_ws_e,u_wn_e,u_wb_e,u_wt_e,u_es_e,u_en_e,u_eb_e,u_et_e,u_bs_e,u_bn_e,u_ts_e,u_tn_e;
	double betaLocal[6][3],coeffWp[4],coeffEp[4],coeffNp[4],coeffSp[4],coeffTp[4],coeffBp[4];
	NI=*(size+0);
	NJ=*(size+1);
	NK=*(size+2);
	
	k=index[0],j=index[1],i=index[2];
	
	//printf("\n%d %d %d",NI,NJ,NK);
	u_p=*(domain+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
	vol_p=*(volumes+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
	
	/*printf("u_p= %12.8lf\n",u_p);*/
		
	if(i!=0)
	{
		u_w=*(domain+(NJ-1)*(NI-1)*k+(NI-1)*j+(i-1));/* cell center*/
		vol_w=*(volumes+(NJ-1)*(NI-1)*k+(NI-1)*j+(i-1));
	}
	else
	{
		u_w =*(slabVariable[1]+(NJ-1)*k+j);
		vol_w =*(slabVolume[1]+(NJ-1)*k+j);
	}
	if(i!=(NI-2))
	{
		u_e=*(domain+(NJ-1)*(NI-1)*k+(NI-1)*j+(i+1));
		vol_e=*(volumes+(NJ-1)*(NI-1)*k+(NI-1)*j+(i+1));
	}
	else
	{
		u_e= *(slabVariable[0]+(NJ-1)*k+j);
		vol_e=*(slabVolume[0]+(NJ-1)*k+j);
	}
	if(j!=0)
	{
		u_s=*(domain+(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i);
		vol_s=*(volumes+(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i);
	}
	else
	{
		u_s=*(slabVariable[3]+(NI-1)*k+i);
		vol_s=*(slabVolume[3]+(NI-1)*k+i);
	}
	if(j!=(NJ-2))
	{
		u_n=*(domain+(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+i);
		vol_n=*(volumes+(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+i);
	}
	else
	{
		u_n=*(slabVariable[2]+(NI-1)*k+i);
		vol_n=*(slabVolume[2]+(NI-1)*k+i);
	}
	if(k!=0)
	{
		u_b=*(domain+(NJ-1)*(NI-1)*(k-1)+(NI-1)*j+i);
		vol_b=*(volumes+(NJ-1)*(NI-1)*(k-1)+(NI-1)*j+i);
	}
	else
	{
		u_b=*(slabVariable[5]+(NI-1)*j+i);
		vol_b=*(slabVolume[5]+(NI-1)*j+i);
	}
	if(k!=(NK-2))
	{
		u_t=*(domain+(NJ-1)*(NI-1)*(k+1)+(NI-1)*j+i);
		vol_t=*(volumes+(NJ-1)*(NI-1)*(k+1)+(NI-1)*j+i);
	}
	else
	{
		u_t=*(slabVariable[4]+(NI-1)*j+i);
		vol_t=*(slabVolume[4]+(NI-1)*j+i);
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
				u_ws=*(domain+(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+(i-1));
				vol_ws=*(volumes+(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+(i-1));
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
				u_wn=*(domain+(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+(i-1));
				vol_wn=*(volumes+(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+(i-1));
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
				u_wb=*(domain+(NJ-1)*(NI-1)*(k-1)+(NI-1)*j+(i-1));
				vol_wb=*(volumes+(NJ-1)*(NI-1)*(k-1)+(NI-1)*j+(i-1));
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
				u_wt=*(domain+(NJ-1)*(NI-1)*(k+1)+(NI-1)*j+(i-1));
				vol_wt=*(volumes+(NJ-1)*(NI-1)*(k+1)+(NI-1)*j+(i-1));
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
				u_es=*(domain+(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+(i+1));
				vol_es=*(volumes+(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+(i+1));
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
				u_en=*(domain+(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+(i+1));
				vol_en=*(volumes+(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+(i+1));
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
				u_eb=*(domain+(NJ-1)*(NI-1)*(k-1)+(NI-1)*j+(i+1));
				vol_eb=*(volumes+(NJ-1)*(NI-1)*(k-1)+(NI-1)*j+(i+1));
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
				u_et=*(domain+(NJ-1)*(NI-1)*(k+1)+(NI-1)*j+(i+1));
				vol_et=*(volumes+(NJ-1)*(NI-1)*(k+1)+(NI-1)*j+(i+1));
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
				u_bs=*(domain+(NJ-1)*(NI-1)*(k-1)+(NI-1)*(j-1)+i);
				vol_bs=*(volumes+(NJ-1)*(NI-1)*(k-1)+(NI-1)*(j-1)+i);
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
				u_bn=*(domain+(NJ-1)*(NI-1)*(k-1)+(NI-1)*(j+1)+i);
				vol_bn=*(volumes+(NJ-1)*(NI-1)*(k-1)+(NI-1)*(j+1)+i);
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
				u_ts=*(domain+(NJ-1)*(NI-1)*(k+1)+(NI-1)*(j-1)+i);
				vol_ts=*(volumes+(NJ-1)*(NI-1)*(k+1)+(NI-1)*(j-1)+i);
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
				u_tn=*(domain+(NJ-1)*(NI-1)*(k+1)+(NI-1)*(j+1)+i);
				vol_tn=*(volumes+(NJ-1)*(NI-1)*(k+1)+(NI-1)*(j+1)+i);
			}
		}
	}
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
	
	coeffWp[0]=vol_n/(vol_wn+vol_p+vol_n+vol_w)
	,coeffWp[1]=vol_s/(vol_ws+vol_p+vol_s+vol_w)
	,coeffWp[2]=vol_t/(vol_wt+vol_p+vol_t+vol_w)
	,coeffWp[3]=vol_b/(vol_wb+vol_p+vol_b+vol_w);
	
	coeffEp[0]=vol_n/(vol_en+vol_p+vol_n+vol_e)
	,coeffEp[1]=vol_s/(vol_es+vol_p+vol_s+vol_e)
	,coeffEp[2]=vol_t/(vol_et+vol_p+vol_t+vol_e)
	,coeffEp[3]=vol_b/(vol_eb+vol_p+vol_b+vol_e); 
	
	coeffSp[0]=vol_e/(vol_es+vol_p+vol_s+vol_e)
	,coeffSp[1]=vol_w/(vol_ws+vol_p+vol_s+vol_w)
	,coeffSp[2]=vol_t/(vol_ts+vol_p+vol_s+vol_t)
	,coeffSp[3]=vol_b/(vol_bs+vol_p+vol_b+vol_s);
	
	coeffNp[0]=vol_e/(vol_en+vol_p+vol_n+vol_e)
	,coeffNp[1]=vol_w/(vol_wn+vol_p+vol_n+vol_w)
	,coeffNp[2]=vol_t/(vol_tn+vol_p+vol_n+vol_t)
	,coeffNp[3]=vol_b/(vol_bn+vol_p+vol_b+vol_n);
	
	coeffTp[0]=vol_e/(vol_et+vol_p+vol_t+vol_e)
	,coeffTp[1]=vol_w/(vol_wt+vol_p+vol_t+vol_w)
	,coeffTp[2]=vol_n/(vol_tn+vol_p+vol_n+vol_t)
	,coeffTp[3]=vol_s/(vol_ts+vol_p+vol_s+vol_t);
	
	coeffBp[0]=vol_e/(vol_eb+vol_p+vol_b+vol_e)
	,coeffBp[1]=vol_w/(vol_wb+vol_p+vol_b+vol_w)
	,coeffBp[2]=vol_n/(vol_bn+vol_p+vol_b+vol_n)
	,coeffBp[3]=vol_s/(vol_bs+vol_p+vol_b+vol_s);
	
	//printf("%12.8lf,%12.8lf,%12.8lf,%12.8lf,%12.8lf,%12.8lf,%12.8lf,%12.8lf,%12.8lf,%12.8lf,%12.8lf,%12.8lf,\n",vol_ws,vol_wn,vol_wb,vol_wt,vol_es,vol_en,vol_eb,vol_et,vol_bs,vol_bn,vol_ts,vol_tn);
	//printf("%12.8lf,%12.8lf,%12.8lf,%12.8lf,%12.8lf,%12.8lf,%12.8lf,%12.8lf,%12.8lf,%12.8lf,%12.8lf,%12.8lf,\n",u_ws_e,u_wn_e,u_wb_e,u_wt_e,u_es_e,u_en_e,u_eb_e,u_et_e,u_bs_e,u_bn_e,u_ts_e,u_tn_e);
	
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

		coeffWp[0]=vol_p/(vol_p+vol_w)
		,coeffWp[1]=vol_p/(vol_p+vol_w)
		,coeffWp[2]=vol_p/(vol_p+vol_w)
		,coeffWp[3]=vol_p/(vol_p+vol_w);
	
		coeffEp[0]=vol_p/(vol_p+vol_e)
		,coeffEp[1]=vol_p/(vol_p+vol_e)
		,coeffEp[2]=vol_p/(vol_p+vol_e)
		,coeffEp[3]=vol_p/(vol_p+vol_e); 
		
		coeffSp[0]=vol_p/(vol_p+vol_s)
		,coeffSp[1]=vol_p/(vol_p+vol_s)
		,coeffSp[2]=vol_p/(vol_p+vol_s)
		,coeffSp[3]=vol_p/(vol_p+vol_s);
	
		coeffNp[0]=vol_p/(vol_p+vol_n)
		,coeffNp[1]=vol_p/(vol_p+vol_n)
		,coeffNp[2]=vol_p/(vol_p+vol_n)
		,coeffNp[3]=vol_p/(vol_p+vol_n);
		
		coeffTp[0]=vol_p/(vol_p+vol_t)
		,coeffTp[1]=vol_p/(vol_p+vol_t)
		,coeffTp[2]=vol_p/(vol_p+vol_t)
		,coeffTp[3]=vol_p/(vol_p+vol_t);
	
		coeffBp[0]=vol_p/(vol_p+vol_b)
		,coeffBp[1]=vol_p/(vol_p+vol_b)
		,coeffBp[2]=vol_p/(vol_p+vol_b)
		,coeffBp[3]=vol_p/(vol_p+vol_b);
	}                                                          
	
	for(iFace=0;iFace<6;iFace++)
		for(comp=0;comp<3;comp++)
			betaLocal[iFace][comp]=*(beta+((NI-1)*(NJ-1)*k+(NI-1)*j+i)*6*3+3*iFace+comp);
	
	DiffFlux[0]/*West*/=	(betaLocal[0][0]*(u_p-u_w)+betaLocal[0][1]*(u_wn_e-u_ws_e)+betaLocal[0][2]*(u_wt_e-u_wb_e));
	
	DiffFlux[1]/*East*/=	(betaLocal[1][0]*(u_e-u_p)+betaLocal[1][1]*(u_en_e-u_es_e)+betaLocal[1][2]*(u_et_e-u_eb_e));
	
	DiffFlux[2]/*South*/=	(betaLocal[2][0]*(u_es_e-u_ws_e)+betaLocal[2][1]*(u_p-u_s)+betaLocal[2][2]*(u_ts_e-u_bs_e));
	
	DiffFlux[3]/*North*/=	(betaLocal[3][0]*(u_en_e-u_wn_e)+betaLocal[3][1]*(u_n-u_p)+betaLocal[3][2]*(u_tn_e-u_bn_e));
	
	DiffFlux[4]/*Bottom*/=	(betaLocal[4][0]*(u_eb_e-u_wb_e)+betaLocal[4][1]*(u_bn_e-u_bs_e)+betaLocal[4][2]*(u_p-u_b));
	
	DiffFlux[5]/*Top*/=	(betaLocal[5][0]*(u_et_e-u_wt_e)+betaLocal[5][1]*(u_tn_e-u_ts_e)+betaLocal[5][2]*(u_t-u_p));
	
	diagonal[0]/*West*/=betaLocal[0][0]+betaLocal[0][1]*(coeffWp[1]-coeffWp[0])+betaLocal[0][2]*(coeffWp[3]-coeffWp[2]);
		    
	diagonal[1]/*East*/=betaLocal[1][0]+betaLocal[1][1]*(coeffEp[0]-coeffEp[1])+betaLocal[1][2]*(coeffEp[2]-coeffEp[3]);
	            
	diagonal[2]/*South*/=betaLocal[2][0]*(coeffSp[1]-coeffSp[0])+betaLocal[2][1]+betaLocal[2][2]*(coeffSp[3]-coeffSp[2]);
		   
	diagonal[3]/*North*/=betaLocal[3][0]*(coeffNp[0]-coeffNp[1])+betaLocal[3][1]+betaLocal[3][2]*(coeffNp[2]-coeffNp[3]);
		      
	diagonal[4]/*Bottom*/=betaLocal[4][0]*(coeffBp[1]-coeffBp[0])+betaLocal[4][1]*(coeffBp[3]-coeffBp[2])+betaLocal[4][2];
	
	diagonal[5]/*Top*/=betaLocal[5][0]*(coeffTp[0]-coeffTp[1])+betaLocal[5][1]*(coeffTp[2]-coeffTp[3])+betaLocal[5][2];	 
	/*printf("%12.8lf %12.8lf %12.8lf %12.8lf %12.8lf %12.8lf %12.8lf\n",DiffFlux[0],DiffFlux[1],DiffFlux[2],DiffFlux[3],DiffFlux[4],DiffFlux[5],diagonal);*/
	return 0;
}
