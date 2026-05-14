#include "APKUpdateBCslabsKunsat.h"	 

/**
This function is used to update the unsaturated relative permeabilty for unsaturated and varablly satuarted region in case of porous media(soil) in boundary face. It update unsatuarated permeabilty based on pressure head and selected model for soils
*/

int UpdateBCslabsKunsat(Base *base,Zone *z,ZoneBoundary *ZoneBc,int nv,int vartag,int direct_index) 
        {       
		
		int  maxn,nki,njk,nij,ni,nj,nk,size,k,j,i,index,*ptrRange;
		double *variable;
		
		variable = FindSlabPtrBC(ZoneBc,vartag);
		 nv = find_nv(base,nv);
		ptrRange = ZoneBc->Point_Range;
				
		ni=*(z->Size+0);
		nj=*(z->Size+1);
		nk=*(z->Size+2); 
		       
		
                nij=ni*nj;
		nki=nk*ni;
		njk=nj*nk;
		
		
		if(nij>njk)
		{
		if(nij>nki)
		maxn=nij;
		else
		maxn=nki;
		}
		else
		{
			if(njk>nki)
			maxn=njk;

			else
			maxn=nki;
		}

		
				
		if (direct_index==1||direct_index==2)
		{
		size=(nk-1)*(nj-1);
		//NueRms=boundary_d_t_s(base,direct_index,z->Flowsol->H,19,z,size,(ZoneBc->Data_Set+nv),((ZoneBc->ZSlab)->H)); 
		
		 switch(direct_index)
		 {
		 case 1: 
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(nj-1)+j;
		                        *(variable+index)=FindKunsat(index,((ZoneBc->ZSlab)->H),z,ZoneBc);
		
					}
					break;
		 case 2: 
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(nj-1)+j;
				         	*(variable+index) =FindKunsat(index,((ZoneBc->ZSlab)->H),z,ZoneBc);
			//printf("\n var=%lf   \tdomain=%lf",*(variable+index),*(domain+(nj-1)*(ni-1)*k+(ni-1)*j+i));
					}
				  break;                 
		}	      
		
		}
					
			
			
		else if (direct_index==3||direct_index==4)
		{
			
			  size=(nk-1)*(ni-1);
				      
			 switch(direct_index) 
			 {
			case 3:
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = k*(ni-1)+i;
						
					   *(variable+index)=FindKunsat(index,((ZoneBc->ZSlab)->H),z,ZoneBc);
						
					}
				break;

			case 4:
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = k*(ni-1)+i;
						
					  *(variable+index)=FindKunsat(index,((ZoneBc->ZSlab)->H),z,ZoneBc);
				
					}
				break;
                          }
		}
		
		else if (direct_index==5||direct_index==6)	
		{
			size=(nj-1)*(ni-1);
			
			
			 switch(direct_index)
			 {
			
			case 5:/*Dir 5 ->Kmax,Top*/
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = j*(ni-1)+i;
						
				    *(variable+index)=FindKunsat(index,((ZoneBc->ZSlab)->H),z,ZoneBc);
					}
				  break;

			case 6:/*Dir 6 ->Kmin,Bottom*/
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = j*(ni-1)+i;
						
				   *(variable+index)=FindKunsat(index,((ZoneBc->ZSlab)->H),z,ZoneBc);
					}
			             break;
			}	      
		}		      
				      
			
			
			return(0);	
		}	
		  
		  
		  
		  		  
		  
		  
	double	FindKunsat(int indexH,double *Head,Zone *z,ZoneBoundary *ZoneBc)
		{
		   
	        double tert,alpha,tert1,tert2,theta_residual,theta_saturation,n,m,Diff,Saturation,arbi;
		double KunSat,theta,A;
		int flagtype;	      
		
		flagtype=z->flagtype_model;
		switch(flagtype)	
		{ 
		       
		
		case 1:
		
		           	    	 alpha=z->Alpha;
					theta_residual=z->Theta_residual;
					theta_saturation=z->Theta_saturation;
					n=z->n;
					m=1.0-1.0/n;
					
					Diff=theta_saturation-theta_residual;	
		     
				        if(*(Head+indexH)>=0)
					KunSat=1.0;
					
					else
					{
					tert=*(Head+indexH)*alpha;
					tert=fabs(tert);
					tert1=pow(tert,n);
					tert2=pow(1+tert1,-m-1);
					Saturation=pow(1+tert1,-m);
					theta=Saturation*(Diff)+theta_residual;
					arbi=1-pow(Saturation,1/m);
					arbi=1-pow(arbi,m);
					arbi=pow(arbi,2.0);
					KunSat=arbi*pow(Saturation,.5);
					}
		   break;
		   
		case 2:
		
		        alpha=z->A;
			theta_residual=z->Theta_residual;
			theta_saturation=z->Theta_saturation;
			m=z->B;
			n=z->D;
			A=z->C;
			
			if(*(Head+indexH)>=0)
			KunSat=1;
			else
			{

			tert=*(Head+indexH);
			tert=fabs(tert)*100;
			tert1=alpha+pow(tert,m);
			arbi=A+pow(tert,n);
			KunSat=A/arbi;


			}   
		    break;
		   
		 }
		 
		return(KunSat);
	}


