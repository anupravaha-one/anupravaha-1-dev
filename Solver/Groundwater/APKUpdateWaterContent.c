#include "APKUpdateWaterContent.h"
/**

This function update the theta(water content) value after every one iteration based on Model selcted i.e VGM or haverlamp. It also update the the unsaturated permeabilty and gradient of theta(watercontent)w.r.t pressure-head 
*/
int UpdateWaterContent(Base *base)
{
       int NI,NJ,NK,i,j,k,*ptrSize,indexH,numZones,iZone,flagtype;
       double *Head,*theta,*KunSat,*Gradient,*thetaPr,*HeadPr;
       double alpha,Saturation,theta_residual,theta_saturation,m,n,t1,t2,h1,h2;// these are all parameter that has to defined there typical value differ from soil to soil 
       Zone *zone;
    
       double arbi,tert,tert1,tert2,Diff,Equi,MN,A;
      
       numZones=base->Number_of_zones;
             
	   
		 
		for(iZone=0;iZone<numZones;iZone++)
		{
		zone=base->CgnsZone+iZone;		
		
		ptrSize=zone->Size;
		
		theta=zone->Flowsol->WaterContent;
		
		thetaPr=zone->PrevFlowsol->WaterContent;
		
		
		Head=zone->Flowsol->H;
		
		HeadPr=zone->PrevFlowsol->H;
		
		ptrSize=zone->Size;
		
		KunSat=zone->Flowsol->UnSatPermeability;
		
		Gradient=zone->Flowsol->GradTheta;
		
		flagtype=zone->flagtype_model;
		
		NI= *(ptrSize+0);
		NJ= *(ptrSize+1);
		NK= *(ptrSize+2);
		
				
		switch(flagtype)	
		{ 
		
			case 1://clay  form Kirkland et al 1992 K=(satuarted)=.131m/day 
			
			alpha=zone->Alpha;
			theta_residual=zone->Theta_residual;
			theta_saturation=zone->Theta_saturation;
			n=zone->n;
			m=1.0-1.0/n;
			
			
		MN=m*n*pow(alpha,n);
		
		Diff=theta_saturation-theta_residual;	
		
		Equi=Diff*MN;
		
				
		for(k=0;k<NK-1;k++)
		for(j=0;j<NJ-1;j++)
		for(i=0;i<NI-1;i++)
		{    
					indexH=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
					if(*(Head+indexH)>=0)
					{ 
					*(KunSat+indexH)=1.00;
					*(theta+indexH)=theta_saturation;
					*(Gradient+indexH)=10e-07;
					}
					else
					{
					
					tert= *(Head+indexH)*alpha; 
					tert=fabs(tert);
					tert1=pow(tert,n);
					tert2=pow(1+tert1,-m-1);
					Saturation=pow(1+tert1,-m);
					*(theta+indexH)=Saturation*(Diff)+theta_residual;
					arbi=1-pow(Saturation,1/m);
					arbi=1-pow(arbi,m);
					arbi=pow(arbi,2.0);
				      *(KunSat+indexH)=arbi*pow(Saturation,0.5);
					t1=*(theta+indexH);
					t2=*(thetaPr+indexH);
					h1=*(Head+indexH);
					h2=*(HeadPr+indexH);
					
					if((fabs(t1-t2)<=10e-8)||(fabs(h1-h2)<=10e-18))
					*(Gradient+indexH)=Equi*pow(fabs(*(Head+indexH)),n-1)*(tert2) ;
					else
					*(Gradient+indexH)=((t1-t2)/(h1-h2));
				
				
				        }
					
					
		}
		
		break;
		
		case 2:
		        
			
			
			alpha=zone->A;
			theta_residual=zone->Theta_residual;
			theta_saturation=zone->Theta_saturation;
			m=zone->B;
			n=zone->D;
			A=zone->C;
				
			MN=m*n*pow(alpha,n);
			
			Diff=theta_saturation-theta_residual;	
			
			Equi=Diff*MN;	
		                	
		
				for(k=0;k<NK-1;k++)
				for(j=0;j<NJ-1;j++)
				for(i=0;i<NI-1;i++)
				{    
					indexH=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
					if(*(Head+indexH)>=0)
					{ 
					*(KunSat+indexH)=1;
					*(theta+indexH)=theta_saturation;//1
					*(Gradient+indexH)=10e-07;
					}
					else
					{
					
					tert=*(Head+indexH);
					tert=fabs(tert)*100;
					tert1=alpha+pow(tert,m);
					*(theta+indexH)=alpha*(Diff)/tert1+theta_residual;
					arbi=A+pow(tert,n);
					*(KunSat+indexH)=A/arbi;
					 
					  t1=*(theta+indexH);
				           t2=*(thetaPr+indexH);
				         h1=*(Head+indexH);
				          h2=*(HeadPr+indexH);
		//*(Gradient+indexH)=m*Diff*alpha*pow(tert1,-2)*pow(tert,m-1);
				         
					 if((fabs(t1-t2)<=10e-8)||(fabs(h1-h2)<=10e-18))
					*(Gradient+indexH)=*(Gradient+indexH);
				         else
					*(Gradient+indexH)=((t1-t2)/(h1-h2));
					}
					
		                  }
		break;
				
		}
	      }
		return 0;
}
