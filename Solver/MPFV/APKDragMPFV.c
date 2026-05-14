#include "APKDragMPFV.h"
#define EPSILON 1e-9
/**
This function calculates the drag force for the case of disperse phase flow. The drag is calculated and stored in the data structure under the source term. This function is called once 
per time step
 */
double DragMPFV(Base *base,int CounterTstep)
{
	double Re,Cd,mu,rho,rho1,v1,v2,dia;
	double modRelV,temp,beta,B,*ptrVolume;
	int cntZ,NI,NJ,NK,*ptrSize,numZones,i,j,k;
	Zone *zone;
	int flagDen,flagMu,flagDia;//Flags
	int totalcells,index,nophase,nop;/*,flagrhoDisperse=0;*/
	double *ptrPoint1X,*ptrPoint1Y,*ptrPoint1Z,sum1=0.0,sum2=0.0,sum3=0.0,C_g=0.0,sumzone=0.0;
	double *coordX,*coordY,*coordZ,centroid[3];
	double V_tot_c_X=0.0,V_tot_c_Y=0.0,V_tot_c_Z=0.0,V_tot_d_X=0.0,V_tot_d_Y=0.0,V_tot_d_Z=0.0,Uavg_c_X=0.0,Uavg_c_Y=0.0,Uavg_c_Z=0.0,Uavg_d_X=0.0,Uavg_d_Y=0.0,Uavg_d_Z=0.0,M_tot_c=0.0,M_tot_d=0.0;
	numZones= base->Number_of_zones;
	FILE *CG;
	CG=fopen("CG.dat","a");
	nop = base->nphase;
	for(nophase=0;nophase<nop;nophase++)
	{	
		flagMu = base->BaseParameters->flagViscoeffPhase[nophase];
		flagDen = base->BaseParameters->flagRhoPhase[nophase];
	}
	for(cntZ=0;cntZ<numZones;cntZ++)
	{
		zone = (base->CgnsZone+cntZ);
		ptrSize = zone->Size;
		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);
		totalcells = (NI-1)*(NJ-1)*(NK-1);
		
		flagDia=base->BaseParameters->flagVariableDia;
		ptrVolume= zone->ZGeom->Volume;
		index=0;
		for(k=0;k<NK-1;k++)
		for(j=0;j<NJ-1;j++)
		for(i=0;i<NI-1;i++)
		{	
			
		
			v1 = *(zone->PrevFlowsol->VarthetaPhase[0]+index);
			v2 = *(zone->PrevFlowsol->VarthetaPhase[1]+index);
			
			mu= *(zone->Flowsol->ViscosityofPhase[0]+index*flagMu);
			rho= *(zone->Flowsol->RhoofPhase[0]+index*flagDen);
			//dia= .002;//*(zone->PrevFlowsol->ParticleDia+index*flagDia);
			dia = *(zone->Flowsol->ParticleDia+index*flagDia);
			
			rho1= *(zone->Flowsol->RhoofPhase[1]+index*flagDen);

			//new part
			V_tot_c_X = V_tot_c_X+(*((zone->PrevFlowsol)->Momentumofphasex[0]+index)*(*(ptrVolume+index)));
			V_tot_c_Y = V_tot_c_Y +(*((zone->PrevFlowsol)->Momentumofphasey[0]+index)*(*(ptrVolume+index)));
			
			if((zone->PrevFlowsol->Momentumofphasez[0])!=NULL)
				V_tot_c_Z = V_tot_c_Z+(*((zone->PrevFlowsol)->Momentumofphasez[0]+index)*(*(ptrVolume+index)));
			
			
			V_tot_d_X = V_tot_d_X +(*((zone->PrevFlowsol)->Momentumofphasex[1]+index)*(*(ptrVolume+index)));
			V_tot_d_Y = V_tot_d_Y+(*((zone->PrevFlowsol)->Momentumofphasey[1]+index)*(*(ptrVolume+index)));
			
			if((zone->PrevFlowsol->Momentumofphasez[1])!=NULL)
				V_tot_d_Z = V_tot_d_Z+(*((zone->PrevFlowsol)->Momentumofphasez[1]+index)*(*(ptrVolume+index)));
			
			
			M_tot_c = M_tot_c+rho*v1*(*(ptrVolume+index));
			M_tot_d = M_tot_d+rho1*v2*(*(ptrVolume+index));
			index++;
		}
			
	}
	
	Uavg_c_X = V_tot_c_X/M_tot_c;//printf("\nUavg_c_X=%lf",Uavg_c_X);
	Uavg_c_Y = V_tot_c_Y/M_tot_c;//printf("\nUavg_c_Y=%lf",Uavg_c_Y);
	if((zone->PrevFlowsol->Momentumofphasez[0])!=NULL)
		Uavg_c_Z = V_tot_c_Z/M_tot_c;//printf("\nUavg_c_Z=%lf",Uavg_c_Z);
		
	Uavg_d_X = V_tot_d_X/M_tot_d;//printf("\nM_tot_d=%lf",M_tot_d);
	Uavg_d_Y = V_tot_d_Y/M_tot_d;//printf("\nM_tot_d=%lf",M_tot_d);
	if((zone->PrevFlowsol->Momentumofphasez[1])!=NULL)
		Uavg_d_Z = V_tot_d_Z/M_tot_d;//printf("\nUavg_d_Z=%lf",Uavg_d_Z);
		
		
	modRelV = (Uavg_c_X-Uavg_d_X)*(Uavg_c_X-Uavg_d_X)+(Uavg_c_Y-Uavg_d_Y)*(Uavg_c_Y-Uavg_d_Y)+(Uavg_c_Z-Uavg_d_Z)*(Uavg_c_Z-Uavg_d_Z);
			
	temp=pow(modRelV,0.5);
	modRelV=temp;
	//printf("\nmodRelV=%e",modRelV);
	Re = (rho*dia*modRelV)/mu;printf("\nRe=%lf\n",Re);
	
	if(base->constant==1)
		Cd=userDefinedCD( Re);
	else
		Cd = base->BaseParameters->ConstantCD;
// 		if(Re>1e-9&&Re<0.49)
// 		{
// 			Cd=16.0/Re;
// 				//Cd=(24.0*(1+0.15*pow(Re,0.687)))/Re;
// 		}
// 		else if(Re>0.49&&Re<100.0)
// 			Cd=20.0/(pow(Re,0.643));
// 		else if(Re>100.0)
// 			Cd = 6.3/(pow(Re,0.385));
// 		else
// 			Cd=0.0;
			
// 		Cd=1.0;
			
	//printf("\ncd=%e,rho=%e,modRelV=%e,dia=%e",Cd,rho,modRelV,dia);
	beta = (0.75*Cd*rho*modRelV)/dia;//printf("\nbeta=%lf",beta);
	printf("\nCd used=%lf\n",Cd);
	for(cntZ=0;cntZ<numZones;cntZ++)
	{
		zone = (base->CgnsZone+cntZ);
		ptrSize = zone->Size;
		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);
		totalcells = (NI-1)*(NJ-1)*(NK-1);
		
		flagDia=base->BaseParameters->flagVariableDia;
		ptrVolume = zone->ZGeom->Volume;	
		index=0;
		for(k=0;k<NK-1;k++)
		for(j=0;j<NJ-1;j++)
		for(i=0;i<NI-1;i++)
		{	
			rho1= *(zone->Flowsol->RhoofPhase[1]+index*flagDen);
			v2 = *(zone->PrevFlowsol->VarthetaPhase[1]+index);
			rho= *(zone->Flowsol->RhoofPhase[0]+index*flagDen);
			if(base->BaseParameters->Bubbleapprox==0)
			{	
				
				*(zone->DragX+index) = beta*v2*Uavg_c_X;
				*(zone->DragY+index) = beta*v2*Uavg_c_Y;
				*(zone->DragZ+index) = beta*v2*Uavg_c_Z;
			}		
			//For small bubble approximation
			else
			{	
				v1 = *(zone->PrevFlowsol->VarthetaPhase[0]+index);
				
				*(zone->DragX+index) = beta*v2*(*((zone->PrevFlowsol)->Momentumofphasex[0]+index)/(rho*v1));
				
				*(zone->DragY+index) = beta*v2*(*((zone->PrevFlowsol)->Momentumofphasey[0]+index)/(rho*v1));
				
				if((zone->PrevFlowsol->Momentumofphasez[0])!=NULL)
					*(zone->DragZ+index) = beta*v2*(*((zone->PrevFlowsol)->Momentumofphasez[0]+index)/(rho*v1));
			}
			//End of small bubble part
			
			*(zone->Dragcoeff+index) = beta/rho1;
			
			B = v2*(rho-rho1);
			*(zone->Bx+index)=B;
				
	
			index++;
		}
		
	}

	//to find Centre of  gravity.//
	for(cntZ=0;cntZ<numZones;cntZ++)
	{
		zone = (base->CgnsZone+cntZ);
		ptrSize = zone->Size;
		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);
		ptrVolume= zone->ZGeom->Volume;
		coordX=(zone->CgnsGrid)->X;
		coordY=(zone->CgnsGrid)->Y; 
		coordZ=(zone->CgnsGrid)->Z;
		index=0;
		sum1=0.0;
		sum2=0.0;
		for(k=0;k<(NK-1);k++)
			for(j=0;j<(NJ-1);j++)
				for(i=0;i<(NI-1);i++)
		{	
			
			v2 = *(zone->PrevFlowsol->VarthetaPhase[1]+index);
			if((j==0)&&(i>=0))
			{	
				ptrPoint1X=coordX+NI*NJ*k+NI*j+i;
				ptrPoint1Y=coordY+NI*NJ*k+NI*j+i;
				ptrPoint1Z=coordZ+NI*NJ*k+NI*j+i;/*point1 is the point with minimum indices in a cell i.e (i,j,k)*/
				FindCellCentre(ptrSize,ptrPoint1X,ptrPoint1Y,ptrPoint1Z,centroid);
				sum1= sum1+v2*centroid[1]*(*(ptrVolume+index));
			
			}
			
			else 
			{	
				sum2= sum2+v2*(centroid[1])*(*(ptrVolume+index));
				
				
			}
			sum3=sum3+v2*(*(ptrVolume+index));
			
			index++;
		}
		//end
		sumzone=sumzone+sum1+sum2;
	}
	
	C_g=sumzone/sum3;	
	if((CounterTstep)%1000==0)
		fprintf(CG,"\n%e\t%e\t%e\t%d",C_g,Uavg_d_Y ,sum3,(CounterTstep)/1000);
	fclose(CG);

	
	return 0;

}
