#include "APKDrag.h"
/**
This function calculates the drag force for the case of disperse phase flow. The drag is calculated and stored in the data structure under the source term. This function is called once 
per time step
*/
double Drag(Base *base)
{
	double Re,Cd,u,v,w=0.0,ud,vd,wd=0.0,mu,rho,dia,Fd_x,Fd_y,Fd_z,FD_X,FD_Y,FD_Z;
        double modRelV,temp,beta,alpha=0.0;
	int cntZ,NI,NJ,NK,*ptrSize,numZones,i,j,k;
	Zone *zone;
	int flagDen,flagMu,flagDia;//Flags
	int totalcells,index;/*,flagrhoDisperse=0;*/
	double *tp,RhoD,modDrag,modDRAG;
	numZones= base->Number_of_zones;
	
	//shash
	double EPSILON;
	EPSILON = base->BaseParameters->EPSILONDrop;
	
	for(cntZ=0;cntZ<numZones;cntZ++)
	{
		zone = (base->CgnsZone+cntZ);
		ptrSize = zone->Size;
		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);
		tp=zone->t_p;
		totalcells = (NI-1)*(NJ-1)*(NK-1);
		flagMu=zone->Flowsol->flagViscoeff;
		flagDen=zone->Flowsol->flagRho;
		flagDia=base->BaseParameters->flagVariableDia;
		index=0;
		for(k=0;k<NK-1;k++)
		for(j=0;j<NJ-1;j++)
		for(i=0;i<NI-1;i++)
		{	
		u = *(zone->PrevFlowsol->u+index);
		v = *(zone->PrevFlowsol->v+index);
		ud = *(zone->PrevFlowsol->U_Disperse+index);
		vd = *(zone->PrevFlowsol->V_Disperse+index);
		
		if((zone->PrevFlowsol->w)!=NULL)
		    {
		    w = *(zone->PrevFlowsol->w+index);
		    wd = *(zone->PrevFlowsol->W_Disperse+index);
		    }
		
		mu= *(zone->PrevFlowsol->Viscoeff+index*flagMu);
		rho= *(zone->PrevFlowsol->Rho+index*flagDen);
		dia= *(zone->PrevFlowsol->ParticleDia+index*flagDia);
		RhoD=*(zone->PrevFlowsol->RhoDisperse);
		modRelV=(u-ud)*(u-ud)+(v-vd)*(v-vd)+(w-wd)*(w-wd);
		//printf("\nmu=%e$rho=%e$dia=%e$RhoD=%e$modRelV=%e",mu,rho,dia,RhoD,modRelV);
		temp=pow(modRelV,0.5);
		modRelV=temp;
		
		Re=(rho*dia*modRelV)/mu;
		//printf("\nRe=%lf",Re);
		if(base->constant==1)
			Cd=userDefinedCDDisperse(Re);
		else
			Cd = base->BaseParameters->ConstantCD;
	
// 		if(Re<=200000.0)
// 			{
// 			if(Re < 0.00000001)
// 			Cd=0;
// 			else if((Re > 0.00000001 )&&(Re<=1.0))
// 			Cd=24.0/Re;
// 			else if((Re > 1.0 )&&(Re<=1000.0))
// 			Cd=(24.0/Re)*(1+0.15*pow(Re,0.687));
// 			else if((Re>1000.0)&&(Re<=200000.0))
//         		Cd=0.44;
//        			}
// 		else
// 			Cd=0.09;
		
		beta=(0.75*Cd*rho*modRelV)/dia;
		alpha=0.0;
		
		Fd_x=beta*(u-ud);
		Fd_y=beta*(v-vd);
		Fd_z=beta*(w-wd);
		
		FD_X=Fd_x*(*(zone->PrevFlowsol->VarthetaDisperse+index));
		FD_Y=Fd_y*(*(zone->PrevFlowsol->VarthetaDisperse+index));
		FD_Z=Fd_z*(*(zone->PrevFlowsol->VarthetaDisperse+index));
		
		modDrag=(Fd_x*Fd_x)+(Fd_y*Fd_y)+(Fd_z*Fd_z);
		modDrag=sqrt(modDrag);
		
		modDRAG=(FD_X*FD_X)+(FD_Y*FD_Y)+(FD_Z*FD_Z);
		modDRAG=sqrt(modDRAG);
		
		*(zone->DragX+index)= Fd_x;
		*(zone->DragY+index)=Fd_y;
		if((zone->PrevFlowsol->w)!=NULL)
		    {
		    *(zone->DragZ+index)=Fd_z;
		    }       
		//Computation of particle response time and is taken as zero when volume fraction is less than the limiting value.   
		//And in general case it doesn't contain the volume fraction term as Fd_x,Fd_y,Fd_z doesn't contain the volume fraction term.
		if(base->flagTurb==1)
		{
			*(tp+index)=(RhoD*dia*dia)/(18.0*mu);
			//printf("\n *(tp+index)=%lf",*(tp+index));
// 			if(modDrag<EPSILON)
// 				*(tp+index)=0.0;
// 			else 
// 				*(tp+index)=RhoD*(*(zone->PrevFlowsol->VarthetaDisperse+index))*modRelV/modDRAG;
		}
		//printf("\nZone-%s$$*(tp+index)=%lf",zone->zonename,*(tp+index));	    								// 	        if(modDrag<EPSILON)
		index++;									// 			    *(tp+index)=0.0;
											// 		else 
											// 			    *(tp+index)=RhoD*(*(zone->PrevFlowsol->VarthetaDisperse+index))*modRelV/modDrag;

		}
		
        }
return 0;

}
