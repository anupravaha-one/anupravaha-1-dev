#include "APKCalcEvapoConst.h"
#define PI 3.141592654
/**
This function CalcEvapoConst() is called in the beginning of the time step from fuction solverSemiCoupled2d() or solverSemiCoupled() if droplet phase flow is choosen. This fuction computes and returns the value of evaporation constant. This also computes the values of two source terms which are the source terms of energy equations and then stores them in appropriate variables.   
*/
double CalcEvapoConst(Base *base)
{
	double dia;
        
	int cntZ,NI,NJ,NK,numZones,i,j,k;
	int *ptrSize,flag,flagrho,flagmassfrac,flagDen,flagMu;
 	Zone *zone;
 	int totalcells,index;/*,flagrhoDisperse=0;*/
 	numZones= base->Number_of_zones;
 	
	
        double BM,BT,Tref,YFref,YAref,Cpg,Cpa,Cpv,Cpf,Kg,Ka,Kv;
	double Tinf,Tsurf,Yfs,LembdaStT,LembdaStM,Pfs,P;
	double Tr,L,RhoF;
	
	double mdotv,Yfin;
	
	double Q,QE;
		
	//shash
	double EPSILON;
	EPSILON = base->BaseParameters->EPSILONDrop;	

	double A,B,C;
	double Mdotv,nmdotv;
	double Red,temp,mu,rho,modRelV,u,ud,v,vd,w,wd;
	double Prg, coeff, factor;
	
	for(cntZ=0;cntZ<numZones;cntZ++)
	{
		zone = (base->CgnsZone+cntZ);
		ptrSize = zone->Size;
		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);
		totalcells = (NI-1)*(NJ-1)*(NK-1);
		
		flagrho=base->BaseParameters->flagRhoDisperse;
		flag=base->BaseParameters->flagVariableDia;
                flagmassfrac = base->BaseParameters->flagDropFuelMassFrac;
		flagDen = zone->Flowsol->flagRho;
		flagMu=zone->Flowsol->flagViscoeff;

		for(k=0;k<NK-1;k++)
		for(j=0;j<NJ-1;j++)
		for(i=0;i<NI-1;i++)
		{	
			index  = (NI-1)*(NJ-1)*k+(NI-1)*j+i;

//Calculating the temperature and pressure which determine evaporation rate
			if(zone->Flowsol->T!=NULL)
			{
				Tinf = *(zone->PrevFlowsol->T+index);
			}
			else
			{
				Tinf = 773.0;
			}
		
		
			P = *(zone->PrevFlowsol->P+index);
			P = (P+101325.0)/1000.0;//kPa
		
			if(flagmassfrac==1)
		   	{
		  		 Yfin = *(zone->Flowsol->DropFuelMassFrac+index);
		  	}
			else
		   	{
		  		 Yfin = 0.0;
		   	}	
		   
			Tsurf = *(zone->PrevFlowsol->Tdrop+index);

//Calculating density of droplet phase:

			RhoF = userDefinedDropRhoF(Tsurf);
		     
// Calculating the parameters required for finding out the evaporation rate:
			
			Yfs = userDefinedDropPfYf(Tsurf,P);			
			BM=(Yfs-Yfin)/(1-Yfs);
		
			YFref=(2.0/3.0)*Yfs+(1.0/3.0)*Yfin;	
			YAref=1-(2.0/3.0)*Yfs-(1.0/3.0)*Yfin; 
			Tref=Tsurf+(Tinf-Tsurf)/3.0;
		
			Cpa = userDefinedDropCpa(Tref);
		
			Cpv = userDefinedDropCpv(Tref);
		
			Cpg = YAref*Cpa+YFref*Cpv;
	
			Ka = userDefinedDropKa(Tref);
		
			Kv = userDefinedDropKv(Tref, Cpv); 
		
			Kg = YAref*Ka+YFref*Kv;

//Calculating the current time step diameter:

			A= base->BaseParameters->InitialDropletDiameter;
			
			B = base->BaseParameters->InitialVolumeFraction;

			C = *(zone->PrevFlowsol->VarthetaDisperse+index);

			*(zone->PrevFlowsol->ParticleDia+index*flag) = A*pow(C,1.0/3.0);
		      
			dia = *(zone->PrevFlowsol->ParticleDia+index*flag);


//Calculating the evaporation rate of a single droplet:
		
	      		mdotv=2.0*(22.0/7.0)*dia*(Kg/Cpg)*log(1+BM);

//Considering the convective effects:

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
		
			modRelV=(u-ud)*(u-ud)+(v-vd)*(v-vd)+(w-wd)*(w-wd);
		
			temp=pow(modRelV,0.5);
			modRelV=temp;
		
			Red=(rho*dia*modRelV)/mu;

			Prg = mu*Cpg/Kg;
			
			factor = pow(Red,0.5)*(pow(Prg,0.33));

			coeff = (1 + 0.3*factor);

//Calculating the total evaporation rate:

			nmdotv =12*(Kg/Cpg)*log(1+BM)*B*pow(C,1.0/3.0)/(pow(A,2))*coeff; // mdotv*(42.0/22.0)*(B)/pow(A,3);//

//Calculating parameters required for heat transfer due to evaporation:
			L = userDefinedDropL(Tsurf);
		
			Cpf = userDefinedDropCpf(Tsurf);
		
			BT=Cpg*(Tinf-Tsurf)/L;

//Calculating total heat transfer rates due to evaporation			
			Q = 12*Kg*((log(1+BM))/BM)*(Tinf-Tsurf)*B*pow(C,1.0/3.0)/(pow(A,2))*coeff;
			QE = nmdotv*L*coeff;	

//Calculating the steady state value of mass and heat transfer rates										
			LembdaStT= 8*Kg*log(1+BT)/(Cpg*RhoF);
			LembdaStM= 8*Kg*log(1+BM)/(Cpg*RhoF);

	        	Mdotv = LembdaStT*(0.7854)*(*(zone->Flowsol->RhoDisperse+index*flagrho)*dia);

		if((i==(int)NI/2)&&(j==(int)NJ/2)) //(int)NI/2
		{
			printf("\n\nFor i =%d, j=%d\n",i,j);
			printf("\nP=%lf \t Yfin = %lf \t Yfs = %lf\n",P,Yfin,Yfs);
			printf("\n**********************************************\n");
			printf("\nReference properties: \n");
			printf("\nTref = %lf \t YFref = %lf \t YAref = %lf\n", Tref,YFref,YAref);
			printf("\n**********************************************\n");
			printf("\nThermal properties:\n");
			printf("\nCpa = %lf \t Cpv = %lf \t Cpg = %lf \t L =%lf \n\n  BM = %lf \t BT = %lf \t LembdaStT = %e \t LembdaStM = %e\n Kg = %lf\n",Cpa,Cpv,Cpg,L,BM,BT,LembdaStT,LembdaStM,Kg);
			printf("\nSource Terms: \n");
			printf("dia = %e \t mdotv=%e \n\n nmdotv = %e \t Mdotv = %e \n\n nQ = %e \t nQE = %e \n",dia,mdotv,nmdotv,Mdotv,Q,QE);
			printf("\n**********************************************\n");
			printf("\nGas-Phase Terms:\n");
			printf("\nTg = %lf \t Ug = %lf \t rhog = %lf\n",Tinf,u,rho);
			printf("\nDroplet-Phase terms:\n");
			printf("\nTsurf = %lf \t Ud = %lf \t varthetad = %lf\n",Tsurf,ud,*(zone->PrevFlowsol->VarthetaDisperse+index));		

		}
		if(*(zone->PrevFlowsol->VarthetaDisperse+index)>EPSILON)
		{
			*(zone->PrevFlowsol->Evaporate+index) = nmdotv;
			*(zone->PrevFlowsol->QDroptoGas+index) = Q;
			*(zone->PrevFlowsol->QLatent+index) = QE;
	       	}
		else
		{
			*(zone->PrevFlowsol->Evaporate+index)=0.0;
		    	*(zone->PrevFlowsol->QDroptoGas+index) = 0.0;
			*(zone->PrevFlowsol->QLatent+index) = 0.0;
	       	}

		}
	
	}

	return 0;	
}



