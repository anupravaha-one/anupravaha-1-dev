#include "APKCheckConservation.h"
/**
This function CheckConservation() is called from function solverSemiCoupled2d() or solverSemiCoupled(). This is used only if the disperse phase flow is choosen. This function is used for conservation check of mass, momentum and energy for gas disperse phase flow. 
*/
int CheckConservation(Root CgnsRoot,int CounterTstep)
{
 
	int numZones,nzone,*ptrSize,a,b,c;
	
	//double fluxsum,fluxsumM,fractal;
	int indexP, NI,NJ,NK;
	
	//double variable,vol,variableT;
	Zone *zone;
	FILE *conserve;
	
	nzone=CgnsRoot.CgnsBase->Number_of_zones;
	
	double *coordX,*coordY,*coordZ,*ptrPoint1X,*ptrPoint1Y,*ptrPoint1Z;
	double centroid[3],centroidfirst[3],tau,tauclus,press,centroidprev;
	int numBc,zBc,index;
	//double *ptrPoint1XX;
	
	//FILE *front,*frontVF,*VelTest,*countWrite;
	//char frontnameT[50],VF_FOU[50],T_FOU[50];
	//char frontnameVF[50];
	//char tempVF[50],tempT[50];
	

	double *domainVolume,*F,volume_flux,vol_E,vol_P,FaceVolFrac,*FD,Moment_fluxDisperseM,Moment_fluxDisperse,Moment_fluxCarrier,FaceVolFracC,volume_fluxC,*FC;
	int indexE;
	
	//conserve=fopen("print/ConserveCheck.dat","w");	
	
	double rhodisperse;
	int flagrhodisperse;
	

	flagrhodisperse=CgnsRoot.CgnsBase->BaseParameters->flagRhoDisperse;
	
	conserve=fopen("SOLUTIONS/Conserve_Test.dat","a");
	fprintf(conserve,"\nTime step = %d\n",CounterTstep);
	fprintf(conserve,"\n##########################################################\n");	tau=0.0;tauclus=0.0;
	for(numZones=0;numZones<nzone;numZones++)
		{
			zone=(CgnsRoot.CgnsBase->CgnsZone+numZones);
			
			coordX=(zone->CgnsGrid)->X;
			coordY=(zone->CgnsGrid)->Y; 
			coordZ=(zone->CgnsGrid)->Z; 
			
			ptrSize=zone->Size;
			NI=*(ptrSize+0);
			NJ=*(ptrSize+1);
			NK=*(ptrSize+2);
			c=0;
			b=0;
			for(a=0;a<NI-1;a++)
			{
			
			ptrPoint1X=coordX+NI*NJ*c+NI*b+a;
			ptrPoint1Y=coordY+NI*NJ*c+NI*b+a;
			ptrPoint1Z=coordZ+NI*NJ*c+NI*b+a;
			FindCellCentre(ptrSize,ptrPoint1X,ptrPoint1Y,ptrPoint1Z,centroid);
			
			indexP=(NI-1)*(NJ-1)*c+(NI-1)*b+a;
			
			tau =tau + (6.0/180.0)*(*(zone->Flowsol->Tcoeff->EffViscoeff+indexP))*(*(zone->Flowsol->u+indexP)/centroid[1]);
			tauclus =tauclus + ((*(ptrPoint1X+1))-(*(ptrPoint1X)))*(*(zone->Flowsol->Rho+indexP*0))*0.3*(*(zone->Flowsol->Scal[0]+indexP));
			//tauclus =tauclus + ((*(ptrPoint1X+1))-(*(ptrPoint1X)))*(*(zone->Flowsol->Viscoeff+0))*(*(zone->Flowsol->u+indexP)/centroid[1]);
			
			fprintf(conserve,"\n %d %4.8e V drag %4.8e %4.8e",a,centroid[0],2*tau,2*tauclus);
			}
		}

		press = 0.0;
for(numZones=0;numZones<nzone;numZones++)
		{
			zone=(CgnsRoot.CgnsBase->CgnsZone+numZones);
			
			numBc=(zone)->Number_of_BC_Patches;
			
			coordX=(zone->CgnsGrid)->X;
			coordY=(zone->CgnsGrid)->Y; 
			coordZ=(zone->CgnsGrid)->Z; 
				
			ptrSize=zone->Size;
			NI=*(ptrSize+0);
			NJ=*(ptrSize+1);
			NK=*(ptrSize+2);
			for(zBc=0;zBc<numBc;zBc++)
			{
			if((((zone)->ZoneBc+zBc)->Direction==1)||(((zone)->ZoneBc+zBc)->Direction==2))
			{
	        	c=0;
			a=0;
			b=0;
			index= c*(NJ-1)+b;
			
			ptrPoint1X=coordX+NI*NJ*c+NI*b+a;
			ptrPoint1Y=coordY+NI*NJ*c+NI*b+a;
			ptrPoint1Z=coordZ+NI*NJ*c+NI*b+a;
			
			FindCellCentre(ptrSize,ptrPoint1X,ptrPoint1Y,ptrPoint1Z,centroidfirst);
			centroidprev = centroidfirst[1];
			indexP=(NI-1)*(NJ-1)*c+(NI-1)*b+a;
			press = *(((zone)->ZoneBc+zBc)->ZSlab->P+index)*(centroidfirst[1]);
			
			for(b=1;b<NJ-2;b++)
			{
			index= c*(NJ-1)+b;
			ptrPoint1X=coordX+NI*NJ*c+NI*b+a;
			ptrPoint1Y=coordY+NI*NJ*c+NI*b+a;
			ptrPoint1Z=coordZ+NI*NJ*c+NI*b+a;
			
			FindCellCentre(ptrSize,ptrPoint1X,ptrPoint1Y,ptrPoint1Z,centroid);
			
			indexP=(NI-1)*(NJ-1)*c+(NI-1)*b+a;
			
			press = press +  (*(((zone)->ZoneBc+zBc)->ZSlab->P+index)*(centroid[1]-centroidprev));
			centroidprev = centroid[1];
			}
			b=NJ-2;
			index= c*(NJ-1)+b;
			indexP=(NI-1)*(NJ-1)*c+(NI-1)*b+a;
			press =  press + *(((zone)->ZoneBc+zBc)->ZSlab->P+index)*(0.1-centroidprev);
			}
		   	}
		fprintf(conserve,"\n total pressure force at inlet = %12.12lf \n",press);
			
		        
			for(a=1;a<NI-2;a++)
			{
			press = 0.0;
			c=0;
			b=0;
			index= c*(NJ-1)+b;
			
			ptrPoint1X=coordX+NI*NJ*c+NI*b+a;
			ptrPoint1Y=coordY+NI*NJ*c+NI*b+a;
			ptrPoint1Z=coordZ+NI*NJ*c+NI*b+a;
			
			FindCellCentre(ptrSize,ptrPoint1X,ptrPoint1Y,ptrPoint1Z,centroidfirst);
			centroidprev = centroidfirst[1];
			indexP=(NI-1)*(NJ-1)*c+(NI-1)*b+a;
			press = (*(zone->Flowsol->P+indexP))*(centroidfirst[1]);
			
			for(b=1;b<NJ-2;b++)
			{
			index= c*(NJ-1)+b;
			ptrPoint1X=coordX+NI*NJ*c+NI*b+a;
			ptrPoint1Y=coordY+NI*NJ*c+NI*b+a;
			ptrPoint1Z=coordZ+NI*NJ*c+NI*b+a;
			
			FindCellCentre(ptrSize,ptrPoint1X,ptrPoint1Y,ptrPoint1Z,centroid);
			
			indexP=(NI-1)*(NJ-1)*c+(NI-1)*b+a;
			
			press = press +  ((*(zone->Flowsol->P+indexP))*(centroid[1]-centroidprev));
			centroidprev = centroid[1];
			}
			b=NJ-2;
			index= c*(NJ-1)+b;
			indexP=(NI-1)*(NJ-1)*c+(NI-1)*b+a;
			press =  press + (*(zone->Flowsol->P+indexP))*(0.1-centroidprev);
		fprintf(conserve,"\npressure at  %d %4.8e",a,press);
			}
		   	
		
		}
		
		for(numZones=0;numZones<nzone;numZones++)
		{
			zone=(CgnsRoot.CgnsBase->CgnsZone+numZones);
			domainVolume=MakeDomainArray(zone,5);
			F = zone->Flowsol->MassfluxDisperse;
	                FC = zone->Flowsol->Massflux;
			ptrSize=zone->Size;
	
		
			NI=*(ptrSize+0);
			NJ=*(ptrSize+1);
			NK=*(ptrSize+2);
		     c=0;
		     for(a=0;a<(NI-1);a=a+1)
		      {
			volume_flux=0.0;
			volume_fluxC=0.0; 
			for(b=0;b<NJ-1;b++)
			{
			if(a<(NI-2))
			{
			indexP=(NI-1)*(NJ-1)*c+(NI-1)*b+a;
			indexE=(NJ-1)*(NI-1)*c+(NI-1)*b+(a+1);
			rhodisperse = *(zone->Flowsol->RhoDisperse+indexP*flagrhodisperse);
			vol_E=*(domainVolume+indexE);
			vol_P=*(domainVolume+indexP);
			
			FaceVolFrac=((*(zone->Flowsol->VarthetaDisperse+indexP))*vol_E)+((*(zone->Flowsol->VarthetaDisperse+indexE))*vol_P);
			
			FaceVolFracC=((*(zone->Flowsol->Rho+indexP*(zone->Flowsol->flagRho)))*vol_E)+((*(zone->Flowsol->Rho+indexE*(zone->Flowsol->flagRho)))*vol_P);
			
			FaceVolFracC = FaceVolFracC/(vol_P+vol_E);
			
			FaceVolFrac = FaceVolFrac/(vol_P+vol_E);
			
			volume_flux = volume_flux + FaceVolFrac*(*(F+((NJ-1)*(NI-1)*c+(NI-1)*b+a)*6+1))*rhodisperse;
			
			volume_fluxC = volume_fluxC + FaceVolFracC*(*(FC+((NJ-1)*(NI-1)*c+(NI-1)*b+a)*6+1));
			}
			else
			{
			indexP=(NI-1)*(NJ-1)*c+(NI-1)*b+a;
			
			rhodisperse = *(zone->Flowsol->RhoDisperse+indexP*flagrhodisperse);
			
			vol_P=*(domainVolume+indexP);
			
			FaceVolFrac=(*(zone->Flowsol->VarthetaDisperse+indexP));
			
			FaceVolFracC=(*(zone->Flowsol->Rho+indexP*(zone->Flowsol->flagRho)));
			
			volume_flux = volume_flux + FaceVolFrac*(*(F+((NJ-1)*(NI-1)*c+(NI-1)*b+a)*6+1))*rhodisperse;
			
			volume_fluxC = volume_fluxC + FaceVolFracC*(*(FC+((NJ-1)*(NI-1)*c+(NI-1)*b+a)*6+1));
			}
			
			}
		      fprintf(conserve,"\n Mass Flux at I = %d Disperse %12.12lf Carrier %12.12lf \n",a,volume_flux,volume_fluxC);
		      }			   	
		}
////////////////////////////////////////////////////////////////////////////////
for(numZones=0;numZones<nzone;numZones++)
		{
			zone=(CgnsRoot.CgnsBase->CgnsZone+numZones);
			FD = zone->Flowsol->MassfluxDisperse;
	                F = zone->Flowsol->Massflux;
	                
		       
			
			ptrSize=zone->Size;
			NI=*(ptrSize+0);
			NJ=*(ptrSize+1);
			NK=*(ptrSize+2);
		     	
			Moment_fluxDisperse=0.0;
			Moment_fluxCarrier=0.0;
			Moment_fluxDisperseM=0.0;
			
		     c=0;
		     for(a=1;a<(NI-1);a=a+1)
		      {
			Moment_fluxDisperse=0.0;
			Moment_fluxCarrier=0.0;
			Moment_fluxDisperseM=0.0;
			for(b=0;b<NJ-1;b++)
			{
			indexP=(NI-1)*(NJ-1)*c+(NI-1)*b+a;
			
			rhodisperse = *(zone->Flowsol->RhoDisperse+indexP*flagrhodisperse);
			
			Moment_fluxDisperse = Moment_fluxDisperse + (*(FD+((NJ-1)*(NI-1)*c+(NI-1)*b+a)*6+1))*(*(zone->Flowsol->VarthetaDisperse+indexP))*(*(zone->Flowsol->U_Disperse+indexP))*rhodisperse;
			
			Moment_fluxDisperseM = Moment_fluxDisperseM + (*(FD+((NJ-1)*(NI-1)*c+(NI-1)*b+a)*6+1))*(*(zone->Flowsol->VarthetaDisperse+indexP))*((*(zone->Flowsol->MomentumDispersex+indexP))/((*(zone->Flowsol->VarthetaDisperse+indexP))*rhodisperse))*rhodisperse;
			
			Moment_fluxCarrier = Moment_fluxCarrier + (*(F+((NJ-1)*(NI-1)*c+(NI-1)*b+a)*6+1))*(*(zone->Flowsol->u+indexP))*(*(zone->Flowsol->Rho+indexP*(zone->Flowsol->flagRho)));
			
	
			}
                        fprintf(conserve,"\n%d MoDisp=%12.12lf  MoDisM=%12.12lf   MoCar=%12.12lf ",a,Moment_fluxDisperse,Moment_fluxDisperseM,Moment_fluxCarrier);
			}
	       }
// for(numZones=0;numZones<nzone;numZones++)
// 		{
// 			
// 			zone=(CgnsRoot.CgnsBase->CgnsZone+numZones);
// 			
// 			ptrSize=zone->Size;
// 			NI=*(ptrSize+0);
// 			NJ=*(ptrSize+1);
// 			NK=*(ptrSize+2);
// 		     	
// 			for(c=0;c<NK-1;c++)
// 			for(b=0;b<NJ-1;b++)
// 			for(a=0;a<NI-1;a++)
// 			 {
// 			 indexP=(NJ-1)*(NI-1)*c+(NI-1)*b+a;
// 			 variable = *(zone->Flowsol->VarthetaDisperse+indexP);
// 			
// 		         Moment_fluxDisperse=0.0;
// 			 Moment_fluxCarrier=0.0;
// 			 Moment_fluxDisperseM=0.0;
// 			 
// 			
// 			Moment_fluxDisperse = Moment_fluxDisperse + (*(zone->Flowsol->DropFuelMassFrac+indexP));
// 			
// 			Moment_fluxDisperseM = Moment_fluxDisperseM + (*(FD+((NJ-1)*(NI-1)*c+(NI-1)*b+a)*6+1))*(*(zone->Flowsol->VarthetaDisperse+indexP))*((*(zone->Flowsol->MomentumDispersex+indexP))/((*(zone->Flowsol->VarthetaDisperse+indexP))*rhodisperse))*rhodisperse;
// 			
// 			Moment_fluxCarrier = Moment_fluxCarrier + (*(F+((NJ-1)*(NI-1)*c+(NI-1)*b+a)*6+1))*(*(zone->Flowsol->u+indexP))*(*(zone->Flowsol->Rho+indexP*(zone->Flowsol->flagRho)));
// 			
// 	
// 			}
//                         fprintf(conserve,"\n%d MoDisp=%12.12lf  MoDisM=%12.12lf   MoCar=%12.12lf ",a,Moment_fluxDisperse,Moment_fluxDisperseM,Moment_fluxCarrier);
// 			
// 		}
fclose(conserve);

return 0;
}



