#include "APKFindgradient.h"
/**
The function Findgradient() calls the functions CalcGradViscoeff() and CalcGradVar() to compute the gradients
of viscosity and that of velocity components respectively, which arise as source terms in the momentum equations.
Please refer Eqn(5.17) - Eqn(5.19) in Project Report No 1, Dec 2005.
 */
double Findgradient(Base *base,Zone *zone,int index[3],int variableTag)
{
	double source;
	double *gradVariable1Slab[6],*gradVariable2Slab[6],*gradVariable3Slab[6];
	double *variable1,*variable2,*variable3;
	double DV[3],DW[3],DMu[3];
	//Turbulent gas-Particle flow
	double *U_g,*V_g,*W_g,*SlabU_g[6],*SlabV_g[6],*SlabW_g[6];
	if(base->flagMPFV==0)
	{
		if(base->flagTurb == 1)
		{
			if((variableTag==62)||(variableTag==63)||(variableTag==64))//disperse phase momentum
			{
				variable3 = zone->Flowsol->Tcoeff->DummyViscoeff_disperse;//
				MakeSlabPtrArray(zone,gradVariable3Slab,726);
			}
			else//Single phase turbulence or gas phase in case of gas-particle flow
			{
				variable3 = zone->Flowsol->Tcoeff->EffViscoeff;//
				MakeSlabPtrArray(zone,gradVariable3Slab,714);
			}
		}
		else
		{
			variable3 = zone->Flowsol->Viscoeff;
			MakeSlabPtrArray(zone,gradVariable3Slab,711);
		}
// 		if(variableTag==62)
// 		{
// 		printf("The eff viscosity for variable is %lf and variable tag is %d [Zone %s]\n",*variable3,variableTag,zone->zonename);
// 		printf("The value in the 0th slab is %e\n",*gradVariable3Slab[0]);
// 		printf("The value in the 1st slab is %e\n",*gradVariable3Slab[1]);
// 		printf("The value in the 2nd slab is %e\n",*gradVariable3Slab[2]);
// 		printf("The value in the 3rd slab is %e\n",*gradVariable3Slab[3]);
// 		printf("The value in the 4th slab is %e\n",*gradVariable3Slab[4]);
// 		printf("The value in the 5th slab is %e\n",*gradVariable3Slab[5]);
// 		}
	}
	else if(base->flagMPFV==1)//Added by chandrasekhar --- Turbulent MPFV
	{
		if((base->flagTurb == 1)&((variableTag==31)||(variableTag==32)||(variableTag==33)||(variableTag==3100)||(variableTag==3200)||(variableTag==3300)))
		{
			variable3 = zone->Flowsol->Tcoeff->EffViscoeff;
			MakeSlabPtrArray(zone,gradVariable3Slab,714);
		}
		else if((base->flagTurb == 1)&((variableTag==44)||(variableTag==45)||(variableTag==46)))
		{
			variable3 = zone->Flowsol->Tcoeff->DummyeffViscoeffofphase[1];
			MakeSlabPtrArray(zone,gradVariable3Slab,731);
		}
// 		if(variableTag==62)
// 		{
// 		printf("The eff viscosity for variable is %lf and variable tag is %d [Zone %s]\n",*variable3,variableTag,zone->zonename);
// 		printf("The value in the 0th slab is %lf\n",*gradVariable3Slab[0]);
// 		printf("The value in the 1st slab is %lf\n",*gradVariable3Slab[1]);
// 		printf("The value in the 2nd slab is %lf\n",*gradVariable3Slab[2]);
// 		printf("The value in the 3rd slab is %lf\n",*gradVariable3Slab[3]);
// 		printf("The value in the 4th slab is %lf\n",*gradVariable3Slab[4]);
// 		printf("The value in the 5th slab is %lf\n",*gradVariable3Slab[5]);
// 		}
	}
	CalcGradViscoeff(base,zone,variable3,gradVariable3Slab,index,DMu);
	if(base->flagMPFV==0)
	{
		if(base->BaseParameters->flag2d==0)
		{
			switch(variableTag)
			{
				case 1:
				{
					variable1 = zone->PrevFlowsol->v;
					variable2 = zone->PrevFlowsol->w;
					MakeSlabPtrArray(zone,gradVariable1Slab,20);
					MakeSlabPtrArray(zone,gradVariable2Slab,30);
	// 				variable1 = zone->Flowsol->v;
	// 				variable2 = zone->Flowsol->w;
	// 				MakeSlabPtrArray(zone,gradVariable1Slab,2);
	// 				MakeSlabPtrArray(zone,gradVariable2Slab,3);
					CalcGradVar(base,zone,variable1,gradVariable1Slab,index,DV);
					CalcGradVar(base,zone,variable2,gradVariable2Slab,index,DW);
					source=DMu[1]*DV[0]+DMu[2]*DW[0]-DMu[0]*DV[1]-DMu[0]*DW[2];
				}
				break;
		
				case 100:
				{
					variable1 = zone->PrevFlowsol->v;
					variable2 = zone->PrevFlowsol->w;
					MakeSlabPtrArray(zone,gradVariable1Slab,20);/* Actually PrevFlowsol Slab value of v has to be passed but current Slab-> v contains the same data*/
					MakeSlabPtrArray(zone,gradVariable2Slab,30);
					CalcGradVar(base,zone,variable1,gradVariable1Slab,index,DV);
					CalcGradVar(base,zone,variable2,gradVariable2Slab,index,DW);
					source=DMu[1]*DV[0]+DMu[2]*DW[0]-DMu[0]*DV[1]-DMu[0]*DW[2];
				}
				break;
		
				case 2:
				{
					variable1 = zone->PrevFlowsol->u;
					variable2 = zone->PrevFlowsol->w;
					MakeSlabPtrArray(zone,gradVariable1Slab,10);
					MakeSlabPtrArray(zone,gradVariable2Slab,30);
	// 				variable1 = zone->Flowsol->u;
	// 				variable2 = zone->Flowsol->w;
	// 				MakeSlabPtrArray(zone,gradVariable1Slab,1);
	// 				MakeSlabPtrArray(zone,gradVariable2Slab,3);
					CalcGradVar(base,zone,variable1,gradVariable1Slab,index,DV);
					CalcGradVar(base,zone,variable2,gradVariable2Slab,index,DW);
					source=DMu[0]*DV[1]+DMu[2]*DW[1]-DMu[1]*(DV[0]+DW[2]);
		
				}
				break;
		
				case 200:
				{
					variable1 = zone->PrevFlowsol->u;
					variable2 = zone->PrevFlowsol->w;
					MakeSlabPtrArray(zone,gradVariable1Slab,10);
					MakeSlabPtrArray(zone,gradVariable2Slab,30);
					CalcGradVar(base,zone,variable1,gradVariable1Slab,index,DV);
					CalcGradVar(base,zone,variable2,gradVariable2Slab,index,DW);
					source=DMu[0]*DV[1]+DMu[2]*DW[1]-DMu[1]*(DV[0]+DW[2]);
		
				}
				break;
		
				case 3:
				{
					variable1 = zone->PrevFlowsol->u;
					variable2 = zone->PrevFlowsol->v;
					MakeSlabPtrArray(zone,gradVariable1Slab,10);
					MakeSlabPtrArray(zone,gradVariable2Slab,20);
	// 				variable1 = zone->Flowsol->u;
	// 				variable2 = zone->Flowsol->v;
	// 				MakeSlabPtrArray(zone,gradVariable1Slab,1);
	// 				MakeSlabPtrArray(zone,gradVariable2Slab,2);
					CalcGradVar(base,zone,variable1,gradVariable1Slab,index,DV);
					CalcGradVar(base,zone,variable2,gradVariable2Slab,index,DW);
					source=DMu[0]*DV[2]+DMu[1]*DW[2]-DMu[2]*(DV[0]+DW[1]);
				}
				break;
		
				case 300:
				{
					variable1 = zone->PrevFlowsol->u;
					variable2 = zone->PrevFlowsol->v;
					MakeSlabPtrArray(zone,gradVariable1Slab,10);
					MakeSlabPtrArray(zone,gradVariable2Slab,20);
					CalcGradVar(base,zone,variable1,gradVariable1Slab,index,DV);
					CalcGradVar(base,zone,variable2,gradVariable2Slab,index,DW);
					source=DMu[0]*DV[2]+DMu[1]*DW[2]-DMu[2]*(DV[0]+DW[1]);
				}
				break;
				
				case 62:
				{
					variable1 = zone->PrevFlowsol->MomentumDispersey;
					variable2 = zone->PrevFlowsol->MomentumDispersez;
					MakeSlabPtrArray(zone,gradVariable1Slab,630);
					MakeSlabPtrArray(zone,gradVariable2Slab,640);
					
					V_g=zone->PrevFlowsol->v;
					MakeSlabPtrArray(zone,SlabV_g,20);
					W_g=zone->PrevFlowsol->w;
					MakeSlabPtrArray(zone,SlabW_g,30);
					
					CalcGradVar_Disperse(base,zone,variable1,gradVariable1Slab,index,DV,V_g,SlabV_g);
					CalcGradVar_Disperse(base,zone,variable2,gradVariable2Slab,index,DW,W_g,SlabW_g);
					source=DMu[1]*DV[0]+DMu[2]*DW[0]-DMu[0]*DV[1]-DMu[0]*DW[2];
		
				}
				
				case 63:
				{
					variable1 = zone->PrevFlowsol->MomentumDispersex;
					variable2 = zone->PrevFlowsol->MomentumDispersez;
					MakeSlabPtrArray(zone,gradVariable1Slab,620);
					MakeSlabPtrArray(zone,gradVariable2Slab,640);
					
					U_g=zone->PrevFlowsol->u;
					MakeSlabPtrArray(zone,SlabU_g,10);
					W_g=zone->PrevFlowsol->w;
					MakeSlabPtrArray(zone,SlabW_g,30);
					
					CalcGradVar_Disperse(base,zone,variable1,gradVariable1Slab,index,DV,U_g,SlabU_g);
					CalcGradVar_Disperse(base,zone,variable2,gradVariable2Slab,index,DW,W_g,SlabW_g);
					source=DMu[0]*DV[1]+DMu[2]*DW[1]-DMu[1]*(DV[0]+DW[2]);;
		
				}
				
				case 64:
				{
					variable1 = zone->PrevFlowsol->MomentumDispersex;
					variable2 = zone->PrevFlowsol->MomentumDispersey;
					MakeSlabPtrArray(zone,gradVariable1Slab,620);
					MakeSlabPtrArray(zone,gradVariable2Slab,630);
					
					U_g=zone->PrevFlowsol->u;
					MakeSlabPtrArray(zone,SlabU_g,10);
					V_g=zone->PrevFlowsol->v;
					MakeSlabPtrArray(zone,SlabV_g,20);
					
					CalcGradVar_Disperse(base,zone,variable1,gradVariable1Slab,index,DV,U_g,SlabU_g);
					CalcGradVar_Disperse(base,zone,variable2,gradVariable2Slab,index,DW,V_g,SlabV_g);
					source=DMu[0]*DV[2]+DMu[1]*DW[2]-DMu[2]*(DV[0]+DW[1]);
		
				}
			}
		}
		else
		{
			switch(variableTag)
			{
				case 1:
				{
					variable1 = zone->PrevFlowsol->v;
					
					MakeSlabPtrArray(zone,gradVariable1Slab,20);
					
					CalcGradVar(base,zone,variable1,gradVariable1Slab,index,DV);
	
					source=DMu[1]*DV[0]-DMu[0]*DV[1];
				}
				break;
		
				case 100:
				{
					variable1 = zone->PrevFlowsol->v;
	
					MakeSlabPtrArray(zone,gradVariable1Slab,20);
					
					/// Actually PrevFlowsol Slab value of v has to be passed but current Slab-> v contains the same data
	
					CalcGradVar(base,zone,variable1,gradVariable1Slab,index,DV);
	
					source=DMu[1]*DV[0]-DMu[0]*DV[1];
				}
				break;
		
				case 2:
				{
					variable1 = zone->PrevFlowsol->u;
	
					MakeSlabPtrArray(zone,gradVariable1Slab,10);
	
					CalcGradVar(base,zone,variable1,gradVariable1Slab,index,DV);
	
					source=DMu[0]*DV[1]-DMu[1]*DV[0];
		
				}
				break;
		
				case 200:
				{
					variable1 = zone->PrevFlowsol->u;
	
					MakeSlabPtrArray(zone,gradVariable1Slab,10);
	
					CalcGradVar(base,zone,variable1,gradVariable1Slab,index,DV);
	
					source=DMu[0]*DV[1]-DMu[1]*DV[0];
		
				}
				break;
				
				case 62:
				{
					variable1 = zone->PrevFlowsol->MomentumDispersey;
					
					MakeSlabPtrArray(zone,gradVariable1Slab,630);
					
					V_g=zone->PrevFlowsol->v;
					
					MakeSlabPtrArray(zone,SlabV_g,20);
					
					CalcGradVar_Disperse(base,zone,variable1,gradVariable1Slab,index,DV,V_g,SlabV_g);
	
					source=DMu[1]*DV[0]-DMu[0]*DV[1];
				}
				break;
				
				case 63:
				{
					variable1 = zone->PrevFlowsol->MomentumDispersex;
					
					MakeSlabPtrArray(zone,gradVariable1Slab,620);
					
					U_g=zone->PrevFlowsol->u;
					
					MakeSlabPtrArray(zone,SlabU_g,10);
					
					CalcGradVar_Disperse(base,zone,variable1,gradVariable1Slab,index,DV,U_g,SlabU_g);
	
					source=DMu[0]*DV[1]-DMu[1]*DV[0];
				}
				break;
			}
		}
	}	
	if(base->flagMPFV==1)//Added by chandrasekhar --- Turbulent MPFV
	{
		if(base->BaseParameters->flag2d==0)
		{
			switch(variableTag)
			{
				case 31:
				case 44:
				{
					variable1 = zone->PrevFlowsol->Momentumy;
					variable2 = zone->PrevFlowsol->Momentumz;
					MakeSlabPtrArray(zone,gradVariable1Slab,320);
					MakeSlabPtrArray(zone,gradVariable2Slab,330);
	// 				variable1 = zone->Flowsol->v;
	// 				variable2 = zone->Flowsol->w;
	// 				MakeSlabPtrArray(zone,gradVariable1Slab,2);
	// 				MakeSlabPtrArray(zone,gradVariable2Slab,3);
					CalcGradVar(base,zone,variable1,gradVariable1Slab,index,DV);
					CalcGradVar(base,zone,variable2,gradVariable2Slab,index,DW);
					source=DMu[1]*DV[0]+DMu[2]*DW[0]-DMu[0]*DV[1]-DMu[0]*DW[2];
				}
				break;
		
				case 3100:
				{
					variable1 = zone->PrevFlowsol->Momentumy;
					variable2 = zone->PrevFlowsol->Momentumz;
					MakeSlabPtrArray(zone,gradVariable1Slab,320);/* Actually PrevFlowsol Slab value of v has to be passed but current Slab-> v contains the same data*/
					MakeSlabPtrArray(zone,gradVariable2Slab,330);
					CalcGradVar(base,zone,variable1,gradVariable1Slab,index,DV);
					CalcGradVar(base,zone,variable2,gradVariable2Slab,index,DW);
					source=DMu[1]*DV[0]+DMu[2]*DW[0]-DMu[0]*DV[1]-DMu[0]*DW[2];
				}
				break;
		
				case 32:
				case 45:
				{
					variable1 = zone->PrevFlowsol->Momentumx;
					variable2 = zone->PrevFlowsol->Momentumz;
					MakeSlabPtrArray(zone,gradVariable1Slab,310);
					MakeSlabPtrArray(zone,gradVariable2Slab,330);
	// 				variable1 = zone->Flowsol->u;
	// 				variable2 = zone->Flowsol->w;
	// 				MakeSlabPtrArray(zone,gradVariable1Slab,1);
	// 				MakeSlabPtrArray(zone,gradVariable2Slab,3);
					CalcGradVar(base,zone,variable1,gradVariable1Slab,index,DV);
					CalcGradVar(base,zone,variable2,gradVariable2Slab,index,DW);
					source=DMu[0]*DV[1]+DMu[2]*DW[1]-DMu[1]*(DV[0]+DW[2]);
		
				}
				break;
		
				case 3200:
				{
					variable1 = zone->PrevFlowsol->Momentumx;
					variable2 = zone->PrevFlowsol->Momentumz;
					MakeSlabPtrArray(zone,gradVariable1Slab,310);
					MakeSlabPtrArray(zone,gradVariable2Slab,330);
					CalcGradVar(base,zone,variable1,gradVariable1Slab,index,DV);
					CalcGradVar(base,zone,variable2,gradVariable2Slab,index,DW);
					source=DMu[0]*DV[1]+DMu[2]*DW[1]-DMu[1]*(DV[0]+DW[2]);
		
				}
				break;
		
				case 33:
				case 46:
				{
					variable1 = zone->PrevFlowsol->Momentumx;
					variable2 = zone->PrevFlowsol->Momentumy;
					MakeSlabPtrArray(zone,gradVariable1Slab,310);
					MakeSlabPtrArray(zone,gradVariable2Slab,320);
	// 				variable1 = zone->Flowsol->u;
	// 				variable2 = zone->Flowsol->v;
	// 				MakeSlabPtrArray(zone,gradVariable1Slab,1);
	// 				MakeSlabPtrArray(zone,gradVariable2Slab,2);
					CalcGradVar(base,zone,variable1,gradVariable1Slab,index,DV);
					CalcGradVar(base,zone,variable2,gradVariable2Slab,index,DW);
					source=DMu[0]*DV[2]+DMu[1]*DW[2]-DMu[2]*(DV[0]+DW[1]);
				}
				break;
		
				case 3300:
				{
					variable1 = zone->PrevFlowsol->Momentumx;
					variable2 = zone->PrevFlowsol->Momentumy;
					MakeSlabPtrArray(zone,gradVariable1Slab,310);
					MakeSlabPtrArray(zone,gradVariable2Slab,320);
					CalcGradVar(base,zone,variable1,gradVariable1Slab,index,DV);
					CalcGradVar(base,zone,variable2,gradVariable2Slab,index,DW);
					source=DMu[0]*DV[2]+DMu[1]*DW[2]-DMu[2]*(DV[0]+DW[1]);
				}
				break;
				
				default:
					printf("Tag not found in Findgradient and exiting\n");
					exit(0);
			}
		}
		else
		{
			switch(variableTag)
			{
				case 31:
				case 44: 
				{
					variable1 = zone->PrevFlowsol->Momentumy;
					
					MakeSlabPtrArray(zone,gradVariable1Slab,320);
					
					CalcGradVar(base,zone,variable1,gradVariable1Slab,index,DV);
	
					source=DMu[1]*DV[0]-DMu[0]*DV[1];
				}
				break;
		
				case 3100:
				{
					variable1 = zone->PrevFlowsol->Momentumy;
	
					MakeSlabPtrArray(zone,gradVariable1Slab,320);
					
					/// Actually PrevFlowsol Slab value of v has to be passed but current Slab-> v contains the same data
					
					CalcGradVar(base,zone,variable1,gradVariable1Slab,index,DV);
	
					source=DMu[1]*DV[0]-DMu[0]*DV[1];
				}
				break;
		
				case 32:
				case 45:
				{
					variable1 = zone->PrevFlowsol->Momentumx;
	
					MakeSlabPtrArray(zone,gradVariable1Slab,310);
	
					CalcGradVar(base,zone,variable1,gradVariable1Slab,index,DV);
	
					source=DMu[0]*DV[1]-DMu[1]*DV[0];//printf("The source is %ld\n",source);
		
				}
				break;
		
				case 3200:
				{
					variable1 = zone->PrevFlowsol->Momentumx;
	
					MakeSlabPtrArray(zone,gradVariable1Slab,310);
					//printf("The value is gradVariable1Slab[0] is %lf",*gradVariable1Slab[0]);
					CalcGradVar(base,zone,variable1,gradVariable1Slab,index,DV);
	
					source=DMu[0]*DV[1]-DMu[1]*DV[0];
		
				}
				break;
				
				default:
					printf("Tag not found in Findgradient and exiting\n");
					exit(0);
			}
 			//printf("\n Vartag==%d",variableTag);
// 			if(variableTag==3100)
// 			{
// 			printf("The value for variable is %lf and variable tag is %d [Zone %s]\n",*variable1,variableTag,zone->zonename);
// 			printf("The value in the 0th slab is %lf\n",*gradVariable1Slab[0]);
// 			printf("The value in the 1st slab is %lf\n",*gradVariable1Slab[1]);
// 			printf("The value in the 2nd slab is %lf\n",*gradVariable1Slab[2]);
// 			printf("The value in the 3rd slab is %lf\n",*gradVariable1Slab[3]);
// 			printf("The value in the 4th slab is %lf\n",*gradVariable1Slab[4]);
// 			printf("The value in the 5th slab is %lf\n",*gradVariable1Slab[5]);
// 			}
		}
	}
	return (source);
}


