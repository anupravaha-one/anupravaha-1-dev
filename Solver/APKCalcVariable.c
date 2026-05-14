#include "APKCalcVariable.h"
/**
The function CalcVariable() is the control function which provides all the inputs necessary to the function UpdateVar() which
performs the Gauss Seidel iterations. Since this function is a generalised function, it takes variable tag as input
and loops over all the zones. If the variable tag corresponds to velocities it skips the zones which are identified as solids.
For all other variable tags it loops over all the zones and finds out the pointer for that variable tag in the domain
using the function MakeDomainArray(). Similarly it finds out the pointers in slabs corresponding to that particular variable tag.
These pointers are sent to the function UpdateVar(). This operation is carried out till the rms value of the difference between
the current iteration values and the previous iteration value over the whole domain is less than 1e-8.
Note: This function selectively loops over those zones which have more than 3 boundaries first. The intention in doing so is
to propagate the boundary information into the domain quickly.
*/
double CalcVariable(Base *base,int variableTag)
{
	Zone *zone;
	int totalcells,SolveSolidTag=1,MPFV_SolveTag=1;
	int iteru=0,numZones,nZones;
	double totalresidual,iterresidual;
	double convcri;
	nZones=(base->Number_of_zones);
	
	//NOTE: This will avoid problems even if one forgets inclusion of this statement before solving for the varialbe through CalcVariable() <<< Narsimha 19th April, 2008
	//CDcoeffCalculation(CgnsRoot.CgnsBase,variableTag);
	if(base->flagMPFV==1)
	{
		if(variableTag>=31 && variableTag<=38)
			MPFV_SolveTag=0;
		else if(variableTag>=41 && variableTag<=55)
			MPFV_SolveTag=0;
	}
		
	if(variableTag==1||variableTag==2||variableTag==3 ||MPFV_SolveTag==0||(base->flagTurb==1 && (variableTag==600||variableTag==601))||variableTag==62||variableTag==63||variableTag==64||variableTag==65||variableTag==710||variableTag==68||variableTag==69) //variableTag==67||
	{
		SolveSolidTag = 0;
		totalcells = base->TotalCellFluid;
		
	}
	else
	{
		SolveSolidTag = 1;
		totalcells = base->TotalCellDomain;
	}

	if(totalcells==0)
		totalcells=1;
	
	convcri =find_convcri(variableTag,base);
	//printf("Value returned = %e\n",convcri);
	if(base->BaseParameters->ConvectionScheme==1)//QUICK
	{
		if(base->BaseParameters->flag2d==1) //2D
		{
			do
			{
				iteru = iteru+1;
				totalresidual = 0.0;
				if(base->flagSolidification==1 && variableTag==7)
				{
					for(numZones=0;numZones<nZones;numZones++)
					update_dele_iteration(base,numZones);
					
					/**Position of update_dele_iteration is changed.  If this is called in zone
					 loop it's not possible to call CalsSource_opti along with this function.  As
					the CalcSource_opti is a 'base' level function.  To make this possible the position of update_dele_iteration is changed.  This did not cause any difference.**/
					
					CalcSource_opti(base,7);
				}

				for(numZones=0;numZones<nZones;numZones++)
				{
					zone = (base->CgnsZone+numZones);
					if(zone->flagSolidFluid==0 || SolveSolidTag == 1)
					{
						ExtractVariables(variableTag,base,zone);
						
						iterresidual = UpdateVarQVK2d(base,zone,variableTag);
						
						totalresidual = totalresidual+iterresidual;
						
						//printf("Zone: %s vartag=%d iterresidual=%e\n", zone->zonename,variableTag,iterresidual); 
					}
				}
				
				totalresidual = sqrt(totalresidual/totalcells);
				

				if(iteru%50 == 0)
					printf("GS IterCount=%d\tResidual of Variable[%d] = %e\n", iteru,variableTag,totalresidual);
			}
			while(totalresidual > convcri);
		}
		else if(base->BaseParameters->flag2d==0) //3D
		{
			do
			{
				iteru = iteru+1;
				totalresidual = 0.0;
				if(base->flagSolidification==1 && variableTag==7)
				{
					for(numZones=0;numZones<nZones;numZones++)
						update_dele_iteration(base,numZones);
					CalcSource_opti(base,7);
				}

				//totalresidual is used to find RMS residual for entire geometry rather than adding RMS residual for individual zones
				for(numZones=0;numZones<nZones;numZones++)
				{
					zone = (base->CgnsZone+numZones);
					if(zone->flagSolidFluid==0 || SolveSolidTag == 1)
					{
						ExtractVariables(variableTag,base,zone);
						
						iterresidual = UpdateVarQVK(base,zone,variableTag);
					
						totalresidual = totalresidual+iterresidual;
					}
				}
				totalresidual = sqrt(totalresidual/totalcells);

				if(iteru%50 == 0)
					printf("GS IterCount=%d\tResidual of Variable[%d] = %e\n", iteru,variableTag,totalresidual);
			}
			while(totalresidual > convcri);
		}
	}
	else if(base->BaseParameters->ConvectionScheme==0)//CD
	{
		if(base->BaseParameters->flag2d==1) //2D
		{
			do
			{
				iteru = iteru+1;
				totalresidual = 0.0;
				if(base->flagSolidification==1 && variableTag==7)
				{
					for(numZones=0;numZones<nZones;numZones++)
						update_dele_iteration(base,numZones);
					CalcSource_opti(base,7);
				}

				//totalresidual is used to find RMS residual for entire geometry rather than adding RMS residual for individual zones
				for(numZones=0;numZones<nZones;numZones++)
				{
					zone = (base->CgnsZone+numZones);
					if(zone->flagSolidFluid==0 || SolveSolidTag == 1)
					{
						ExtractVariables(variableTag,base,zone);
						
						iterresidual = UpdateVarCD2d(base,zone,variableTag);
						
						totalresidual = totalresidual+iterresidual;
					}
				}
				totalresidual = sqrt(totalresidual/totalcells);
				
				if(iteru%50 == 0)
					printf("GS IterCount=%d\tResidual of Variable[%d] = %e\n", iteru,variableTag,totalresidual);
			}
			while(totalresidual > convcri);
		}
		else if(base->BaseParameters->flag2d==0) //3D
		{
			do
			{
				iteru = iteru+1;
				totalresidual = 0.0;
				if(base->flagSolidification==1 && variableTag==7)
				{
					for(numZones=0;numZones<nZones;numZones++)
						update_dele_iteration(base,numZones);
					CalcSource_opti(base,7);
				}

				// totalresidual is used to find RMS residual for entire geometry rather than adding RMS residual for individual zones*/
				for(numZones=0;numZones<nZones;numZones++)
				{
					zone = (base->CgnsZone+numZones);
					if(zone->flagSolidFluid==0 || SolveSolidTag == 1)
					{
						ExtractVariables(variableTag,base,zone);
						
						iterresidual = UpdateVarCD(base,zone,variableTag);
					
						totalresidual = totalresidual+iterresidual;
					}
				}
				totalresidual = sqrt(totalresidual/totalcells);

				if(iteru%50 == 0)
					printf("GS IterCount=%d\tResidual of Variable[%d] = %e\n", iteru,variableTag,totalresidual);
			}
			while(totalresidual > convcri);
		}
	}
	
	//Is it necessary to call RefreshSurfaceCells() here even after calling ExtractVariables() above ??? <<< Narsimha 5th Feb, 2008 
	//Yes it is.  After calculation is over for all the zones, all the interfaces are to be refreshed once again such that most recent value will be available at the interfaces ...<<< Bhanesh 11th May, 2008.
	RefreshSurfaceCells(variableTag,base);
	
	return(totalresidual);
}
