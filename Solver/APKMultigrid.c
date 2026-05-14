#include "APKMultigrid.h"

double Multigrid(Base *base,int gdi)
{
	Zone *zone;
	int totalBaseCells=0,numZones,iZone;
//	int j,i,k,NI,NJ,NK,*ptrSize;
	int itr,pressureItercase,a,g;
	double residZone,residBase;
	
	numZones=base->Number_of_zones;
	totalBaseCells = base->TotalCellFluid;
	
	if(base->flagVariableDensityProblem == 0)
		pressureItercase = 1; /*case for calling UpdatePressure, this is for normal 3d*/
		
	if(base->flagVariableDensityProblem == 1)
		pressureItercase = 3; /*case for calling UpdatePressureVarDensity, this is for 3d variable Density case*/
	
	if(gdi==(base->CoarseGridIndex))
	{
		do
		{
			residBase=0.0;
			for(iZone=0;iZone<numZones;iZone++)
			{
				zone=base->CgnsZone+iZone;
				if(zone->flagSolidFluid==0)
				{
					InterZoneVariableTransfer(4,base,zone,0);
					
					for(g=0;g<gdi;g++)
						DummyVariableRestrict(4,base,zone,g);

					MultigridExtractVariables(4,base,zone,gdi);
					
					switch(pressureItercase)
					{
						case 1:
							residZone=Smoothing(base,zone,gdi);
						break;
					}
					residBase=residBase+residZone;
					FillVariable2Dummy(4,base,zone,gdi);
					for(g=gdi;g>0;g--)
						DummyVariableProlongate(4,base,zone,g);
				}
			}
			residBase=(residBase/numZones);
			printf("residual for corsest grid=%12.8e\n",residZone);
		}
		while(residZone > 1e-5);
	}
	else
	{
		if(gdi == 0)
		{
			for(itr=0;itr<3;itr++)
			{
				residBase=0.0;
				for(iZone=0;iZone<numZones;iZone++)
				{
					zone=base->CgnsZone+iZone;
					if(zone->flagSolidFluid==0)
					{
						MultigridExtractVariables(4,base,zone,gdi);
						
						switch(pressureItercase)
						{
							case 1:
									residZone=UpdatePressure(base,zone);
							break;
						}
						residBase=residBase+residZone;
					}
				}
				residBase=(residBase/numZones);
			}
			for(iZone=0;iZone<numZones;iZone++)
			{
				zone=base->CgnsZone+iZone;
				if(zone->flagSolidFluid==0)
				{
					MultigridExtractVariables(4,base,zone,gdi);
				
					ResidFind0(base,zone,gdi);
				}
			}
		}
		else
		{
			for(itr=0;itr<3;itr++)
			{
				residBase=0.0;
				for(iZone=0;iZone<numZones;iZone++)
				{
					zone=base->CgnsZone+iZone;
					if(zone->flagSolidFluid==0)
					{
						InterZoneVariableTransfer(4,base,zone,0);
							
						for(g=0;g<gdi;g++)
							DummyVariableRestrict(4,base,zone,g);

						MultigridExtractVariables(4,base,zone,gdi);
						switch(pressureItercase)
						{
							case 1:
									residZone=Smoothing(base,zone,gdi);
							break;
						}
						residBase=residBase+residZone;
						FillVariable2Dummy(4,base,zone,gdi);
						for(g=gdi;g>0;g--)
							DummyVariableProlongate(4,base,zone,g);
					}
				}
				residBase=(residBase/numZones);
			}
			for(iZone=0;iZone<numZones;iZone++)
			{
				zone=base->CgnsZone+iZone;
				if(zone->flagSolidFluid==0)
				{
					MultigridExtractVariables(4,base,zone,gdi);
				
					ResidFind(base,zone,gdi);
				}
			}
		}
		
		for(iZone=0;iZone<numZones;iZone++)
		{
			zone=base->CgnsZone+iZone;
			if(zone->flagSolidFluid==0)
			{
				Restriction(base,zone,gdi);
			}
		}
		
		gdi++;
		
		if(gdi!=0)
		{
			for(iZone=0;iZone<numZones;iZone++)
			{
				zone=base->CgnsZone+iZone;
				if(zone->flagSolidFluid==0)
				{
					Initialize_CoarseGrid((zone->grid+gdi)->V,(zone->grid+gdi)->GSize,0.0,gdi);
					Initialize_CoarseGrid((zone->grid+gdi)->DummyV,(zone->grid+gdi)->GSize,0.0,gdi);
					
					for(a=0;a<6;a++)
						Initialize_CoarseGridSlabs(a,(zone->grid+gdi)->SlabV[a],(zone->grid+gdi)->GSize,0.0,gdi);
						Initialize_CoarseGridSlabs(a,(zone->grid+gdi)->DummySlabV[a],(zone->grid+gdi)->GSize,0.0,gdi);
				}
			}
		}
		
//		for(i=0;i<=base->VWcycleflag;i++)
			Multigrid(base,gdi);
		
		for(iZone=0;iZone<numZones;iZone++)
		{
			zone=base->CgnsZone+iZone;
			if(zone->flagSolidFluid==0)
			{
				MultigridExtractVariables(4,base,zone,gdi);
				
				Prolongation(base,zone,gdi);
			}
		}
		
		gdi--;
		
		for(iZone=0;iZone<numZones;iZone++)
		{
			zone=base->CgnsZone+iZone;
			if(zone->flagSolidFluid==0)
			{
				Correction(base,zone,gdi);
			}
		}
		
		if(gdi!=0)
		{
			for(itr=0;itr<3;itr++)
			{
				residBase=0.0;
				for(iZone=0;iZone<numZones;iZone++)
				{
					zone=base->CgnsZone+iZone;
					if(zone->flagSolidFluid==0)
					{
						InterZoneVariableTransfer(4,base,zone,0);
							
						for(g=0;g<gdi;g++)
							DummyVariableRestrict(4,base,zone,g);

						MultigridExtractVariables(4,base,zone,gdi);
						switch(pressureItercase)
						{
							case 1:
								residZone=Smoothing(base,zone,gdi);
							break;
						}
						residBase=residBase+residZone;
						FillVariable2Dummy(4,base,zone,gdi);
						for(g=gdi;g>0;g--)
							DummyVariableProlongate(4,base,zone,g);
					}
				}
				residBase=(residBase/numZones);
			}
		}
	}
	return 0;
}

double Initialize_CoarseGrid(double *P,int *Size,double constant,int gdi)
{
	int i,j,k,NI,NJ,NK;
	
	NI=*(Size+0);
	NJ=*(Size+1);
	NK=*(Size+2);
	
	for(k=0;k<NK-1;k++)
	for(j=0;j<NJ-1;j++)
	for(i=0;i<NI-1;i++)
	{
		*(P+(NI-1)*(NJ-1)*k+(NI-1)*j+i)=constant;
	}
	
	return(0);
}

double Initialize_CoarseGridSlabs(int f,double *P,int *Size,double constant,int gdi)
{
	int i,j,k,NI,NJ,NK;
	
	NI=*(Size+0);
	NJ=*(Size+1);
	NK=*(Size+2);
	
	switch(f)
	{
		case 0:
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			{
				*(P+(NJ-1)*k+j) = constant;
			}
		break;
		case 1:
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			{
				*(P+(NJ-1)*k+j) = constant;
			}
		break;
		case 2:
			for(k=0;k<NK-1;k++)
			for(i=0;i<NI-1;i++)
			{
				*(P+(NI-1)*k+i) = constant;
			}
		break;
		case 3:
			for(k=0;k<NK-1;k++)
			for(i=0;i<NI-1;i++)
			{
				*(P+(NI-1)*k+i) = constant;
			}
		break;
		case 4:
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				*(P+(NI-1)*j+i) = constant;
			}
		break;
		case 5:
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				*(P+(NI-1)*j+i) = constant;
			}
		break;
	}
	return(0);
}
