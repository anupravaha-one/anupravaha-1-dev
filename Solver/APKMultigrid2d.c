#include "APKMultigrid2d.h"

double Multigrid2d(Base *base,int gdi)
{
	Zone *zone;
	int totalBaseCells=0,numZones,iZone;
//	int j,i,k,NI,NJ,NK,*ptrSize;
	int itr,pressureItercase,a,g;
	double residZone,residBase;
	
	numZones=base->Number_of_zones;
	totalBaseCells = base->TotalCellFluid;
	
	if(base->flagVariableDensityProblem == 0)
		pressureItercase = 0; /*case for calling UpdatePressure2d*/
		
	if(base->flagVariableDensityProblem == 1)
		pressureItercase = 2; /*case for calling UpdatePressure2dVarDensity*/
		
	if(gdi==(base->CoarseGridIndex))
	{
		do
		{
//		for(i=0;i<10;i++)
//		{
			residBase=0.0;
			for(iZone=0;iZone<numZones;iZone++)
			{
				zone=base->CgnsZone+iZone;
				if(zone->flagSolidFluid==0)
				{
					InterZoneVariableTransfer2d(4,base,zone,0);
					
					for(g=0;g<gdi;g++)
						DummyVariableRestrict2d(4,base,zone,g);

					MultigridExtractVariables(4,base,zone,gdi);
					
					switch(pressureItercase)
					{
						case 0:
								residZone=Smoothing2d(base,zone,gdi);
						break;
					}
					
					residBase=residBase+residZone;
					FillVariable2Dummy2d(4,base,zone,gdi);
					for(g=gdi;g>0;g--)
						DummyVariableProlongate2d(4,base,zone,g);
				}
			}
			residBase=(residBase/numZones);
//			printf("residual for corsest grid=%12.8e\n",residBase);
//		}
		}
		while(residBase > 1e-5);
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
							case 0:
									residZone=UpdatePressure2d(base,zone);
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
				
					ResidFind2d0(base,zone,gdi);
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
						InterZoneVariableTransfer2d(4,base,zone,0);
						
						for(g=0;g<gdi;g++)
							DummyVariableRestrict2d(4,base,zone,g);

						MultigridExtractVariables(4,base,zone,gdi);
						switch(pressureItercase)
						{
							case 0:
									residZone=Smoothing2d(base,zone,gdi);
							break;
						}
						residBase=residBase+residZone;
						FillVariable2Dummy2d(4,base,zone,gdi);
						for(g=gdi;g>0;g--)
							DummyVariableProlongate2d(4,base,zone,g);
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
				
					ResidFind2d(base,zone,gdi);
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
					Initialize_CoarseGrid2d((zone->grid+gdi)->V,(zone->grid+gdi)->GSize,0.0,gdi);
					Initialize_CoarseGrid2d((zone->grid+gdi)->DummyV,(zone->grid+gdi)->GSize,0.0,gdi);
					
					for(a=0;a<6;a++)
					{
						Initialize_CoarseGridSlabs2d(a,(zone->grid+gdi)->SlabV[a],(zone->grid+gdi)->GSize,0.0,gdi);
						Initialize_CoarseGridSlabs2d(a,(zone->grid+gdi)->DummySlabV[a],(zone->grid+gdi)->GSize,0.0,gdi);
					}
				}
			}
		}
//		for(i=0;i<=base->VWcycleflag;i++)
			Multigrid2d(base,gdi);
		
		for(iZone=0;iZone<numZones;iZone++)
		{
			zone=base->CgnsZone+iZone;
			if(zone->flagSolidFluid==0)
			{
				MultigridExtractVariables(4,base,zone,gdi);
				
				Prolongation2d(base,zone,gdi);
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
						InterZoneVariableTransfer2d(4,base,zone,0);
						
						for(g=0;g<gdi;g++)
							DummyVariableRestrict2d(4,base,zone,g);

						MultigridExtractVariables(4,base,zone,gdi);
						switch(pressureItercase)
						{
							case 0:
								residZone=Smoothing2d(base,zone,gdi);
							break;
						}
						residBase=residBase+residZone;
						FillVariable2Dummy2d(4,base,zone,gdi);
						for(g=gdi;g>0;g--)
							DummyVariableProlongate2d(4,base,zone,g);
					}
				}
				residBase=(residBase/numZones);
			}
		}
	}
	return 0;
}

double Initialize_CoarseGrid2d(double *P,int *Size,double constant,int gdi)
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

double Initialize_CoarseGridSlabs2d(int f,double *P,int *Size,double constant,int gdi)
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
