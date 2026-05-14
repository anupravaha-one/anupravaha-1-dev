#include "APKRead_Solution.h"

/**
This function Reads the restarted solution and starts from that point.
*/

int Read_Solution(Base *base)
{
 
	int i,Coor_size,index_zone;
	int *ptrSize,NI,NJ,NK,solvable_variables,inner;
	Zone *zone;
	FILE *solution;
	double *variable;
	
	solvable_variables = base->total_var;
	solution = fopen("SOLUTIONS/RESTART_DATA/solution.dat","r");
	if(solution==NULL)
	{
		printf("File opening failed in Read_Solution.c\n\n");
		printf("Could not open file SOLUTIONS/RESTART_DATA/solution.dat\n");
		printf("Exiting ....\n");
		exit(12);
	}
	/* Loop Over Zone Nodes for reading the nodes under each Zone Node */ 
	for(index_zone=0;index_zone<base->Number_of_zones;index_zone++)
	{
		zone = base->CgnsZone+index_zone;
		ptrSize= zone->Size;
		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);
		Coor_size = (NI-1)*(NJ-1)*(NK-1);
		
		for(inner=0;inner<solvable_variables;inner++)
		{
			variable = MakeDomainArray(zone,base->VARTAG[inner]);
			for(i=0;i<Coor_size;i++)
			fscanf(solution,"%lf",(variable+i));
			
			if((base->VARTAG[inner]==1) || (base->VARTAG[inner]==2) || (base->VARTAG[inner]==3))
			{
				variable = MakeDomainArray(zone,(base->VARTAG[inner])*100);
				for(i=0;i<Coor_size;i++)
					fscanf(solution,"%lf",(variable+i));
			}
		}
		
		
		if(base->flagTurb==1)
		{
			variable = zone->Flowsol->Tcoeff->EffViscoeff;
			for(i=0;i<Coor_size;i++)
				fscanf(solution,"%lf\n",(variable+i));
		}
		
		variable = zone->Flowsol->Massflux;
		Coor_size = Coor_size *6;
		for(i=0;i<Coor_size;i++)
			fscanf(solution,"%lf",(variable+i));
	}

	fclose(solution);
	return (0);
}



