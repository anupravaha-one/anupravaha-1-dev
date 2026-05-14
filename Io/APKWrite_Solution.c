#include "APKWrite_Solution.h"

/**
This function Writes the restarted solution and starts from that point.
*/

int Write_Solution(Base *base,char *filename)
{
 
	int i,Coor_size,index_zone;
	int *ptrSize,NI,NJ,NK,solvable_variables,inner;
	Zone *zone;
	FILE *solution;
	double *variable;

	solvable_variables = base->total_var;
	
	solution = fopen(filename,"w");
	if (solution == NULL) {
	   printf("\nCould not open file SOLUTIONS/RESTART_DATA/solution.dat to write \n Probably missing SOLUTIONS/RESTART_DATA directory \n--- trying to create SOLUTIONS/RESTART_DATA directory.... \n");
	   system("mkdir SOLUTIONS/RESTART_DATA");
	   solution = fopen(filename,"w");
	   if (solution==NULL) {printf("\nError opening file SOLUTIONS/RESTART_DATA/solution.dat\n"); exit(1);
	   } else { printf("...Created 'SOLUTIONS/RESTART_DATA' directory\n");}
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
			fprintf(solution,"%.20e\n",*(variable+i));
			
			if((base->VARTAG[inner]==1) || (base->VARTAG[inner]==2) || (base->VARTAG[inner]==3))
			{
				variable = MakeDomainArray(zone,(base->VARTAG[inner])*100);
				for(i=0;i<Coor_size;i++)
					fprintf(solution,"%.20e\n",*(variable+i));
			}
		}
		
		
		if(base->flagTurb==1)
		{
			variable = zone->Flowsol->Tcoeff->EffViscoeff;
			for(i=0;i<Coor_size;i++)
				fprintf(solution,"%.20e\n",*(variable+i));
		}
		
		variable = zone->Flowsol->Massflux;
		Coor_size = Coor_size *6;
		for(i=0;i<Coor_size;i++)
			fprintf(solution,"%.20e\n",*(variable+i));
	
	}

	fclose(solution);
	return (0);
}



