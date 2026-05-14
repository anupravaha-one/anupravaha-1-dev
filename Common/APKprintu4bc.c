#include "APKprintu4bc.h"
/**
This function,printu4bc() prints the slab values of the required variable ( e.g., u ) into a data file, which goes into 'SOLUTIONS' 
folder. To this any other variable you wish to print into a data file can be added. 
*/

int printu4bc(Base *CgnsBase)
{
  FILE *solsu;    

        solsu=fopen("SOLUTIONS/printu4bc.dat","w");

	int NI,NJ,NK;
	int i,j,k;
	Zone *zone;
	int direction,index,zn;
        printf("\n Enter the zone number \n");
	 scanf("%d",&zn);
         printf("\n Enter the Bc index \n");
	 scanf("%d",&index);

 		zone =(CgnsBase->CgnsZone+zn); 
                NI=*(zone->Size);
                NJ=*(zone->Size+1);
		NK=*(zone->Size+2);
 	 		
 		{
 			direction= (zone->ZoneBc+index)->Direction;
 			if (direction==1 || direction==2)
			  {for(k=0;k<NK-1;k++)
						{
							for(j=0;j<NJ-1;j++)
							  {	fprintf(solsu,"%lf\n",*(((zone->ZoneBc+index)->ZSlab)->u+k*(NJ-1)+j));
							  }
						}
                          }
                          
 			else if (direction==3 || direction==4)
 		         {for(k=0;k<NK-1;k++)
						{
							for(i=0;i<NI-1;i++)
							  {	fprintf(solsu,"%lf\n",*(((zone->ZoneBc+index)->ZSlab)->u+k*(NI-1)+i));
							  }
						}
                          }
 			else if (direction==5 || direction==6)
 			 {for(j=0;j<NJ-1;j++)
						{
							for(i=0;i<NI-1;i++)
							  {	fprintf(solsu,"%lf\n",*(((zone->ZoneBc+index)->ZSlab)->u+j*(NI-1)+i));
							  }
						}
                          }
 		
 		
 		}
      

	fclose(solsu);
	return 0;
}


