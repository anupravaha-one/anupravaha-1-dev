#include "APKPoint_to_Coord.h"

/**
This function is used to find coordinates that span the boundary using Point Range.
In other words the boundary face varies from this x,y,z coordinates to x',y',z' coordinates respectively.
*/ 
int point_coord(Zone *zone, int *point_range)
{
	int NI,NJ,NK,i,p[6];
	double c1,c2,c3,c4,c5,c6;
	double *x,*y,*z;
	char coordname1[33],coordname2[33],coordname3[33];
	NI=*(zone->Size);
	NJ=*(zone->Size+1);
	NK=*(zone->Size+2);
	for(i=0;i<6;i++)
	{
		p[i]=*(point_range+i);
	}
	for(i=0;i<33;i++)
	{
		coordname1[i]=zone->CgnsGrid->Array_Name1[i];
		coordname2[i]=zone->CgnsGrid->Array_Name2[i];
		coordname3[i]=zone->CgnsGrid->Array_Name3[i];
	}
	x=zone->CgnsGrid->X;
	y=zone->CgnsGrid->Y;
	z=zone->CgnsGrid->Z;
	c1= *(x+(NI*NJ)*(p[2]-1)+NI*(p[1]-1)+(p[0]-1));
	c2= *(y+(NI*NJ)*(p[2]-1)+NI*(p[1]-1)+(p[0]-1));
	c3= *(z+(NI*NJ)*(p[2]-1)+NI*(p[1]-1)+(p[0]-1));
	c4= *(x+(NI*NJ)*(p[5]-1)+NI*(p[4]-1)+(p[3]-1));
	c5= *(y+(NI*NJ)*(p[5]-1)+NI*(p[4]-1)+(p[3]-1));
	c6= *(z+(NI*NJ)*(p[5]-1)+NI*(p[4]-1)+(p[3]-1));
	printf("\nPoint Range is:%d %d %d-> %d %d %d",p[0],p[1],p[2],p[3],p[4],p[5]);
	printf("\nThe Span of the Boundary Patch is \n From %s=%lf \t %s=%lf \t %s=%lf  \n   To %s=%lf \t %s=%lf \t %s=%lf \n"
	       ,coordname1,c1,coordname2,c2,coordname3,c3,coordname1,c4,coordname2,c5,coordname3,c6);
	       return(0);
}

