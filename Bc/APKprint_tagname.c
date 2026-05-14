#include "APKprint_tagname.h"

/**
This function is used to convert BC tag number specified according to CGNS convention 
to corresponding BC tag name like 'inflow' etc, described in the header file cgnslib.h.
*/
int print_tagname(int n)
{
	switch(n)
	{
		case 0:
		printf("\n~~~~~~~~THIS BOUNDARY HAS BEEN TAGGED AS - NULL~~~~~~~~\n");
		break;
		case 1:
		printf("\n~~~~~~~~THIS BOUNDARY HAS BEEN TAGGED AS - USER DEFINED~~~~~~~~\n"); 
		break;
		case 2:
		printf("\n~~~~~~~~THIS BOUNDARY HAS BEEN TAGGED AS - AxisymmetricWedge~~~~~~~~\n");	
		break;
		case 3:
		printf("\n~~~~~~~~THIS BOUNDARY HAS BEEN TAGGED AS - DegenerateLine~~~~~~~~\n ");	
		break;
		case 4:
		printf("\n~~~~~~~~THIS BOUNDARY HAS BEEN TAGGED AS - DegeneratePoint~~~~~~~~\n");	
		break;
		case 5:
		printf("\n~~~~~~~~THIS BOUNDARY HAS BEEN TAGGED AS - Dirichlet~~~~~~~~\n");	
		break;
		case 6:
		printf("\n~~~~~~~~THIS BOUNDARY HAS BEEN TAGGED AS - Extrapolate~~~~~~~~\n");
		break;
		case 7:
		printf("\n~~~~~~~~THIS BOUNDARY HAS BEEN TAGGED AS - Farfield~~~~~~~~\n");
		break;
		case 8:
		printf("\n~~~~~~~~THIS BOUNDARY HAS BEEN TAGGED AS - General~~~~~~~~\n");	
		break;
		case 9:
		printf("\n~~~~~~~~THIS BOUNDARY HAS BEEN TAGGED AS - Inflow~~~~~~~~\n");
		break;
		case 10:
		printf("\n~~~~~~~~THIS BOUNDARY HAS BEEN TAGGED AS - InflowSubsonic~~~~~~~~\n");
		break;
		case 11:
		printf("\n~~~~~~~~THIS BOUNDARY HAS BEEN TAGGED AS - InflowSupersonic~~~~~~~~\n");
		break;
		case 12: 
		printf("\n~~~~~~~~THIS BOUNDARY HAS BEEN TAGGED AS  - Neumann~~~~~~~~\n");
		break;
		case 13:
		printf("\n~~~~~~~~THIS BOUNDARY HAS BEEN TAGGED AS - Outflow~~~~~~~~\n");
		break;
		case 14:
		printf("\n~~~~~~~~THIS BOUNDARY HAS BEEN TAGGED AS - OutflowSubsonic~~~~~~~~\n");
		break;
		case 15:
		printf("\n~~~~~~~~THIS BOUNDARY HAS BEEN TAGGED AS - OutflowSupersonic~~~~~~~~\n");
		break;
		case 16:
		printf("\n~~~~~~~~THIS BOUNDARY HAS BEEN TAGGED AS - SymmetryPlane~~~~~~~~\n");	
		break;
		case 17:
		printf("\n~~~~~~~~THIS BOUNDARY HAS BEEN TAGGED AS - SymmetryPolar~~~~~~~~\n");
		break;
		case 18:
		printf("\n~~~~~~~~THIS BOUNDARY HAS BEEN TAGGED AS - TunnelInflow~~~~~~~~\n");
		break;
		case 19:
		printf("\n~~~~~~~~THIS BOUNDARY HAS BEEN TAGGED AS - TunnelOutflow~~~~~~~~\n");
		break;
		case 20:
		printf("\n~~~~~~~~THIS BOUNDARY HAS BEEN TAGGED AS - Wall~~~~~~~~\n");
		break;
		case 21:
		printf("\n~~~~~~~~THIS BOUNDARY HAS BEEN TAGGED AS - WallInviscid~~~~~~~~\n");
		break;
		case 22:
		printf("\n~~~~~~~~THIS BOUNDARY HAS BEEN TAGGED AS - WallViscous~~~~~~~~\n");
		break;
		case 23:
		printf("\n~~~~~~~~THIS BOUNDARY HAS BEEN TAGGED AS - WallViscousHeatFlux~~~~~~~~\n");
		break;
		case 24:
		printf("\n~~~~~~~~THIS BOUNDARY HAS BEEN TAGGED AS - WallViscousIsothermal~~~~~~~~\n");
		break;
		case 25:
		printf("\n~~~~~~~~THIS BOUNDARY HAS BEEN TAGGED AS - FamilySpecified~~~~~~~~\n");
		break;
		}
		//printf("iam at end of tagname");
	return 0;
}	
	

