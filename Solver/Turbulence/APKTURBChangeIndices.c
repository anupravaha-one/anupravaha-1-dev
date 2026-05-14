#include"APKTURBChangeIndices.h"

int ChangeIndices(Base *base,Zone *zone,int indices[6],int Do,int BlockCellCount)
{
	int NI,NJ,NK,direction;
	int countBc,numBc;
	ZoneBoundary *zBc;

	NI=*(zone->Size+0);
	NJ=*(zone->Size+1);
	NK=*(zone->Size+2);
	numBc=zone->Number_of_BC_Patches;

	indices[0]=0;
	indices[1]=NI-1;
	indices[2]=0;
	indices[3]=NJ-1;
	indices[4]=0;
	indices[5]=NK-1;

	if(Do == 1)
	{
		for(countBc=0;countBc<numBc;countBc++)
		{
			zBc = zone->ZoneBc+countBc;

			if(zBc->Boundary_Tag==22)
			{
				direction=zBc->Direction;
				switch(direction)
				{
					case 1://IMax face
						if(indices[1] >= BlockCellCount)
							indices[1] = indices[1]-BlockCellCount;
					break;
					case 2://IMin face
						indices[0] = indices[0]+BlockCellCount;
					break;
					case 3://JMax face
						if(indices[3] >= BlockCellCount)
							indices[3] = indices[3]-BlockCellCount;
					break;
					case 4://JMin face
						indices[2] = indices[2]+BlockCellCount;
					break;
					case 5://KMax face
						if(indices[5] >= BlockCellCount)
							indices[5] = indices[5]-BlockCellCount;
					break;
					case 6://KMin face
						indices[4] = indices[4]+BlockCellCount;
					break;

				}
			}
		}
	}

	return 17;
}



