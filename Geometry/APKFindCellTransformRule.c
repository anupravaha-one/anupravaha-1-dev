#include "APKFindCellTransformRule.h"

/*    * ******************************************************************************************************* *
      + ------------------------------------------------------------------------------------------------------- +
      * Name of the function : FindCellTransformRule()                                                          *
      * Return Type          : Nil                                                                              *
      * Number of Parameters : 1                                                                                *
      * Datatype of Param.   : (Interface *)                                                                    *
      * Parameters           : 1) Pointer to Interface                                                          *
      * Calls the Functions  :-->FindMinIndices()                                                               *
      *                       -->GetDonorRange()                                                                * 
      *                       -->FindInterfaceCellVertices1                                                     *
      *                       -->FindInterfaceCellVertices2                                                     *  
      *                       -->FindAddConstant()                                                              *
      * Description          : The necessity of this function arises from the fact that the Transformation Rule *
      *                        is just a rotation matrix.We need to find Add Constant so as to tranform the     *
      *                        indices in language of current zone into language of adjoining zone.             *
      * Functions Description: FindInterfaceCellVertices1():                                                    * 
      *                                  Finds all the 8 Vertices of the starting Cell of an interface in the   *
      *                                  Adjoining  Zone in the notations of curent zone.                       * 
      *                                  Inputs: starting vertex of the interface in notations of the current zone*
      *                                  ( point_s ) and direction of the interface.                            *                            *
      *                        FindInterfaceCellVertices2():                                                    *
      *                                  Finds all the 8 Vertices of the starting Cell of an interface in the   *
      *                                  Adjoining  Zone in the notations of the Adjoining Zone.                *
      *                                  Inputs: All the 8 vertices in the notation of current zone(*v1 )       *
      *                                    and Trasformation Matrix (*trans).                                   *
      *                                 Output: All the 8 vertices in the notation of the Adjoining zone(*v2 )  *
      *                        FindMinIndices():                                                                 *
      *                                  The function is used to find cell index which is the minimum of all the*
      *                                  8 vertices which make a cell.                                          *
      *                                  Input:8 vertices of cell (* vertices)                                  *
      *                                  Output:Cell Index(*cell)                                               *
      *                        FindAddConstant():                                                               *
      *                                  Input: vertices of a cell in the notation of current zone(cell1)       *
      *                                  vertices of the same cell in the notation of the adjoining zone(*cell2)*
      *                                  tranformation matrix(*trans)                                           *
      *                                  Output:Additional constant(*addconst)                                  *
      * !!! IMPORTANT !!!--->  GetDonorRange():                                                                 *
      *                                  This function is necessiated because the Point Donor Range given by    *
      *                                  GridZ is not according to the standard definition given in             *
      *                                  CGNS Recommended Practice AIAA R-101-2002 Page.No.65.                  *
      *                                  This function corrects this anamoly.                                   *
      + ------------------------------------------------------------------------------------------- + */


/**
This function is needed because the Point Donor Range given by GridZ is not according to the standard 
definition given in CGNS Recommended Practice AIAA R-101-2002 Page.No.65. This function corrects this anamoly.
*/
int GetDonorRange(int transform[3],int pDR[6])
{
	int m,n,count,wMin[3],wMax[3];
	for(count=0;count<3;count++)
	{
		if(pDR[count]>pDR[count+3])
		{
			wMax[count]=pDR[count];
			wMin[count]=pDR[count+3];
		}
		else
		{
			wMax[count]=pDR[count+3];
			wMin[count]=pDR[count];
		}	
	}
	for(count=0;count<3;count++)
	{
		m=transform[count];
		n=fabs(m);
		if(m>0)
		{
			pDR[n-1]=wMin[n-1];
			pDR[n-1+3]=wMax[n-1];
		}
		else
		{
			pDR[n-1]=wMax[n-1];
			pDR[n-1+3]=wMin[n-1];
		}
	}
	return 0;
}

int GetPointRange(int transform[3],int pDR[6])
{
	int m,n,count,wMin[3],wMax[3];
	for(count=0;count<3;count++)
	{
		if(pDR[count]>pDR[count+3])
		{
			wMax[count]=pDR[count];
			wMin[count]=pDR[count+3];
		}
		else
		{
			wMax[count]=pDR[count+3];
			wMin[count]=pDR[count];
		}	
	}
	for(count=0;count<3;count++)
	{
		m=transform[count];
		n=fabs(m);
		
			pDR[n-1]=wMin[n-1];
			pDR[n-1+3]=wMax[n-1];
		
	
	}
	return 0;
}


/**
This function is used to find the minimum indices(i_min, j_min, k_min) out of a set of indices(i, j, k) of 
eight vertices of a cell.
*/
int FindMinIndices(int vertices[][3],int *cell)
{
	int pts[8][3];
	int i,j;
	
	for(i=0;i<8;i++)
	{
		for(j=0;j<3;j++)
		{
			pts[i][j]=vertices[i][j];
		}
	}
	//*finding minimum index*//
	for(i=0;i<3;i++)
		*(cell+i)=pts[1][i];
	
	for(j=0;j<8;j++)
	{
		for(i=0;i<3;i++)
		{
			if(*(cell+i) >= pts[j][i])
				*(cell+i) = pts[j][i];
		}
	}
	return 0;
}
/**
This function is used to produce the indices of the eight vertices of a cell that is the first cell of an interface 
and lies in the adjoining Zone. 

The index of the starting point of the Interface and the direction of the interface 
is given. The indices are found in the indexing convention followed in the present Zone. Depending on the direction 
of the interface, the minimum index is found. The other inidces are found by adding the indices by one appropriately.
*/
int FindInterfaceCellVertices1(int point_s[3],int Vertices[8][3],int dir)
{
	
	int point[3];
	int i;
	
	if((dir==1 || dir==3 )|| dir==5) /* 1-Imax, 3-Jmax, 5-Kmax */
	{
		for(i=0;i<3;i++)
		{
			point[i]=point_s[i];
		}
	}
	if(dir==2)/*2-Imin*/
	{
		point[0]= point_s[0]-1;
		point[1]= point_s[1];
		point[2]= point_s[2];
	}

	if(dir==4)/*4-Jmin*/
	{
		point[0]= point_s[0];
		point[1]= point_s[1]-1;
		point[2]= point_s[2];
	}
	if(dir==6)/*6-Kmin*/
	{
		point[0]= point_s[0];
		point[1]= point_s[1];
		point[2]= point_s[2]-1;
	}
	
	Vertices[0][0]=point[0];     	Vertices[0][1]=point[1];	Vertices[0][2]=point[2];
	Vertices[1][0]=point[0]+1;	Vertices[1][1]=point[1];	Vertices[1][2]=point[2];
	Vertices[2][0]=point[0];  	Vertices[2][1]=point[1]+1;	Vertices[2][2]=point[2]; 
	Vertices[3][0]=point[0]+1;    	Vertices[3][1]=point[1]+1;	Vertices[3][2]=point[2];  
	Vertices[4][0]=point[0];   	Vertices[4][1]=point[1];	Vertices[4][2]=point[2]+1;
	Vertices[5][0]=point[0]+1; 	Vertices[5][1]=point[1];	Vertices[5][2]=point[2]+1;
	Vertices[6][0]=point[0]; 	Vertices[6][1]=point[1]+1;	Vertices[6][2]=point[2]+1;
	Vertices[7][0]=point[0]+1;   	Vertices[7][1]=point[1]+1;	Vertices[7][2]=point[2]+1;
	
	return 0;
}      
/**
This function transforms the the indices of the cell VERTICES(not cell indices) found by FindInterfaceCellVertices1(). 
The transformation invloves the rotation of the index system by mutliplication with the transformation matrix and then 
translation of the index system. The transformation rule is given by I2=[T]*I1+(B2-[T]*B1) where I1 is the index 
vector in the notation of the present Zone, I2 is the index vector for the same point in the index notation of 
the adjoining Zone, [T] is the transforamtion matrix (rotation matrix), B1 is the index vector(in the notations 
of the present Zone) of the starting point of the interface and B2 is the index vector of the same point in the 
notation system of the adjoining Zone. 
*/

int FindInterfaceCellVertices2(int v1[8][3],int *trans,int v2[8][3],int *begin1,int *begin2)
{
	int k,j,i;
	for (k=0;k<8;k++)
	{
		for(j=0;j<3;j++)
		{ 
			v2[k][j]=0;
			for(i=0;i<3;i++)
				v2[k][j]=v2[k][j]+ *(trans+j*3+i)*(v1[k][i]-*(begin1+i));
			
			v2[k][j]=v2[k][j]+*(begin2+j);

		}
	}
	return 0;
}

/*int FindInterfaceCellVertices2(int v1[8][3],int trans[3][3],int v2[8][3],int *begin1,int *begin2)
{
	int k,j,i;
	for (k=0;k<8;k++)
	{
		for(j=0;j<3;j++)
		{ 
			v2[k][j]=0;
			for(i=0;i<3;i++)
				v2[k][j]=v2[k][j]+trans[j][i]*(v1[k][i]-*(begin1+i));
			
			v2[k][j]=v2[k][j]+*(begin2+j);

		}
	}
	return 0;
}*/

/**
This function calculates the Additive constant (translation vector for cell transformations).
 
As the input arguments it takes the indices of the starting cell of the interface which lies in the neighboring Zone in the
notation of the present Zone as well as in the notation of the neighboring Zone. Given the indices of the same 
cell in both the notation systems and the tarnsformation (rotation) matrix the translation vector can be calculated.
*/
int FindAddConstant(int cell1[3],int cell2[3],int *trans,int addconst[3])
{
	int temp[3];
	int i,j;
	for(i=0;i<3;i++)
	{
		temp[i]=0;
		for(j=0;j<3;j++)
			temp[i]=temp[i]+*(trans+i*3+j)*cell1[j];
		addconst[i]=cell2[i]-temp[i];
	}
	return 0;
}
/*int FindAddConstant(int cell1[3],int cell2[3],int trans[3][3],int addconst[3])
{
	int temp[3];
	int i,j;
	for(i=0;i<3;i++)
	{
		temp[i]=0;
		for(j=0;j<3;j++)
			temp[i]=temp[i]+trans[i][j]*cell1[j];
		addconst[i]=cell2[i]-temp[i];
	}
	return 0;
}*/
/**
The transformation rule for the cell indices is different from the vertex-indices transformation. The 
transformation matrix (rotation matrix) used remains the same except the translation vector changes, 
as the number of cells are not the same as the number of cells and moreover the index of a cell is the 
index of a vertex which has the minimum indices in the cell. This function is used to find the Additive 
constant (translation vector) for tranformation of the cell indices from the notation system of the present 
Zone to the notation system of the neighboring Zone.
 
The procedure used is as follows-

Knowing the starting point index (from the point range) of the interface in the notation of the present Zone, the 
indices of all the vertices of the starting cell (which lies in the neighboring Zone) of the interface are calculated
(in present notation system) using FindInterfaceCellVertices1().

Using the vertex indices transformation rule, the cell vertex indices found above are transformed into the notation 
system of the neighboring Zone using FindInterfaceCellVertices2().

The cell indices of the cell mentioned above are found by taking the minimum of the indices of the vertices using 
FindMinIndices(). Taking the minimum of the indices in both the notation systems gives us the indices of the cell 
in both the notation systems.

Finally the additional constant is found using FindAddConstant().
*/
  int FindCellTransformRule(Interface *face)
{
	int *pointStart1,*pointStart2,direction,vertices1[8][3],vertices2[8][3],*transformMatrix;
	int cell_1[3],cell_2[3];
	
	GetPointRange(face->Transform,face->Point_Range);

	pointStart1=face->Point_Range;

	GetDonorRange(face->Transform,face->Point_Donor_Range);
		
	pointStart2=face->Point_Donor_Range;
		
	direction=face->Direction;
		
	FindInterfaceCellVertices1(pointStart1,vertices1,direction);
		
	transformMatrix=face->Transformation_Matrix;

	FindInterfaceCellVertices2(vertices1,transformMatrix,vertices2,pointStart1,pointStart2);

	FindMinIndices(vertices1,cell_1);

	FindMinIndices(vertices2,cell_2);

	FindAddConstant(cell_1,cell_2,transformMatrix,face->addConstant);

	return 0;
}


