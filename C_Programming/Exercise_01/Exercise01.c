#include<stdio.h>
int main()
{
	int m=3,n=3,p=3,q=3,i,j,k; /* Assigning the values of the variables used */
	float first[3][3] = {{3,3,3}, /* Assigning the values of the two matrices used */
		{0,4,4},
		{0,0,5}}
	,second[3][3]={{0.33,-0.25,0},
		{0,0.25,-0.2},
		{0,0,0.2}},multiply[3][3],sum=0.00;
	
	if (n==p) /* Matrix multiplication condition checking */
	{
		for(i=0;i<m;i++){ /* Calculation of Matrix multiplication */
			for(j=0;j<q;j++){
				for(k=0;k<p;k++){
					sum=sum+first[i][k]*second[k][j];
				}
				multiply[i][j]=sum;
				sum=0;
				}}
				
		printf("Product of the matrices:\n"); /* Printing the values of the multiplication of two matrices */
 
	    for (i = 0; i < m; i++) {
	      for (j = 0; j < q; j++){
	        printf("%0.2f\t", multiply[i][j]);
	      }
	      printf("\n");			
		}}
	else /* Matrix multiplication condition fails */
	{
		printf("The size of the matrix is not correct \n");
	}	
	
}
