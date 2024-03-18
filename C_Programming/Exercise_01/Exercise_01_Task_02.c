#include<stdio.h>
int main()
{
	int m,n,p,q,i,j,k; /* variables used */
	float first[10][10],second[10][10],multiply[10][10],sum=0.00; /* Matrix sizes assigned*/
	printf("Enter the number of rows and columns of first matrix \n");
	scanf("%d %d",&m,&n); /* getting the size of first matrix from user */
	printf("Enter the number of rows and columns of second matrix \n");
	scanf("%d %d",&p,&q); /* getting the size of second matrix from user */
	if (n==p) /* checking the matrix multiplication condition */
	{
		printf("Enter the elements of first matrix \n"); /* getting the elements of first matrix from user */
		for(i=0;i<m;i++){
			for(j=0;j<n;j++){
				scanf("%f",&first[i][j]);
			}}
			
		printf("Enter the elements of second matrix \n"); /* getting the elements of second matrix from user */
		for(i=0;i<p;i++){
			for(j=0;j<q;j++){
				scanf("%f",&second[i][j]);
			}}
  		
	  	printf("First matrix:\n"); /* printing the elements of first matrix */
 
	    for (i = 0; i < m; i++) {
	      for (j = 0; j < n; j++){
	        printf("%0.2f\t", first[i][j]);
	      }
	      printf("\n");			
		}
		
		printf("Second matrix:\n"); /* printing the elements of second matrix */
 
	    for (i = 0; i < p; i++) {
	      for (j = 0; j < q; j++){
	        printf("%0.2f\t", second[i][j]);
	      }
	      printf("\n");			
		}
		
		 /* Calculating matrix multiplication */
		for(i=0;i<m;i++){
			for(j=0;j<q;j++){
				for(k=0;k<p;k++){
					sum=sum+first[i][k]*second[k][j];
				}
				multiply[i][j]=sum;
				sum=0;
				}}
				
		printf("Product of the matrices:\n"); /* printing the multiplication of two matrices */
 
	    for (i = 0; i < m; i++) {
	      for (j = 0; j < q; j++){
	        printf("%0.2f\t", multiply[i][j]);
	      }
	      printf("\n");			
		}}
	else /* matrix multiplication condition fails */
	{
		printf("The size of the matrix is not correct \n");
	}	
	
}
