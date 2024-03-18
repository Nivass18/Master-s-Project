#include<stdio.h>
#include <stdlib.h>
#include<math.h>
void matrix_multiplication(int m, int n, int p, int q,double a[m][n],double b[p][q],double c[m][q]);
int main()
{
    int n;
    printf("Enter the size of the matrixs, n value for n*n matrix\t");
    scanf("%d",&n,"\n");
    printf("The calculation of A-matrixs is being done here\n");
    double A[n][n];
    int a,b;
    double c;
    int D[n][n];
    double L_U[n][n];
    for(int i=0; i<n; i++)
    {
        for(int j=0; j<n; j++)
        {
            if(i==j)
            {
                A[i][j]=2,
                D[i][j]=A[i][j],
                L_U[i][j]=0;
            }
            else
            {   b=abs(i-j);
                a=-pow(2,b);
                c=pow(2,a);
                A[i][j]=-c,
                D[i][j]=0,
                L_U[i][j]=A[i][j];
            }
            
        }
    }

    printf("\n");
    double D_inv[n][n];
    for(int i=0; i<n;i++)
    {
        for(int j=0; j<n; j++)
        {
            if(i==j)
            {
                //printf("inverse");
                D_inv[i][j]= pow(D[i][j],-1);
            }
            else
            {
                D_inv[i][j]=0;
            }
            
        }
    }
    printf("L_U matrix is\n");
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
              printf("%f ",L_U[i][j]);
        }
        printf("\n");
    }
    double T_1[n][n];
    matrix_multiplication(n,n,n,n,D_inv,L_U,T_1);
    double b_matrix[n][1];
    for(int i=0; i<n;i++)
    {
        b_matrix[i][0]=1;
    }
    printf("\n");
    double T_2[n][1];
    int m=1;
    matrix_multiplication(n,n,n,m,D_inv,b_matrix,T_2);
    printf("\n");
    double u[n][1];
    double r[n][1];
    double norm_r;
    norm_r=1;
    double r_1;
    double u_new[n][1];
    for(int i=0; i<n; i++)
    {
        u[i][0]=0;
    }
    for(int i=0; i<n;i++)
    {
        printf("%lf\t",u[i][0]);
    }
    printf("\n");
    int count=0;
    double  T_3[n][1];
    double T_4[n][1];
    while(norm_r>(1*pow(10,-3)))
    {   count= count + 1;
        //printf("Entered the while loop");
        matrix_multiplication(n,n,n,m,T_1,u,T_3);
        //printf("T_3 is %d\n ",&T_3);
        for(int i=0;i<n;i++)
        {
            u_new[i][0]= - T_3[i][0] + T_2[i][0];
            //printf("u new is %d",u_new[i][0]); 
        }
        printf("\n");
        matrix_multiplication(n,n,n,m,A,u_new,T_4);
        for(int i=0; i<n; i++)
        {
            r[i][0]=b_matrix[i][0] - T_4[i][0];
            
        }
        r_1=0;
        for(int i=0; i<n; i++)
        {
            r_1= r_1 + r[i][0]*r[i][0];  
            //printf("r1 matrix is %lf",r[i][0]); 
        }
        norm_r=sqrt(r_1);
        printf("The r-value is:%lf\t",norm_r);
        printf("\n");
        printf("New u value\n");
        for(int i=0; i<n;i++)
        {
            printf("%lf\t",u[i][0]);
        }
        printf("\n");
        for(int i=0;i<n;i++)
        {
            u[i][0]=u_new[i][0];
        }
        //if(count>50)
        //{
        //    break;
        //}
    }
    return 0;
}
void matrix_multiplication(int m, int n, int p, int q,double a[m][n],double b[p][q],double c[m][q])
{
    double x=0;
    int k=0;
    //printf("Matrixs C\n");
    for(k=0;k<m;k++)
    {
        for(int i=0;i<q;i++)
        {  
            for(int j=0;j<p;j++)
            {   
                x=x+a[k][j]*b[j][i];
            }
            c[k][i]=x;
            x=0;
        }
    }
}
