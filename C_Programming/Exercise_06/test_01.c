#include<stdio.h>
#include<math.h>
#include<time.h>


#define size 1000

double A[size][size];
double Dinv[size][size];
double x[size][1];
double L_U[size][size];//declare the relevant matrices
double Result[size][1];
double b[size][1];
double temp[size][1];
double first_row[size];
double residual[size][1];
double x_new[size][1];
int i,j,k;
double norm_residual,residual_1;

void multiplication_of_matrices(int m, int n, int p,int q,double matrixA[m][n],double matrixB[p][q],double Result[m][q])
{
    double a=0;
    int k=0;
 
    for(k=0;k<m;k++){
        for(i=0;i<q;i++){
            for(j=0;j<p;j++){
                a=a+matrixA[k][j]*matrixB[j][i];
            }
            Result[k][i]=a;
            a=0;
        }
            
    }
   
}

int main()
{

    clock_t g;
	g=clock();
    //Matrix_A_Creation
    for(i=0;i<size;i++){
        if(i==0){
            first_row[i]=1;
        }
        else{
            double a=pow(2,i);
            double denominator=pow(2,a);
            first_row[i]=-1/denominator;
        }
      
    }

    for (i = 0; i< size; i++){
        for (j = 0; j < size; j++){
            if (i == j)                        // special case for the main diagonal
                A[i][j] = first_row[0];
            else if (i + j < size)             // normal case for small indexes
                A[i][j] = first_row[i+j];
            else                                   // special case for large indexes
                A[i][j] = first_row[2*(size-1) - (i + j)];
        }
    }



    //x_Value
    for(i=0;i<size;i++){
        x[i][0]=0;
    }

    //B_Vector
    for(i=0;i<size;i++){
        b[i][0]=1;
    }

    //D_inv_Matrix
    for(i=0;i<size;i++){
        for(j=0;j<size;j++){
            if(i==j){
                Dinv[i][j]=1;
            }
            else{
                Dinv[i][j]=0;
            }
        }
    }


    //L_U Matrix
    for(i=0;i<size;i++){
        for(j=0;j<size;j++){
            if(i==j){
                L_U[i][j]=0;
            }
            else{
                L_U[i][j]=A[i][j];
            }
        }
    }
    double Result_1[size][size];
    multiplication_of_matrices(size,size,size,size,Dinv,L_U,Result_1);
    double Result_2[size][1];
    int m=1;
    multiplication_of_matrices(size,size,size,m,Dinv,b,Result_2);
    norm_residual=1;
    double Result_3[size][1];
    double Result_4[size][1];
    while(norm_residual>(1*pow(10,-3)))
    {
        multiplication_of_matrices(size,size,size,m,Result_1,x,Result_3);//multiply L+U and the approximation

        for(i=0;i<size;i++){
            x_new[i][0]=-Result_3[i][0]+Result_2[i][0];
        }

        multiplication_of_matrices(size,size,size,m,A,x_new,Result_4);
        for(i=0;i<size;i++){
            residual[i][0]=b[i][0]-Result_4[i][0];
            
        }
        residual_1=0;
        for(i=0;i<size;i++){
            residual_1=residual_1 + residual[i][0]*residual[i][0];
        }
        norm_residual=sqrt(residual_1);
        printf("\nNew X value\t");
        for(i=0;i<size;i++){
            printf("%lf\t",x[i][0]);
        }

        for(i=0;i<size;i++){
            x[i][0]=x_new[i][0];
        }
        
    
    }
    g=clock()-g;
    double f=((double)g)/CLOCKS_PER_SEC;
    printf("The function took %f seconds to execute \n",f);
    return 0;

}