#include<stdio.h>
#include<math.h>
#include<omp.h>
#include<time.h>


#define size 4000

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
double Result_1[size][size];
double Result_2[size][1];
int i,j,k;
double norm_residual,residual_1;
double Result_3[size][1];
double Result_4[size][1];

void multiplication_of_matrices(int m, int n, int p,int q,double matrixA[m][n],double matrixB[p][q],double Result[m][q])
{
    double a=0;
    int k=0;
 
    #pragma omp parallel for ordered private(i,j,k) shared(matrixA,matrixB,Result) reduction (+:a)

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

    //Matrix_A_Creation with the help of first row of the matrix
    #pragma omp parallel for ordered private(i) shared(first_row) 
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
    #pragma omp parallel for ordered private(i,j) shared(A) 
    for (i = 0; i< size; i++){
        for (j = 0; j < size; j++){
            if (i == j)                        // for main diagonal elements
                A[i][j] = first_row[0];
            else if (i + j < size)             // for index less than size
                A[i][j] = first_row[i+j];
            else                                   // for index larger than size
                A[i][j] = first_row[2*(size-1) - (i + j)];
        }
    }

  

    //initialization o x_Value
    #pragma omp parallel for ordered private(i,j) shared(x) 
    for(i=0;i<size;i++){
        x[i][0]=0;
    }

    //initialization of B_Vector
    //#pragma omp parallel for ordered private(i,j) shared(b) 
    for(i=0;i<size;i++){
        b[i][0]=1;
    }

    //D_inv_Matrix creation from the Matrix A
    //#pragma omp parallel for ordered private(i,j) shared(Dinv) 
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


    //L_U Matrix using lower and upper triangular matrix of A
    #pragma omp parallel for ordered private(i,j) shared(L_U) 
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
    
    multiplication_of_matrices(size,size,size,size,Dinv,L_U,Result_1); //initial value of D^-1(L_U)
    
    int m=1;
    multiplication_of_matrices(size,size,size,m,Dinv,b,Result_2);//initial value of D^-1(b)
    norm_residual=1;

    #pragma omp single
    while(norm_residual>(1*pow(10,-3))) //converging condition
    {
        multiplication_of_matrices(size,size,size,m,Result_1,x,Result_3); //(D^-1(L_U)*X)

        for(i=0;i<size;i++){
            x_new[i][0]=-Result_3[i][0]+Result_2[i][0]; //calculating -D^-1(L_U)*X+D^-1*b
        }

        multiplication_of_matrices(size,size,size,m,A,x_new,Result_4);
        for(i=0;i<size;i++){
            residual[i][0]=b[i][0]-Result_4[i][0]; //calculation of the residual b-A*x
            
        }
        residual_1=0;
        for(i=0;i<size;i++){
            residual_1=residual_1 + residual[i][0]*residual[i][0]; //norm of the residual
        }
        norm_residual=sqrt(residual_1);
        printf("\nNew X value\t");
        for(i=0;i<size;i++){
            printf("%lf\t",x[i][0]); //printing the new value of the X
        }

        for(i=0;i<size;i++){
            x[i][0]=x_new[i][0];
        }
        
    
    }
    //Execution time calculation
    g=clock()-g;
    double f=((double)g)/CLOCKS_PER_SEC;
    printf("The function took %f seconds to execute \n",f);
    return 0;

}