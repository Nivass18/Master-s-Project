#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
 
 
void main(){
    int i,x=50;

    printf("Parallel loop with private(x)\n\n");
    #pragma omp parallel for private(x)
    for(i=0;i<=5;i++){
        x=i;
        printf("The Value of X inside the parallel loop with private variable is %d\n",x);
    }
    printf("The value of the X outside the parallel loop is %d\n", x);
    printf("\n");
    
    printf("Parallel loop with lastprivate(x)\n\n");
    #pragma omp parallel for lastprivate(x)
    for(i=0;i<5;i++){
        x=i;
        printf("The Value of X inside the parallel loop with last private variable is %d\n",x);
    }
    printf("The value of X outside the parallel loop is %d\n",x);
 
}