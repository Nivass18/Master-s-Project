#include<stdio.h>
#include<conio.h>

double b_values(i_start,i_end,B,n){
    double b[B];
    for(B=0;B<n;B++){
        b[B]=(i_start+(B/n)*(i_end-i_start));  
    }
    return *b;
}

int main()
{
    int r_start,r_end,i_start,i_end,i,A,B;
	double a,b,x,y,c,v,m=80,n=49;
	printf("MandelBrot set \n");
	printf("Enter the upper limit of real part ");
	scanf("%d",&r_start);
	printf("Enter the lower limit of real part ");
	scanf("%d",&r_end);
	printf("Enter the upper limit imaginary part ");
	scanf("%d",&i_start);
	printf("Enter the upper limit imaginary part ");
	scanf("%d",&i_end);
    printf("The values of b is %f \n");
}