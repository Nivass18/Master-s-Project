#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>


double madelbrotfunction(double *x, double *y, double *a, double *b){
	double t=*x;
	*x=(((*x)*(*x))-((*y)*(*y)))+(*a);
	*y=(2*t*(*y))+(*b);
	return *x,*y;
}

double magnitude(double *x, double *y){
	return abs(((*x)*(*x)+(*y)*(*y)));
}
int main()
{
	int r_start,r_end,i_start,i_end,i,A,B;
	double a,b,x,y,t,m=80,n=49;
	printf("MandelBrot set \n");
	printf("Enter the upper limit of real part ");
	scanf("%d",&r_start);
	printf("Enter the lower limit of real part ");
	scanf("%d",&r_end);
	printf("Enter the upper limit imaginary part ");
	scanf("%d",&i_start);
	printf("Enter the upper limit imaginary part ");
	scanf("%d",&i_end);
    clock_t g;
	g=clock();
    //calculation of constant outside the for loop

    double temp_imaginary=i_end-i_start; 
    double temp_real=r_end-r_start;

	for(B=0;B<n;B++)
	{
		for(A=0;A<m;A++)
		{
			b=(i_start+(B/n)*(/*i_end-i_start*/ temp_imaginary));
			a=(r_start+(A/m)*(/*r_end-r_start*/ temp_real));
			x=0;
			y=0;
			for(i=0;i<=100000;i++)
			{
				madelbrotfunction(&x,&y,&a,&b);
				if(magnitude(&x,&y)>4)
					break;	
			}
			if(i==100001)
				printf("*");
			else
				printf("-");	
		}
		printf("\n");
	}
	g=clock()-g;
    double f=((double)g)/CLOCKS_PER_SEC;
    printf("Mandelbrot set took %f seconds to execute \n",f);

	
}