#include<stdio.h>
#include<stdlib.h>
#include<math.h>


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
	

	for(B=0;B<49;B++)
	{
		for(A=0;A<80;A++)
		{
			b=(i_start+(B/n)*(i_end-i_start));
			a=(r_start+(A/m)*(r_end-r_start));
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
	

	
}