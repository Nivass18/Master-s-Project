#include<stdio.h>
#include<conio.h>

double madelbrotfunction(double *x, double *y, double *a, double *b){
	double t=*x;
	*x=(((*x)*(*x))-((*y)*(*y)))+(*a);
	*y=(2*t*(*y))+(*b);
	return *x,*y;
}

double magnitude(double *x, double *y){
	return sqrt(((*x)*(*x)+(*y)*(*y)));
}
int main()
{
	int r_start,r_end,i_start,i_end,i,A,B;
	double a,b,x,y,t,m=80,n=49;
	printf("MandelBrot set \n");
	printf("Enter the upper limit of real part ");
	scanf("%d",r_start);
	printf("Enter the lower limit of real part ");
	scanf("%d",r_end);
	printf("Enter the upper limit imaginary part ");
	scanf("%d",i_start);
	printf("Enter the upper limit imaginary part ");
	scanf("%d",i_end);
	for(B=0;B<=n;B++)
	{
		b=(i_start+(B/n)*(i_end-i_start));
		for(A=0;A<=m;A++)
		{
			a=(r_start+(A/m)*(r_end-r_start));
			x=0;
			y=0;
			for(i=0;i<=1000000;i++)
			{
				madelbrotfunction(&x,&y,&a,&b);
				if(magnitude(&x,&y)>2)
					break;	
			}
			if(i==1000001)
				printf("*");
			else
				printf("-");	
		}
		printf("\n");
	}
	
	
}