#include<stdio.h>
#include<omp.h>

double f(double *x){
	return 1/(1+((*x)*(*x)));
}
int main(){
	int n,i;
	double a,b,h,x,sum=0,integral;
	printf("Enter the value of the sub-intervals ");
	scanf("%d",&n);
	printf("Enter the value of the upper limit ");
	scanf("%lf",&a);
	printf("Enter the value of the lower limit ");
	scanf("%lf",&b);
	h=(b-a)/(n);
	
	#pragma omp parallel for ordered private(i) shared(n,a,h)
	for(i=1;i<n;i++){
		x=a+i*h;
		sum=sum+f(&x);
	}
	 integral=(h/2)*(f(&a)+f(&b)+2*sum);
  /*Print the answer */
  #pragma omp ordered
  printf("\nThe integral is: %lf\n",integral);
}