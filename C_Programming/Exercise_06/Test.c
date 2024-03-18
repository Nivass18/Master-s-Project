#include <stdlib.h>
#include <stdio.h>
 
// Define an array outside the scope of any function.
double x[20000][20000];
int elements = sizeof(x) / sizeof(*x);
 
int main(void)
{
	int i,j;
	double *p;
 
	// Do something with the array x.
 
	printf("Filling the x array...\n");
 
	for (i = 0; i < elements; i++)
	{
        for(j=0;j<elements;j++){
            	x[i][j] = 3.14;
        }
	
	}
 
	// Dynamically allocate an array from the heap.
 
	p = malloc(sizeof(double) * elements);
 
	if (p)			// test for success
	{
		// Do something with the array pointed to by p.
 
		printf("Filling the dynamically-allocated array...\n");
 
		for (i = 0; i < elements; i++)
		{
			p[i] = 3.14;
		}
 
		// Send the block back to the heap from whence it came.
 
		free(p);
	}
 
	return 0;
}