#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include "mpi.h"

/*Note:
	No of processors = No of vertices – 2 
(No of processor should satisfy the equation e.g if No of processor = 5, then No of vertices to be entered is 7)*/


#define infinity 9999 //fixing the value of infinity and to be used in cost matrix 
#define max 50 //size of the arrays
/// The below function ia for the assignment of the cost_matrix(the weighted-undirectional graph)
void cost_matrix_fn(int N,int min , int cost_matrix[max][max])
{
    //getting the vertex number and corresponding weight from the random and assign it to the cost_matrix
    int i,j;
    int w;
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        { 
            if (i<j)
            {
                w = rand() % min;
                //for undirectional graph
                cost_matrix[i][j]=w;  
            }
            else 
            {
                cost_matrix[i][j]=0;
            }
        }
    }
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            cost_matrix[j][i]=cost_matrix[i][j];
        }
    } 

    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            if(cost_matrix[i][j]==0)
            {
                cost_matrix[i][j]=infinity;
            }
                
        }
    }
    //Priting the cost_matrix
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            printf("%d \t",cost_matrix[i][j]);
        }
        printf("\n");
    }
}

/// The function for calculation of the dijksthra algorithm
void  dijkstra_algorithm(int N,int cost_matrix[max][max],int shortest_distance[max],int start, int min, int visited_status[N],int predecor[max])
    /* Applying the shortest path algorithm and getting the
    shortest path from the source node to the each other nodes*/
{
    int position,n_time,d;
    // The nodes unvisited are given a status of 0 and the adjacent node distances from start point is initiallized to the shortest_distance.
   
     /// the distance matrixs first value should be less than this, a random assumption.
        for(int j=1;j<N;j++)
        {
            if(shortest_distance[j]<min && visited_status[j]==0)
            {
                min=shortest_distance[j]; // The minimum distances is found using a search logic with a<b and it should be unvisited.
                position=j;// The minimums index is stored to have the way to move forward to the next step and to update the status.
            }
        }
        visited_status[position]=1; // Sinces the node is visited the status is updated, so that they are not altered again.
        for(int j=1;j<N;j++) // The minimum value found is updated in the table/array shortest_distance
        {   d=min + cost_matrix[position][j]; // min=shortest_distance[position]
            if((d<shortest_distance[j]) && visited_status[j]==0)
                {
                    predecor[j] = position,
                    shortest_distance[j]=d;
                   
                }
        }
}
/// The below function is for printing the shortest path
/// The main function for the start of the dijkstra program 
int main(int argc, char** argv)
{
    srand(time(NULL));
    int size_Of_Cluster //size of the processors
        ,process_Rank, //rank of the processor
        min,N, //minimum value for random creation, and size of the array
        source; //start node
    const int root_rank=0;
    double st1,st2,st3,st4;
    int cost_matrix[max][max], //cost_matrix with random values
        distance[max],predecor[max], //distance array
        visited_status[max]; //status of the corresponding node
    //MPI initialization
    MPI_Init(&argc, &argv);
    //size and the rank of the processor
    MPI_Comm_size(MPI_COMM_WORLD, &size_Of_Cluster);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_Rank);
    //getting source, size of array and the minimum value for random creation and broadcasting to all processors
    if(process_Rank==root_rank)
    {    
        printf("Enter the Number of node\n");
        scanf("%d",&N);
    }
    MPI_Bcast(&N, 1 , MPI_INT, root_rank, MPI_COMM_WORLD);
    if(process_Rank==root_rank)
    {    
        printf("Enter the max value of the range of the weights for cost matrix\n");
        scanf("%d",&min);
    }
    MPI_Bcast(&min, 1 , MPI_INT, root_rank, MPI_COMM_WORLD);
    if(process_Rank==root_rank)
    {    
        printf("Enter the source node to start the calculation of the shortest path\n");
        scanf("%d",&source);
    }
    MPI_Bcast(&source, 1 , MPI_INT, root_rank, MPI_COMM_WORLD);
    //calculating cost matrix in processor rank zero and broadcasting to all processors
    st1=MPI_Wtime();
    if(process_Rank==0)
    {
        cost_matrix_fn(N,min ,cost_matrix);
    }
    MPI_Bcast(&cost_matrix,N*N,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    //initializing distance, visited status and predecor initial value in all processors
    st2=MPI_Wtime();
    for(int i=1;i<N;i++)
    {
        distance[i]=cost_matrix[source][i];     
        visited_status[i]=0;
        predecor[i]=source;
    } 
    MPI_Bcast(&distance,N,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&visited_status,N,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&predecor,N,MPI_INT,0,MPI_COMM_WORLD);
    //calculating minimum distance in all processors
    for(int i=0;i<N-2;i++)
    {
        if(process_Rank==i)
        {
            dijkstra_algorithm(N,cost_matrix,distance,source,min,visited_status,predecor);
            MPI_Bcast(&distance,N,MPI_INT,i,MPI_COMM_WORLD);
            MPI_Bcast(&visited_status,N,MPI_INT,i,MPI_COMM_WORLD);
        }
    }
    //printing the minimum distance and the corresponding path
    st3=MPI_Wtime();
    if(process_Rank==0)
    {
        for(int i=0;i<N;i++)
        {
            if(i!=source)
            {
                printf("\nDistance of Node %d = %d ",i,distance[i]);
                printf("\nCorresponding path is = %d ",i);
                
                int j=i;
                do
                {
                    j=predecor[j];
                    printf("<---%d",j);
                }while(j!=source);
            }
        }
    printf("\n");
    printf("End of the program \n");
    }
    st4=MPI_Wtime();
    if(process_Rank==0)
    {
        printf("Parallel Code with size of %d and took %f seconds to execute serial part \n",N,((st2-st1)+(st4-st3)));
    }
    MPI_Finalize();
    return 0;
}