#include<stdio.h>
#include<conio.h>
#include<stdlib.h>
#include<time.h>

#define infinity 9999

void create_adj_matrix(int n,int adj_matrix[n][n]/*struct Edge_Vertex Edges[]*/){
    /* getting the vertex number and corresponding weight from the randomly generated value and 
    assign it to the adj_matrix*/
    int i,j;
    for(i=0;i<n;i++){
        for(j=0;j<n;j++)
        { 
            if (i<j){
                int w = rand() % 20; //creating random values
                //for undirectional graph
                adj_matrix[i][j]=w;
            }
            else {
                adj_matrix[i][j]=0;
            }
        }
    }
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            adj_matrix[j][i]=adj_matrix[i][j];
        }
    }   
}

void print_matrix(int n,int adj_matrix[][n])
    /*Printing the adj_matrix*/
{
    printf("\nThe Adjacency matrix representation of the graph \n\n");
     int i,j;
     for(i=0;i<n;i++)
     {
        for(j=0;j<n;j++)
         {
            printf("%d \t",adj_matrix[i][j]);
         }
         printf("\n");
     }
}

void  dijkstra_algorithm(int n,int cost_matrix[n][n],int distance[n],int predesor[n],int adj_matrix[][n],int start_node)
    /* Applying the shortest path algorithm and getting the
    shortest path from the source node to the each other nodes*/
{
    int visited_node[n],count,next_node,minimum_dist,i,j;
    ////create the cost matrix and assign the weight (if there is a connection between edges) and infinity to other
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            if(adj_matrix[i][j]==0)
                cost_matrix[i][j]=infinity;
            else
                cost_matrix[i][j]=adj_matrix[i][j];
        }
    }

    for(i=0;i<n;i++){
        visited_node[i]=0; //assigning visited node to 0
        distance[i]=cost_matrix[start_node][i]; //assigning dist to the startnode row of cost matrix
        predesor[i]=start_node; //to get the path of the corresponding minimum dist

    }
    visited_node[start_node]=1; //visited status of the source index to 1
    distance[start_node]=0; //distance of the source is zero
    count=1;

    while(count<n-1){
        minimum_dist=infinity;
        //setting the minimum dist to infiinity

        //check all the values in the distance array with minium and staus of the index
        for(i=0;i<n;i++)
            if(distance[i]<minimum_dist&&!visited_node[i]){
                minimum_dist=distance[i];
                next_node=i;
            }
        
            visited_node[next_node]=1;
            //calculating the minium distance
            for(i=0;i<n;i++)
                if(!visited_node[i])
                    if(minimum_dist+cost_matrix[next_node][i]<distance[i]){
                        distance[i]=minimum_dist+cost_matrix[next_node][i];
                        predesor[i]=next_node;
                    }
        count++;

    }

}

void print_cost_matrix(int n,int cost_matrix[n][n])
    /*Printing the adj_matrix*/   
{
     int i,j;
     for(i=0;i<n;i++)
     {
        for(j=0;j<n;j++)
         {
            printf("%d \t",cost_matrix[i][j]);
         }
         printf("\n");
     }
}

void print_shortest_path(int n,int distance[n],int predesor[n],int start_node){
    //Printing the shortest path from the source to other nodes
    int i,j;
    for(i=0;i<n;i++)
   
		if(i!=start_node)
		{
			printf("\nDistance of Node %d = %d ",i,distance[i]);
			printf("\nCorresponding path is = %d ",i);
			
			j=i;
			do
			{
				j=predesor[j];
				printf("<---%d",j);
			}while(j!=start_node);
	}
    printf("\n");
}


void main(){
    int n;
    printf("Enter the size of the array "); //getting the size of the array from the user
    scanf("%d",&n);
    printf("\n");
    printf("The value of the maximum unknown distance is %d\n", infinity);
    int adj_matrix[n][n];//initialization of adj matrix
    int cost_matrix[n][n],start_node,distance[n],predesor[n]; //initialization of cost matrix, distance and precedor 
  
    printf("\nEnter the starting node "); //getting the startnode from the user
    scanf("%d",&start_node);

    clock_t g;
	g=clock();
    create_adj_matrix(n,adj_matrix); //function call adj matrix with zeros and ones
    print_matrix(n,adj_matrix);

    dijkstra_algorithm(n,cost_matrix,distance,predesor,adj_matrix,start_node);  //function call to calculate the minimum distance

    printf("\nThe cost matrix for the given graph is \n\n");
    print_cost_matrix(n,cost_matrix); //priting cost matrix

    printf("\nThe shortest path using Dijkstra Algorithm \n\n");
    print_shortest_path(n,distance,predesor,start_node);//printing shortest path

    g=clock()-g;
    double f=((double)g)/CLOCKS_PER_SEC;
    printf("Serial Code with size of %d and took %f seconds to execute \n",n,f);
}
