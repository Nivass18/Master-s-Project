#include<stdio.h>
#include<conio.h>
#include<stdlib.h>

#define vertices 5 //fixing the number of vertices
#define infinity 9999 //fixing the value of infinity and to be used in cost matrix 
#define max 10 //size of the arrays

struct Edge_Vertex {
    int u,v,w; // variables to access the edges array
}; 


//Function declaration
void adj_matrix_init(int adj_matrix[][vertices]); //Initialization of adj matrix with zeros
void create_adj_matrix(int adj_matrix[][vertices],struct Edge_Vertex Edges[],int n); //Creating adj matrix with zeros(no connection between vertices) and ones (if there is a connection)
void print_matrix(int adj_matrix[][vertices]);//For printing the adj matrix
void dijkstra_algorithm(int cost_matrix[max][max],int distance[max],int predesor[max],int adj_matrix[][vertices], int start_node);//Algorithm to calculate distance from cost matrix
void print_cost_matrix(int cost_matrix[max][max]);//Printing cost matrix
void print_shortest_path(int distance[max],int predesor[max],int start_node);//Shortest path 

//Main function 
void main(){
    printf("\n");
    printf("The value of the maximum unknown distance is %d\n", infinity);
    int n,adj_matrix[vertices][vertices]; //initialization of adj matrix
    int cost_matrix[max][max],start_node,distance[max],predesor[max]; //initialization of cost matrix, distance and precedor 


    struct Edge_Vertex Edges[] = {
        {0,1,4},{0,2,2},{1,2,3},{2,1,1},{1,3,2},{1,4,3},{2,4,5},{2,3,4},{4,3,1}
    }; //Values from the example graph
    n=sizeof(Edges)/sizeof(Edges[0]); //calculating memory allocated for the struct Edge_Vertex data_type (2 int= 8, total 7 ele, 56/7-->8)

    
    adj_matrix_init(adj_matrix); // function call adj matrix with all zeros
    create_adj_matrix(adj_matrix,Edges,n); //function call adj matrix with zeros and ones
    print_matrix(adj_matrix); //printing the adj  matrix

    printf("\nEnter the starting node "); //getting the startnode from the user
    scanf("%d",&start_node);

    dijkstra_algorithm(cost_matrix,distance,predesor,adj_matrix,start_node); //function call to calculate the minimum distance

    printf("\nThe cost matrix for the given graph is \n\n");
    print_cost_matrix(cost_matrix);//priting cost matrix

    printf("\nThe shortest path using Dijkstra Algorithm \n\n");
    print_shortest_path(distance,predesor,start_node);//printing shortest path
}

void adj_matrix_init(int adj_matrix[][vertices]){
// Return Adj matrix with zeros
    int u,v;
    for(u=0;u<vertices;u++){
        for(v=0;v<vertices;v++){
            adj_matrix[u][v]=0;
         }
     }
}


void create_adj_matrix(int adj_matrix[][vertices],struct Edge_Vertex Edges[],int n){
    /* getting the vertex number and corresponding weight from the example graph and 
    assign it to the adj_matrix*/
    int i;
    for(i=0;i<n;i++){
       int u=Edges[i].u;
       int v=Edges[i].v;
       int w=Edges[i].w;
       adj_matrix[u][v]=w;

    }
     
}

void print_matrix(int adj_matrix[][vertices])
{ //Printing the adj matrix with zeors and the corresponding weight from the struct datatype
    printf("\nThe Adjacency matrix representation of the graph \n\n");
     int i,j;
     for(i=0;i<vertices;i++)
     {
        for(j=0;j<vertices;j++)
         {
        printf("%d \t",adj_matrix[i][j]);
         }
         printf("\n");
     }
}


void  dijkstra_algorithm(int cost_matrix[max][max],int distance[max],int predesor[max],int adj_matrix[][vertices],int source)
    /* Applying the shortest path algorithm and getting the
    shortest path from the source node to the each other nodes*/
{
    int visited_node[max],count,next_node,minimum_dist,i,j;
    //create the cost matrix and assign the weight (if there is a connection between edges) and infinity to other 
    for(i=0;i<vertices;i++){
        for(j=0;j<vertices;j++){
            if(adj_matrix[i][j]==0)
                cost_matrix[i][j]=infinity;
            else
                cost_matrix[i][j]=adj_matrix[i][j];
        }
    }

    for(i=0;i<vertices;i++){
        visited_node[i]=0; //assigning visited node to 0
        distance[i]=cost_matrix[source][i]; //assigning dist to the startnode row of cost matrix
        predesor[i]=source; //to get the path of the corresponding minimum dist

    }
    visited_node[source]=1; //visited status of the source index to 1
    distance[source]=0; //distance of the source is zero
    count=1;

    while(count<vertices-1){
        minimum_dist=infinity; //setting the minimum dist to infiinity

        //check all the values in the distance array with minium and staus of the index
        for(i=0;i<vertices;i++)
            if(distance[i]<minimum_dist&&!visited_node[i]){
                minimum_dist=distance[i];
                next_node=i;
            }
        
            visited_node[next_node]=1;
        //calculating the minium distance
            for(i=0;i<vertices;i++)
                if(!visited_node[i])
                    if(minimum_dist+cost_matrix[next_node][i]<distance[i]){
                        distance[i]=minimum_dist+cost_matrix[next_node][i];
                        predesor[i]=next_node;
                    }
        count++;

    }

}

void print_cost_matrix(int cost_matrix[max][max])
    /*Printing the cost matrix*/   
{
     int i,j;
     for(i=0;i<vertices;i++)
     {
        for(j=0;j<vertices;j++)
         {
            printf("%d \t",cost_matrix[i][j]);
         }
         printf("\n");
     }
}

void print_shortest_path(int distance[max],int predesor[max],int source){
    //Printing the shortest path from the source to other nodes
    int i,j;
    for(i=0;i<vertices;i++)
   
		if(i!=source)
		{
			printf("\nDistance of Node %d = %d ",i,distance[i]);
			printf("\nCorresponding path is = %d ",i);
			
			j=i;
			do
			{
				j=predesor[j];
				printf("<---%d",j);
			}while(j!=source);
	}
    printf("\n");
}



