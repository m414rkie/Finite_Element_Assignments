#include <stdio.h>
#include <string.h>
#include <mpi.h>

const int MAXstring = 100;

int main(void){

char Greeting[MAXstring];
int q;
int comm_sz;
int world_rank;
	
MPI_Init(NULL, NULL);

MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

if (world_rank != 0){
		
	printf("Hello World from Process %d of %d \n",world_rank,comm_sz);
	
	
	MPI_Send(Greeting, strlen(Greeting)+1,MPI_CHAR, 0, 0, MPI_COMM_WORLD);
	
} else {
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	printf("Greetings from process %d of %d\n", world_rank, comm_sz-1);

	for (q=1; q<comm_sz;q++){
		MPI_Recv(Greeting, MAXstring, MPI_CHAR, q, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		printf("%s", Greeting);
		printf("SubProcess\n");
	}
}
	
MPI_Finalize();

}