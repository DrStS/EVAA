// C99
// Start program: mpirun -np 1 server
#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdbool.h>
#include <unistd.h> // needed for sleep() on POSIX system

#define MAX_DATA 100
int main( int argc, char **argv ) 
{ 
    int providedThreadSupport;
    bool terminateListening = false;
    char portName[MPI_MAX_PORT_NAME];
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &providedThreadSupport);
    if (MPI_THREAD_MULTIPLE != providedThreadSupport) {
        printf( "Requested MPI thread support is not guaranteed.\n");
    }
    MPI_Open_port(MPI_INFO_NULL, portName);
    printf("Server available at port:%s\n", portName); 
    #pragma omp parallel num_threads(2) shared(portName,terminateListening)
    {
        // Use OpemMP section construct for function parallelism
        #pragma omp sections
        {	
            #pragma omp section
            {
            // Do some work
            sleep(15);
            // Connect to yourself in order to terminate listening
            terminateListening = true;
            MPI_Comm dummy;
            MPI_Comm_connect(portName, MPI_INFO_NULL, 0, MPI_COMM_WORLD, &dummy);
            printf("Server is connected to itself.\n");
            MPI_Comm_disconnect(&dummy);
            printf("Server is disconnected.\n");
            MPI_Close_port(portName); 
            }            
            #pragma omp section
            {
            // Listening section
            while (1) {
                MPI_Comm interClient = MPI_COMM_NULL;
                MPI_Comm_accept(portName, MPI_INFO_NULL, 0, MPI_COMM_WORLD, &interClient);
                if (terminateListening == true) {
                    break;
                }
                MPI_Status status;
                char clientName[MAX_DATA];
                MPI_Recv(clientName, MAX_DATA, MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG, interClient, &status);
                printf("Client is connected with name: %s\n", clientName);
                MPI_Comm_disconnect(&interClient);
                printf("Client is disconnected.\n");
            }
            }
        } // End of sections
    } // End of parallel section
    MPI_Finalize();
    return (0);
}