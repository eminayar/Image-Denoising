/*
 Muhammed Emin Ayar
 Compile Status: Compiling
 Program Status: Working
 */
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>
#include <string>
#include <math.h>
#include <time.h>

int main(int argc, char** argv){
    
    MPI_Init(NULL,NULL);
    
    //initialize global variables
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int P=world_size-1;
    int i,j,r,t;
    int T=500000/(P);
    int M=200/sqrt(P);
    int N = 200;
    int KOK=sqrt(P);
    
    if( world_rank == 0 ){
        FILE* fp;
        char* line = NULL;
        char* token;
        size_t len = 0;
        ssize_t read;
        
        //get beta and pi from arguments
        double beta = atof(argv[3]);
        double pi = atof(argv[4]);
        
        //open input file
        fp = fopen(argv[1],"r");
        if( fp == NULL ){
            printf("No such file");
            MPI_Finalize();
            return 0;
        }
        
        //allocate memory for the main array
        int **arr = NULL;
        arr = (int **)malloc(sizeof(int*)*N);
        for(i = 0 ; i < N ; i++)
            arr[i] = (int *)malloc(sizeof(int) * N);
        
        //read the input into the array
        int row=0;
        while((read = getline(&line, &len, fp)) != -1 ){
            int column=0;
            while ((token = strsep(&line, " "))){
                arr[row][column]=atoi(token);
                column++;
            }
            row++;
        }
        
        fclose(fp);
        if( line ) free(line);
        
        //send beta and pi values to the slaves
        double beta_and_pi[2] ={beta,pi};
        for(i=1 ; i<=P ; i++ ){
            MPI_Send(beta_and_pi, 2, MPI_DOUBLE, i , 1 , MPI_COMM_WORLD);
        }
        
        //send array data to slaves
        int process=1;
        for(i=0;i<KOK;i++)
            for(j=0;j<KOK;j++){
                for( int k=0 ; k<M ; k++ )
                    MPI_Send( arr[i*M+k]+j*M , M , MPI_INT, process , k+2 , MPI_COMM_WORLD);
                process++;
            }
        
        //collect final state from slaves
        int* feedback = NULL;
        feedback = (int *)malloc(sizeof(int)*M*M);
        process=1;
        for(i=0;i<KOK;i++)
            for(j=0;j<KOK;j++){
                MPI_Recv(feedback, M*M, MPI_INT, process, T+M+10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                for(r=0;r<M*M;r++)
                    arr[i*M+r/M][j*M+r%M] = feedback[r];
                process++;
            }
        
        //write to the output file
        fp = fopen(argv[2],"w");
        for(i=0 ; i<N ; i++, fprintf(fp,"\n") )
            for(j=0; j<N ; j++)
                fprintf(fp,"%d ",arr[i][j]);
        
        fclose(fp);
        
    }else{
        
        //define variables
        double* beta_and_pi = NULL;
        int* subarr = NULL;
        int* initial = NULL;
        int* temp_arr = NULL;
        int* top_row = NULL;
        int* bottom_row = NULL;
        int* left_column = NULL;
        int* right_column = NULL;
        int* left_top = NULL;
        int* rigt_top = NULL;
        int* left_bottom = NULL;
        int* right_bottom = NULL;
        srand(time(0)*world_rank);
        
        //get beta and pi from master
        beta_and_pi = (double *)malloc(sizeof(double)*2);
        MPI_Recv(beta_and_pi, 2, MPI_DOUBLE, 0 , 1 , MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        double beta=beta_and_pi[0],pi=beta_and_pi[1];
        double gamma=0.5*log((1-pi)/pi);
        
        //allocate memory
        subarr = (int *)malloc(sizeof(int)*M*M);
        initial = (int *)malloc(sizeof(int)*M*M);
        top_row = (int *)malloc(sizeof(int)*M);
        temp_arr = (int *)malloc(sizeof(int)*M);
        bottom_row = (int *)malloc(sizeof(int)*M);
        left_column = (int *)malloc(sizeof(int)*M);
        right_column = (int *)malloc(sizeof(int)*M);
        left_top = (int *)malloc(sizeof(int));
        left_bottom = (int *)malloc(sizeof(int));
        right_bottom = (int *)malloc(sizeof(int));
        rigt_top = (int *)malloc(sizeof(int));
        
        //get initial data from master
        for(i=0 ; i<M ; i++)
            MPI_Recv(subarr+M*i, M , MPI_INT, 0 , i+2 , MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for(i=0 ; i<M*M ; i++)
            initial[i]=subarr[i];
        
        //define iteration variables
        int flag=0,rand_row=0,rand_column=0,val=0;
        double delta_E=0.0;
        float prob=0.0;
        
        for(t=0 ; t<T ; t++){
            flag=t+M+3; //message flag
            //first send messages to avoid deadlock
            if( world_rank > KOK ) //topmost
                MPI_Send(subarr, M , MPI_INT, world_rank-KOK, flag, MPI_COMM_WORLD);
            if( world_rank <= P-KOK ) //bottommost
                MPI_Send(subarr+M*(M-1), M, MPI_INT, world_rank+KOK , flag, MPI_COMM_WORLD);
            if( (world_rank-1)%KOK > 0 ){ //leftmost
                for(i=0;i<M;i++)
                    temp_arr[i]=subarr[i*M];
                MPI_Send(temp_arr, M , MPI_INT , world_rank-1, flag , MPI_COMM_WORLD);
            }
            if( (world_rank-1)%KOK < KOK-1 ){ //rightmost
                for(i=0;i<M;i++)
                    temp_arr[i]=subarr[i*M+M-1];
                MPI_Send(temp_arr, M , MPI_INT , world_rank+1, flag , MPI_COMM_WORLD);
            }
            if( world_rank > KOK && (world_rank-1)%KOK > 0 ) //upper_left
                MPI_Send(subarr, 1 , MPI_INT , world_rank-KOK-1 , flag , MPI_COMM_WORLD);
            if( world_rank > KOK && (world_rank-1)%KOK < KOK-1 ) //upper_rigt
                MPI_Send(subarr+M-1, 1 , MPI_INT , world_rank-KOK+1 , flag , MPI_COMM_WORLD);
            if( world_rank <= P-KOK && (world_rank-1)%KOK > 0 ) //bottom_left
                MPI_Send(subarr+(M-1)*M, 1, MPI_INT, world_rank+KOK-1, flag , MPI_COMM_WORLD);
            if( world_rank <= P-KOK && (world_rank-1)%KOK < KOK-1) //bottom_right
                MPI_Send(subarr+M*M-1, 1, MPI_INT, world_rank+KOK+1, flag, MPI_COMM_WORLD);
            
            //Now receive messages
            if( world_rank > KOK ) //topmost
                MPI_Recv(top_row, M , MPI_INT, world_rank-KOK, flag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if( world_rank <= P-KOK ) //bottommost
                MPI_Recv(bottom_row, M, MPI_INT, world_rank+KOK , flag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if( (world_rank-1)%KOK > 0 ) //leftmost
                MPI_Recv(left_column, M , MPI_INT , world_rank-1, flag , MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if( (world_rank-1)%KOK < KOK-1 ) //rightmost
                MPI_Recv(right_column, M , MPI_INT , world_rank+1, flag , MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if( world_rank > KOK && (world_rank-1)%KOK > 0 ) //upper_left
                MPI_Recv(left_top, 1 , MPI_INT , world_rank-KOK-1 , flag , MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if( world_rank > KOK && (world_rank-1)%KOK < KOK-1 ) //upper_rigt
                MPI_Recv(rigt_top, 1 , MPI_INT , world_rank-KOK+1 , flag , MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if( world_rank <= P-KOK && (world_rank-1)%KOK > 0 ) //bottom_left
                MPI_Recv(left_bottom, 1, MPI_INT, world_rank+KOK-1, flag , MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if( world_rank <= P-KOK && (world_rank-1)%KOK < KOK-1) //bottom_right
                MPI_Recv(right_bottom, 1, MPI_INT, world_rank+KOK+1, flag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            //select a point
            rand_row=rand()%M;
            rand_column=rand()%M;
            val = 0;
            //calculate adjacent values
            for(i=rand_row-1; i<=rand_row+1; i++)
                for(j=rand_column-1; j<=rand_column+1; j++ ){
                    if( i<0 && world_rank <= KOK ) continue;
                    if( i>=M && world_rank > P-KOK ) continue;
                    if( j<0 && (world_rank-1)%KOK <= 0 ) continue;
                    if( j>=M && (world_rank-1)%KOK >= KOK-1 ) continue;
                    if( i<0 && j<0 )
                        val+=left_top[0];
                    else if( i<0 && j>=0 && j<M )
                        val+=top_row[j];
                    else if( i<0 && j>= M )
                        val+=rigt_top[0];
                    else if( i>=0 && i<M && j<0 )
                        val+=left_column[i];
                    else if( i>=0 && i<M && j>=0 && j<M )
                        val+=subarr[i*M+j];
                    else if( i>=0 && i<M && j>=M )
                        val+=right_column[i];
                    else if( i>=M && j<0 )
                        val+=left_bottom[0];
                    else if( i>=M && j>=0 && j<M )
                        val+=bottom_row[j];
                    else if( i>=M && j>=M )
                        val+=right_bottom[0];
                    else
                        printf("DEBUG");
                }
            val -= subarr[rand_row*M+rand_column];
            delta_E= -2*gamma*initial[rand_row*M+rand_column]*subarr[rand_row*M+rand_column]-2*beta*subarr[rand_row*M+rand_column]*val;
            prob = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
            if( log(prob) < delta_E )
                subarr[rand_row*M+rand_column] = -subarr[rand_row*M+rand_column];
        }
        //send the final array to master
        MPI_Send(subarr, M*M , MPI_INT, 0 , T+M+10 , MPI_COMM_WORLD);
    }
    
    MPI_Finalize();
    
    return 0;
}

