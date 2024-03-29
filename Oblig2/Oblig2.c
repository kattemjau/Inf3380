#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>
#include <omp.h>

OMP_NUM_THREADS=4;

void compute(double* matrix1, double* matrix2, double* matrix3, int acol, int bcols, int brow, int arow){
	double sum=0;
	int k, d, c, teller=0;
	// printf("arow: %d, bcols %d, brow %d, acols %d\n", arow, bcols, brow, acol);
	#pragma omp for
	for (c = 0; c < arow; c++) {
		for (d = 0; d < bcols; d++){
			for (k = 0; k < brow; k++) {
				//1d array ikke 2d
				sum += matrix1[c*acol+k]*matrix2[brow*d+k]; //TODO
			}
			// printf("%d\n",teller );
			matrix3[teller++]=sum;
			sum = 0;
		}
		// printf("\n");
	}
}
void deallocate(double** matrix){
	free(*matrix);
	free(matrix);
}
void read_matrix_binaryformat (char* filename, double*** matrix, int* num_rows, int* num_cols){
	int i;
	FILE* fp = fopen (filename,"rb");
	fread (num_rows, sizeof(int), 1, fp);
	fread (num_cols, sizeof(int), 1, fp);
	/* storage allocation of the matrix */
	*matrix = (double**)malloc((*num_rows)*sizeof(double*));
	(*matrix)[0] = (double*)malloc((*num_rows)*(*num_cols)*sizeof(double));
	for (i=1; i<(*num_rows); i++)

	(*matrix)[i] = (*matrix)[i-1]+(*num_cols);
	/* read in the entire matrix */
	fread ((*matrix)[0], sizeof(double), (*num_rows)*(*num_cols), fp);
	fclose (fp);
}
void write_matrix_binaryformat (char* filename, double* matrix, int num_rows, int num_cols){
	FILE *fp = fopen (filename,"wb");
	fwrite (&num_rows, sizeof(int), 1, fp);
	fwrite (&num_cols, sizeof(int), 1, fp);
	fwrite (matrix, sizeof(double), num_rows*num_cols, fp);
	fclose (fp);
}
void skrivOmB(double** matrix2, double* matrixB, int rows, int cols){
	int d, k, teller=0;
	// printf("Rows: %d, COls: %d\n",rows, cols );
	for (k = 0; k < cols; k++) {
		for (d = 0; d < rows; d++) {
			matrixB[teller++]=matrix2[d][k];
		}
	}
}



int main(int argc, char *argv[]) {
	char* file1;
	char* file2;
	char* outputfile;
	double* matrix1;
	double* matrix2;
	double* matrix3;
	int num_rows1, num_cols1, num_rows2, num_cols2, num_rows3, num_cols3;
	int num_procs, my_rank, sqr;
	MPI_Status status;
	double** matrixA1;
	double** matrixB2;

	//start open mpi
	MPI_Init (&argc, &argv);
	MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size (MPI_COMM_WORLD, &num_procs);

	// printf("My rank %d\n", my_rank );


	if(my_rank==0){

		if(argc!=5){
			printf("Usage: ./'program' file1 file2 outputfile squareroot\n");
		}
		file1=argv[1];
		file2=argv[2];
		outputfile=argv[3];
		printf("\nFile1: %s File 2: %s File3: %s\n",file1, file2, outputfile);

		//read matrix 1 and 2
		read_matrix_binaryformat(file1, &matrixA1, &num_rows1, &num_cols1);
		printf("read matrix1 with rows: %d, cols_: %d, sizeof: %d\n",num_rows1, num_cols1,num_cols1*num_rows1);

		read_matrix_binaryformat(file2, &matrixB2, &num_rows2, &num_cols2);
		printf("read matrix2 with rows: %d, cols_: %d, sizeof: %d\n",num_rows2, num_cols2,num_cols2*num_rows2);


		sqr=atoi(argv[4]);
		if((sqr*sqr+1)!=num_procs){
			printf("Sqr is wrong!\n");
			return EXIT_FAILURE;
		}
		// printf("sqr: %d\n",sqr );

	}

	MPI_Bcast (&sqr, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast (&num_cols1, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast (&num_rows1, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast (&num_cols2, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast (&num_rows2, 1, MPI_INT, 0, MPI_COMM_WORLD);

	num_rows3=num_rows1;
	num_cols3=num_cols2;
	// printf("matrix 3: Rows: %d Cols: %d\n",num_rows3, num_cols3 );

	// int pos=num_rows3*litenRuteRow;
	// int hos=num_cols3*litenRuteCol;
	//TODO: sjekke dette kobo har skrevet kode før
	// if (restRad%my_rank > restRad) {
	// 	litenRuteRow++;
	// }
	// int restCol = num_cols3%sqr;
	// //TODO: sjekke dette kobo har skrevet kode før
	// if (restCol%my_rank > restCol) {
	// 	litenRuteCol++;
	// }



	//dele opp arbeids oppgaver
	if(my_rank==0){
		double* matrixb=malloc(num_rows2*num_cols2*sizeof(double));
		// printf("starter flipping\n");
		skrivOmB(matrixB2, matrixb, num_rows2, num_cols2);
		// printf("ferdig\n");
		deallocate(matrixB2);

		//orginal strl

		int startPos=0, start=0;
		int i, j;
		// printf("Sender matrixA\n" );
		for(i=0; i<sqr; i++){
			int litenRuteRow=num_rows3/sqr;
			int restRow=num_rows3%sqr;

			if(i<restRow){
				litenRuteRow++; //alle gjor dette
			}
			for(j=0; j<sqr; j++){
				// rows=(i*sqr+j)/sqr;
				// printf("Sending matrixA1, start: %d, size: %d, to: %d\n",startPos,litenRuteRow*num_cols1 ,i*sqr+j+1 );
				MPI_Send(matrixA1[startPos], litenRuteRow*num_cols1, MPI_DOUBLE, i*sqr+j+1, 1, MPI_COMM_WORLD);
			}
			startPos+=litenRuteRow;
		}
		deallocate(matrixA1);

		// printf("Sender matrixB\n" );
		for(i=0; i<sqr; i++){
			int litenRuteCol=num_cols3/sqr;
			int restCol = num_cols3%sqr;

			if(i<restCol){
				litenRuteCol++;
			}
			for(j=0; j<sqr; j++){
				// cols=(j*sqr+i)%sqr;
				// printf("Sending matrixb, start: %d, size: %d, to: %d\n",start*num_rows2,litenRuteCol*num_rows2 ,(j*sqr+i+1));
				MPI_Send(&matrixb[start*num_rows2], litenRuteCol*num_rows2, MPI_DOUBLE, (j*sqr+i+1), 1, MPI_COMM_WORLD);
			}
			start+=litenRuteCol;
		}
		free(matrixb);

	}
	else{
		int litenRuteRow=num_rows3/sqr;
		int restRow=num_rows3%sqr;
		int litenRuteCol=num_cols3/sqr;
		int restCol = num_cols3%sqr;

		int col=(my_rank-1)%sqr;
		// printf("my rank: %d, col: %d, rest: %d\n",my_rank, col, restCol);

		if((my_rank-1)<restRow*sqr){
			litenRuteRow++; //alle gjor dette
		}
		if(col<restCol){
			litenRuteCol++;
		}


		//allocate matrix
		matrix1=(double*)malloc(litenRuteRow*num_cols1*sizeof(double));
		matrix2=(double*)malloc(litenRuteCol*num_rows2*sizeof(double));

		// printf("My ranK: %d, reading matrix with size: %d * %d\n",my_rank, litenRuteRow, litenRuteCol);
		MPI_Recv(&matrix1[0],litenRuteRow*num_cols1 , MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);

		MPI_Recv(&matrix2[0],litenRuteCol*num_rows2, MPI_DOUBLE, 0, 1 , MPI_COMM_WORLD, &status);
		// printf("My ranKK: %d, recieved matrix2 with cols: %d, cols_: %d, sizeof: %d\n",my_rank, litenRuteCol, num_rows2,litenRuteCol*litenRuteRow );



		//allocate resultarray as 2d array or 1d array
		matrix3=(double*)malloc((litenRuteRow)*(litenRuteCol)*sizeof(double));
		//
		// double** matrixC= (double**)malloc(litenRuteRow*sizeof(double*));
		// int i;
		// for (i = 0; i < litenRuteRow; i++) {
		//   if ((matrixC[i] = (double*)malloc(litenRuteCol*sizeof(double))) == NULL) {
		//     perror("malloc(): N");
		//     free(matrixC);
		// 		free(matrix1);
		// 		free(matrix2);
		//     exit(EXIT_FAILURE);
		//   }
		// }


		// printf("Compute\n");
		// int arow, int bcols, int brow TODO: needs to be changed
		compute(matrix1, matrix2, matrix3, num_cols1, litenRuteCol, num_rows2, litenRuteRow);

		//send array back
		// printf("my ranK: %d, SENDING DATA\n", my_rank);
		MPI_Send(&matrix3[0], litenRuteRow*litenRuteCol, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);

		// printf("Sending complete for: %d\n",my_rank );
		free(matrix1);
		free(matrix2);
		free(matrix3);

	}



	if(my_rank==0){
		// printf("Allocate\n" );
		printf("matrix 3: Rows: %d Cols: %d\n",num_rows3, num_cols3 );
		matrix3 = (double*)malloc((num_rows3)*(num_cols3)*sizeof(double));
		//recieve data from children
		int i, startPos=0, teller, bredde_til_neste;
		for(i=0; i<num_procs-1; i++){
			int litenRuteRow=num_rows3/sqr;
			int restRow=num_rows3%sqr;
			int litenRuteCol=num_cols3/sqr;
			int restCol = num_cols3%sqr;
			int col=i%sqr;
			// printf("my rank: %d, col: %d, rest: %d\n",i, col, restCol);
			if(i<restRow*sqr){
				litenRuteRow++; //alle gjor dette
			}
			if(col<restCol){
				litenRuteCol++;
			}

			// rows=i/sqr; //hvilken firkant den er i
			// cols=i%sqr;

			// int index=0;
			double tempArr[litenRuteRow*litenRuteCol];
			memset(tempArr, 0, sizeof(tempArr));
			// printf("PROC NUM %d liten rute col: %d, litenRuteRow, %d\n",i,litenRuteCol, litenRuteRow );
			// motta data TODO: dosent wait for data
			MPI_Recv(tempArr, litenRuteCol*litenRuteRow, MPI_DOUBLE, i+1, 1 , MPI_COMM_WORLD, &status);
			// sette sammendelene
			// for (j = 0; j < litenRuteRow; j++) {
			// 	for (k = 0; k < litenRuteCol; k++) {
			// 		printf("Legger inn i matrix3: %d, fra temp: %d\n",j*litenRuteCol+k+start, k*j+k);
			// 		matrix3[j*litenRuteCol+k+start]=tempArr[k*j+k];
			//
			// 	}
			// }
			// printf("sizeof: %d\n",sizeof(tempArr) );

			bredde_til_neste = num_cols3 - litenRuteCol;
			// printf(" bruhsni %d , ", bredde_til_neste );
			// printf("%d thinkg: %d\n", i, i/sqr+1);
      teller = startPos;
      //legge inn i res-array
      int j;
      for(j = 0; j < litenRuteCol*litenRuteRow; j++){
        matrix3[teller] = tempArr[j];
				// printf("%d ", j);
        if (((j+1) % litenRuteCol )== 0 && j!=0){
					// printf("%d teller: %d\n",i,  teller);
          teller += bredde_til_neste+1;
        } else{
          teller++;
        }
        // printf("TELLER: %d\n", teller);
      }

			if ((i+1) % sqr == 0) {
		        //  startPos += litenRuteRow*litenRuteCol; //- ((sqr-1)*litenRuteCol);
						startPos = teller - (num_cols3 - litenRuteCol);
		       } else {
		         startPos += litenRuteCol;
		       }

					//  printf("\n\n\n");
			// free(tempArr);
		}



		printf("Write to file\n");
		//write matrix
		write_matrix_binaryformat(outputfile, matrix3, num_rows3, num_cols3);

		printf("deallocate\n" );
		free(matrix3);
	}

	// printf("Killing thread: %d\n", my_rank);
	MPI_Finalize();
	return 0;
}
