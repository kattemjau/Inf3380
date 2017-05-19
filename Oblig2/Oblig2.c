#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
// #include <math.h>

void compute(double* matrix1, double* matrix2, double** matrix3, int arow, int bcols, int brow){
	double sum=0;
	int k, d, c;
	for (c = 0; c < arow; c++) {
		for (d = 0; d < bcols; d++){
			for (k = 0; k < brow; k++) {
				//1d array ikke 2d
				sum += matrix1[c*arow+k]*matrix2[d*bcols+k]; //TODO
			}
			// printf("%d %d\n",c, d );
			matrix3[c][d] = sum;
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
			printf("Rows: %d, COls: %d\n",rows, cols );
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

				if(argc!=4){
					printf("Usage: ./'program' file1 file2 outputfile\n");
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



				num_rows3=num_rows1;
				num_cols3=num_cols2;
				sqr=2;
				printf("matrix 3: Rows: %d Cols: %d\n",num_rows3, num_cols3 );
				// printf("sqr: %d\n",sqr );

			}

			MPI_Bcast (&sqr, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast (&num_cols1, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast (&num_rows3, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast (&num_cols3, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast (&num_rows2, 1, MPI_INT, 0, MPI_COMM_WORLD);

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
				printf("starter flipping\n");
				skrivOmB(matrixB2, matrixb, num_rows2, num_cols2);
				printf("ferdig\n");
				deallocate(matrixB2);

				int rows, cols;
				int startPos=0, start=0;

				int i, j;
		    for(i=0; i<sqr; i++){
					int litenRuteRow=num_rows3/sqr;
					int restRow=num_rows3%sqr;

					int litenRuteCol=num_cols3/sqr;
					int restCol = num_cols3%sqr;

					if(i<restRow*sqr){
						litenRuteRow++; //alle gjor dette
					}
					if(i<restCol*sqr){
						litenRuteCol++;
					}

					rows=i/sqr; //hvilken firkant den er i
					cols=i%sqr;

					for(j=0; j<sqr; j++){
						//send matrix 1		//array out of bound exception
						printf("Sending matrixA1, start: %d, size: %d, to: %d\n",rows*startPos,litenRuteRow*num_cols1 ,i*sqr+j+1 );
						MPI_Send(matrixA1[rows*startPos], litenRuteRow*num_cols1, MPI_DOUBLE, i*sqr+j+1, 1, MPI_COMM_WORLD);
						//send matrix 2
						// MPI_Send(&matrixb[cols*start*num_rows2], litenRuteCol*num_rows2, MPI_DOUBLE, j*sqr+i+1, 1, MPI_COMM_WORLD);

					}
					startPos+=litenRuteRow;
					start+=litenRuteCol;

				}
				deallocate(matrixA1);
				free(matrixb);

			}
			else{
				//allocate matrix
				int litenRuteRow=num_rows3/sqr;
				int restRow=num_rows3%sqr;

				int litenRuteCol=num_cols3/sqr;
				// printf("litenRuteCol %d\n", litenRuteCol );
				int restCol = num_cols3%sqr;

				if(my_rank<restRow*sqr){
					litenRuteRow++;
				}
				if(my_rank<restCol*sqr){
					litenRuteCol++;
				}
				// rows=my_rank/sqr;

				matrix1=(double*)malloc(litenRuteRow*num_cols1*sizeof(double));
				matrix2=(double*)malloc(litenRuteCol*num_rows2*sizeof(double));

				printf("My ranK: %d, recieved matrix1 with rows: %d, cols_: %d, sizeof: %d\n",my_rank, litenRuteRow, num_cols1, litenRuteRow*num_cols1);
				MPI_Recv(matrix1,litenRuteRow*num_cols1 , MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);

				// printf("My ranKK: %d, recieved matrix2 with cols: %d, cols_: %d, sizeof: %d\n",my_rank, litenRuteCol, num_rows2, litenRuteRow*num_rows2);
				// MPI_Recv(&matrix2[0],litenRuteCol*num_rows2, MPI_DOUBLE, 0, 1 , MPI_COMM_WORLD, &status);



				//allocate resultarray as 2d array or 1d array
				// matrix3=(double*)malloc((litenRuteRow)*(litenRuteCol)*sizeof(double));
				//
				// double** matrixC= (double**)malloc(den sin lille firkant * sizeof(double*));
				// int i
			  // for (i = 0; i < m; i++) {
			  //   if ((u->image_data[i] = (double*)malloc(sizeof(double)*n)) == NULL) {
			  //     perror("malloc(): N");
			  //     // free(u->image_data);
			  //     free(u);
			  //     exit(EXIT_FAILURE);
			  //   }
			  // }


				// printf("Compute\n");
				// int arow, int bcols, int brow TODO: needs to be changed
				// compute(matrix1, matrix2, matrixC, num_rows1, num_cols2, num_rows2);

				//send array back
				// MPI_Send(matrix3[0], pos*hos, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);

				free(matrix1);
				free(matrix2);

			}

			/*


			if(my_rank==0){
				// printf("Allocate\n" );
				printf("matrix 3: Rows: %d Cols: %d\n",num_rows3, num_cols3 );
				matrix3 = (double*)malloc((num_rows3)*(num_cols3)*sizeof(double));
				//recieve data from children
				int i;
				for(i=0; i<num_procs; i++){
					// motta data
					recv(matrix3[start], lengde, MPI_DOUBLE, 0, 1 , MPI_COMM_WORLD, &status);
					// sette sammendelene

				}




				printf("Write to file\n");
				//write matrix
				write_matrix_binaryformat(outputfile, matrix3, num_rows3, num_cols3);

				printf("deallocate\n" );
				free(matrix3);
			// }

*/
			MPI_Finalize();
			return 0;
		}
