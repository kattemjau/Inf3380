#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

void compute(double** matrix1, double** matrix2, double* matrix3, int arow, int bcols, int brow){
	double sum=0;
	int k, d, c;
	for (c = 0; c < arow; c++) {
		for (d = 0; d < bcols; d++){
			for (k = 0; k < brow; k++) {
				//1d array ikke 2d
				sum += (matrix1[c][k])*(matrix2[k][d]);
			}
			// printf("%d %d\n",c, d );
			matrix3[(c*arow)+d] = sum;
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

		int main(int argc, char *argv[]) {
			char* file1;
			char* file2;
			char* outputfile;
			double** matrix1;
			double** matrix2;
			double* matrix3;
			int num_rows1, num_cols1, num_rows2, num_cols2, num_rows3, num_cols3;
			int num_procs, my_rank;


			if(argc!=4){
				printf("Usage: ./'program' file1 file2 outputfile\n");
			}
			file1=argv[1];
			file2=argv[2];
			outputfile=argv[3];
			printf("\nFile1: %s File 2: %s File3: %s\n",file1, file2, outputfile);

			//start open mpi
			MPI_Init (&argc, &argv);
		  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
		  MPI_Comm_size (MPI_COMM_WORLD, &num_procs);

			if(my_rank==0){
			//read matrix 1 and 2
				read_matrix_binaryformat(file1, &matrix1, &num_rows1, &num_cols1);
				read_matrix_binaryformat(file2, &matrix2, &num_rows2, &num_cols2);
				//sets values
				if(num_rows1>=num_rows2){	num_rows3=num_rows1;
				}else{num_rows3=num_rows2;
				}if(num_cols1>=num_cols2){num_cols3=num_cols1;
				}else{num_cols3=num_cols2;}

				printf("Allocate\n" );
				printf("matrix 3: Rows: %d Cols: %d\n",num_rows3, num_cols3 );
				matrix3 = (double*)malloc((num_rows3)*(num_cols3)*sizeof(double));
			}

			MPI_Bcast (&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast (&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

			//dele opp arbeids oppgaver

			printf("Compute\n");
			// int arow, int bcols, int brow
			compute(matrix1, matrix2, matrix3, num_rows1, num_cols2, num_rows2);


			if(my_rank==0){
				printf("Write to file\n");
				//write matrix
				write_matrix_binaryformat(outputfile, matrix3, num_rows3, num_cols3);

				printf("deallocate\n" );
				deallocate(matrix1);
				deallocate(matrix2);
				free(matrix3);				
			}

			return 0;
		}
