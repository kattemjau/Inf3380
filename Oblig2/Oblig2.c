#include <stdio.h>
#include <stdlib.h>

void compute(double** matrix1, double** matrix2, double** matrix3, int arow, int bcols, int acols){
	int sum=0, k, d, c;
	for (c = 0; c < arow; c++) {
		for (d = 0; d < bcols; d++){
			for (k = 0; k < acols; k++) {
				//1d array ikke 2d
				sum += (*matrix1[c]+k)*(*matrix2[k]+d);
			}

			matrix3[c][d] = sum;
			printf("%x ",sum );
			sum = 0;
		}
		printf("\n");
	}
}
void allocateMatrix(double*** matrix, int num_rows, int num_cols){
	int i;
	*matrix = (double**)malloc((num_rows)*sizeof(double*));
	(*matrix)[0] = (double*)malloc((num_rows)*(num_cols)*sizeof(double));
	for (i=1; i<(num_rows); i++)
	(*matrix)[i] = (*matrix)[i-1]+(num_cols);
}
void deallocate(double** matrix){
	int i;
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
	void write_matrix_binaryformat (char* filename, double** matrix, int num_rows, int num_cols){
			FILE *fp = fopen (filename,"wb");
			fwrite (&num_rows, sizeof(int), 1, fp);
			fwrite (&num_cols, sizeof(int), 1, fp);
			fwrite (matrix[0], sizeof(double), num_rows*num_cols, fp);
			fclose (fp);
		}

		int main(int argc, char *argv[]) {
			char* file1;
			char* file2;
			char* outputfile;
			double** matrix1;
			double** matrix2;
			double** matrix3;
			int num_rows1;
			int num_cols1;
			int num_rows2;
			int num_cols2;
			int num_rows3;
			int num_cols3;


			if(argc!=4){
				printf("Usage: ./'program' file1 file2 outputfile\n");
			}
			file1=argv[1];
			file2=argv[2];
			outputfile=argv[3];
			printf("\nFile1: %s File 2: %s File3: %s\n",file1, file2, outputfile);

			//start open mpi

			//read matrix 1 and 2
			read_matrix_binaryformat(file1, &matrix1, &num_rows1, &num_cols1);
			read_matrix_binaryformat(file2, &matrix2, &num_rows2, &num_cols2);
			//sets values
			if(num_rows1>=num_rows2){	num_rows3=num_rows1;
			}else{num_rows3=num_rows2;
			}if(num_cols1>=num_rows2){num_cols3=num_cols1;
			}else{num_cols3=num_cols2;}

			printf("Allocate\n" );
			allocateMatrix(&matrix3, num_rows3, num_cols3);

			printf("Compute\n");
			// int arow, int bcols, int acols
			compute(matrix1, matrix2, matrix3, num_rows1, num_cols2, num_cols1);

			printf("Write to file\n");
			//write matrix
			write_matrix_binaryformat(outputfile, matrix3, num_rows3, num_cols3);

			printf("deallocate\n" );
			// deallocate(matrix1);
			// deallocate(matrix2);
			// deallocate(matrix3);

			return 0;
		}
