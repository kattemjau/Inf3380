#include <stdio.h>
#include <stdlib.h>

// void compute(double** matrix1, double** matrix2, double** matrix3, row, cols){
	// int sum=0, k, d, c;
	// for (c = 0; c < m; c++) {
	// 	for (d = 0; d < row; d++) {
	// 		for (k = 0; k < cols; k++) {
	// 			sum = sum + matrix1[c][k]*matrix2[k][d];
	// 		}
	//
	// 		matrix3[c][d] = sum;
	// 		sum = 0;
	// 	}
	// }
// }
void allocateMatrix(double*** matrix, int num_rows, int num_cols){
	int i;
	*matrix = (double**)malloc((num_rows)*sizeof(double*));
	(*matrix)[0] = (double*)malloc((num_rows)*(num_cols)*sizeof(double));
	for (i=1; i<(num_rows); i++)
	(*matrix)[i] = (*matrix)[i-1]+(num_cols);
}
void deallocate(double* matrix){
	int i;
	// free(*matrix);
	free(matrix);
}
void read_matrix_binaryformat (char* filename, double* matrix, int* num_rows, int* num_cols){
	int i;
	FILE* fp = fopen (filename,"rb");
	fread (num_rows, sizeof(int), 1, fp);
	fread (num_cols, sizeof(int), 1, fp);
	/* storage allocation of the matrix */
	// *matrix = (double**)malloc((*num_rows)*sizeof(double*));
	// (*matrix)[0] = (double*)malloc((*num_rows)*(*num_cols)*sizeof(double));
	// for (i=1; i<(*num_rows); i++)
	// (*matrix)[i] = (*matrix)[i-1]+(*num_cols);
	/* read in the entire matrix */

	matrix=(double*)malloc((*num_rows)*(*num_cols)*sizeof(double)); //allokerer plass til hele filen

	printf("allokert\n");
	fread(matrix, sizeof(double), (*num_rows)*(*num_cols), fp);
	printf("allokert\n");

	for(i=0; i<num_rows; i++){

		printf("%x ", matrix[i]);

		}


	fclose (fp);
}
void write_matrix_binaryformat (char* filename, double** matrix, int num_rows, int num_cols){
	FILE *fp = fopen (filename,"wb");
	fwrite (&num_rows, sizeof(int), 1, fp);
	fwrite (&num_cols, sizeof(int), 1, fp);
	fwrite (matrix[0], sizeof(double), num_rows*num_cols, fp);
	fclose (fp);
// }void printMat(double* matrix, int num_rows, int num_cols){
	int k;
	// for(i=0; i<num_rows; i++){
	printf("%d ", sizeof(matrix));
		// for(k=0;k<(num_cols*num_rows);k++){
		// }
		printf("\n");
	// }
}

int main(int argc, char *argv[]) {
	char* file1;
	char* file2;
	char* outputfile;
	double* matrix1;
	double* matrix2;
	// double** matrix3;
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
	// read_matrix_binaryformat(file1, matrix1, &num_rows1, &num_cols1);
	// read_matrix_binaryformat(file2, &matrix2, &num_rows2, &num_cols2);

	FILE* fp = fopen (file1,"rb");
	fread (&num_rows1, sizeof(int), 1, fp);
	fread (&num_cols1, sizeof(int), 1, fp);

	int total=(num_rows1*num_cols1);
	printf("Total %d\n", total );
	matrix1=malloc(total*sizeof(double));
	fread(matrix1, sizeof(double), total, fp);

	printf("matrix 1: Rows: %d Cols: %d\n",num_rows1, num_cols1 );

	int i;
	for (i = 0; i < total; i++) {
		printf("%d ", matrix1[i] );
	}
	printf("\n" );


	// printf("matrix 2: Rows: %d Cols: %d\n",num_rows2, num_cols2 );

	//sets values
	if(num_rows1>=num_rows2){	num_rows3=num_rows1;
	}else{num_rows3=num_rows2;
	}if(num_cols1>=num_cols2){num_cols3=num_cols1;
	}else{num_cols3=num_cols2;}

	// printf("matrix 3: Rows: %d Cols: %d\n",num_rows3, num_cols3 );

	// allocateMatrix(&matrix3, num_rows3, num_cols3);

	// printMat(matrix1, num_rows1, num_cols1);

	//solve matrix
	// solveMatrix(matrix3, matrix1, matrix2,num_rows3, num_cols3);
	//write matrix
	// write_matrix_binaryformat(outputfile, matrix3, num_rows3, num_cols3);
	// deallocate(matrix1);
	// deallocate(matrix2);
	// deallocate(matrix3);

	return 0;
}
