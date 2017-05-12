#include <stdio.h>
#include <stdlib.h>

void computeNvidia(){
	// http://www.uio.no/studier/emner/matnat/ifi/INF3380/v17/undervisningsmateriale/inf3380-gpu-2015.pdf
/*	We want to compute Q=M×N, assuming
	Q, M, N are all square
	matrices of same size Width×Width
	Each matrix has a 1D contiguous data storage
	Naive kernel implementation
	__global__ void MatrixMulKernel(float* Md, float* Nd, float* Qd, int Width)
	{
		int Row = blockIdx.y
		*TILE_WIDTH + threadIdx.y;
		int Col = blockIdx.x
		*TILE_WIDTH + threadIdx.x;
		float Qvalue = 0;
		for (int k=0; k<Width; ++k)
		Qvalue += Md[Row
		*Width+k]
		* Nd[k
		*Width+Col];
		Qd[Row
		*Width+Col] = Qvalue;
	} */
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

			allocateMatrix(&matrix3, num_rows3, num_cols3);

			//solve matrix
			// solveMatrix(matrix3, matrix1, matrix2,num_rows3, num_cols3);
			//write matrix
			// write_matrix_binaryformat(outputfile, matrix3, num_rows3, num_cols3);
			deallocate(matrix1);
			deallocate(matrix2);
			deallocate(matrix3);

			return 0;
		}
