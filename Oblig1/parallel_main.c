#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

typedef struct{
  float** image_data; /* a 2D array of floats */
  int m; /* # pixels in y-direction */
  int n; /* # pixels in x-direction */
}image;

void allocate_image(image *u, int m, int n);
void deallocate_image(image *u);
void convert_jpeg_to_image(const unsigned char* image_chars, image *u);
void convert_image_to_jpeg(const image *u, unsigned char* image_chars);
void iso_diffusion_denoising(const image* u, const image *u_bar, float kappa, int iters, int start, int slutt);


void import_JPEG_file(const char *filename, unsigned char **image_chars,
  int *image_height, int *image_width,
  int *num_components);
void export_JPEG_file(const char *filename, unsigned char *image_chars,
  int image_height, int image_width,
  int num_components, int quality);


int main(int argc, char *argv[]){
  int m, n, c, iters;
  int my_m, my_n, my_rank, num_procs;
  float kappa;
  image u, u_bar, whole_image;
  unsigned char *image_chars, *my_image_chars;
  char *input_jpeg_filename, *output_jpeg_filename;

  //mpi initialization
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &num_procs);

  
    /* read from command line: kappa, iters, input_jpeg_filename, output_jpeg_filename */
    /* ... */
  if (my_rank==0) {
    import_JPEG_file(input_jpeg_filename, &image_chars, &m, &n, &c);
    allocate_image (&whole_image, m, n);
  }
  MPI_Bcast (&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    /* divide the m x n pixels evenly among the MPI processes */
  //size on max ruter
  my_m = m/num_procs;
  if(my_rank==(num_procs-1)){
    my_m+=m%num_procs;
  }
  my_n = n;

  allocate_image (&u, my_m, my_n);
  allocate_image (&u_bar, my_m, my_n);



  if(my_rank==0){
    //parrent sends
    int i, k=0;
    for(i=1; i<num_procs; i++)
    {
      k = m/num_procs*i;
      if(i==(num_procs-1)){
        k+=m%num_procs;
      }
      //drops id 0
      MPI_Send(&image_chars[k],k*n,MPI_UNSIGNED_CHAR, i, 1, MPI_COMM_WORLD);

      
    }

  }else{
    MPI_Recv(&my_image_chars[0], my_m*my_n ,MPI_UNSIGNED_CHAR, 0, 1, MPI_COMM_WORLD);
    convert_jpeg_to_image (my_image_chars, &u);
    iso_diffusion_denoising_parallel (&u, &u_bar, kappa, iters);
    convert_image_to_jpeg(&u, my_image_chars);
    //endres
    MPI_Send(&my_image_chars[0],my_m*my_n, MPI_UNSIGNED_CHAR, my_rank, 1, MPI_COMM_WORLD);
    
  }

  if(my_rank==0){
    int i;
    for(i=1; i<num_procs; i++){
      MPI_Recv(&image_chars[i*m/num_procs], i*n*m/num_procs ,MPI_UNSIGNED_CHAR, i, 1, MPI_COMM_WORLD);
      
    }
  }


  //sette sammen bildene

    /* each process sends its resulting content of u_bar to process 0 */
    /* process 0 receives from each process incoming values and */
    /* copy them into the designated region of struct whole_image */
    /* ... */
  if (my_rank==0) {
    convert_image_to_jpeg(&whole_image, image_chars);
    export_JPEG_file(output_jpeg_filename, image_chars, m, n, c, 75);
    deallocate_image (&whole_image);
  }
}




void allocate_image(image *u, int m, int n){
  //m er loddrett n er bortover
  u->n=n;
  u->m=m;

  u->image_data=(float**)malloc(m*sizeof(float*));
  int i;
  for(i=0; i<n; i++){
    u->image_data[i]=(float*)malloc(sizeof(float)*n);
  }

}
void deallocate_image(image *u){
  int i;
  for(i=0; i<u->m; i++){
    free(u->image_data[i]);
  }
  free(u->image_data);

}
void convert_jpeg_to_image(const unsigned char* image_chars, image *u){
  //convert 1d array to 2d array
  printf("%d, %d\n",u->m, u->n );
  int i, k;
  for (i = 0; i < u->n; i++) {
    for (k = 0; k < u->m; k++) {
      u->image_data[i][k]=image_chars[i*u->n+k];
    }
  }
}
void convert_image_to_jpeg(const image *u, unsigned char* image_chars){
  //make 2d array(struct) to 1d array
  int i, k;
  for (i = 0; i < u->m; i++) {
    for (k = 0; k < u->n; k++) {
      image_chars[i*u->n+k]=u->image_data[i][k];
    }
  }
}


void iso_diffusion_denoising(const image* u, const image *u_bar, float kappa, int iters, int start, int slutt){
  //smoothening function
  int i, p, j, k;

  //iterations
  for (p = 0; p < iters; p ++) {

    //start to stop m rows
    printf("start: %d, slutt: %d\n",start, slutt );
    for (i = start; i < slutt; i ++) {
      //n pos
      for (j = 0; j < u->n; j ++) {
        if(i==0){
          //top
          if(j==0){
            u_bar->image_data[i][j] = u->image_data[i][j] + kappa * (u->image_data[i + 1][j] +  u->image_data[i][j + 1] - 2 * u->image_data[i][j]);
          }else if(j==u->n-1){
            u_bar->image_data[i][j] = u->image_data[i][j] + kappa * (u->image_data[i + 1][j] +  u->image_data[i][j - 1] - 2 * u->image_data[i][j]);
          }else{
            u_bar->image_data[i][j] = u->image_data[i][j] + kappa * (u->image_data[i + 1][j]  + u->image_data[i][j - 1] + u->image_data[i][j + 1] - 3 * u->image_data[i][j]);
          }
        }else if(i==u->m-1){
          //bunnen
          if(j==0){
            u_bar->image_data[i][j] = u->image_data[i][j] + kappa * ( u->image_data[i - 1][j] + u->image_data[i][j + 1] - 2 * u->image_data[i][j]);
          }else if(j==u->n-1){
            u_bar->image_data[i][j] = u->image_data[i][j] + kappa * ( u->image_data[i - 1][j]  + u->image_data[i][j - 1]  - 2 * u->image_data[i][j]);
          }else{
            u_bar->image_data[i][j] = u->image_data[i][j] + kappa * ( u->image_data[i - 1][j]  + u->image_data[i][j - 1] + u->image_data[i][j + 1] - 3 * u->image_data[i][j]);
          }
        }else if(j==0){
          //venstre
          u_bar->image_data[i][j] = u->image_data[i][j] + kappa * (u->image_data[i - 1][j] + u->image_data[i + 1][j]  + u->image_data[i][j + 1] - 3 * u->image_data[i][j]);


        }else if(j==u->n-1){
          //hoyre
          u_bar->image_data[i][j] = u->image_data[i][j] + kappa * (u->image_data[i - 1][j] + u->image_data[i + 1][j]  + u->image_data[i][j - 1] - 3 * u->image_data[i][j]);

        }else{
          u_bar->image_data[i][j] = u->image_data[i][j] + kappa * (u->image_data[i - 1][j] + u->image_data[i + 1][j]  + u->image_data[i][j - 1] + u->image_data[i][j + 1] - 4 * u->image_data[i][j]);
        }



      }
    }

  //copy new to old to start agein
    for (i = 0; i < u->m; i++) {
      for (k = 0; k < u->n; k++) {
        u->image_data[i][k]=u_bar->image_data[i][k];
      }
    }
  }
}
