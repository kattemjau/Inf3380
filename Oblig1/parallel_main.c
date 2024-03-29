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
void iso_diffusion_denoising(const image* u, const image *u_bar, float kappa);


void import_JPEG_file(const char *filename, unsigned char **image_chars,
  int *image_height, int *image_width,
  int *num_components);
void export_JPEG_file(const char *filename, unsigned char *image_chars,
  int image_height, int image_width,
  int num_components, int quality);



int getRank(int m, int num_procs, int my_rank){
  int rank= m/(num_procs-1);
  if(my_rank==(num_procs-1)){
    rank+=m%(num_procs-1);
  }

  return rank;
}

int main(int argc, char *argv[]){
  int m, n, c, iters;
  int my_m, my_n, my_rank, num_procs;
  float kappa;
  image u, u_bar, whole_image;
  unsigned char *image_chars, *my_image_chars;
  char *input_jpeg_filename, *output_jpeg_filename;
  MPI_Status status;

  if(argc<5){
    printf("Wrong args! Usage: ./filename <kappa> <iters> <inutfile> <outputfile>\n");
  }
  kappa = atof(argv[1]);
  iters = atoi(argv[2]);
  input_jpeg_filename=argv[3];
  output_jpeg_filename=argv[4];


  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &num_procs);
  // printf("deler opp i: %d\n", num_procs-1);

  // printf("My rank is:%d\n",my_rank );
  if (my_rank==0) {
    import_JPEG_file(input_jpeg_filename, &image_chars, &m, &n, &c);
    allocate_image (&whole_image, m, n);
    // printf("hele bildet er %d x %d\n", m, n);
  }
  MPI_Bcast (&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
  //size on max ruter
  my_m =getRank(m, num_procs, my_rank);
  my_n=n;

  if(my_rank==0){
    //parrent sends
    int i, k=0, rows=0;
    for(i=1; i<num_procs; i++)
    {
      k = getRank(m, num_procs, i);
      // printf("I AM %d GOT %d FROM getRank\n", i, k);
      // printf("%d indeks: Sender size: %dx%d aka %d til childersn NO %d\n", rows*n, my_m, my_n, my_m*my_n, i);
      if(i==1 || i==num_procs-1){
        if(i==1){
          MPI_Send(&image_chars[rows*n],(k*n)+n,MPI_UNSIGNED_CHAR, i, 1, MPI_COMM_WORLD);  }
        else{
          MPI_Send(&image_chars[(rows-1)*n],k*n+n,MPI_UNSIGNED_CHAR, i, 1, MPI_COMM_WORLD);  }
      }else{
        MPI_Send(&image_chars[(rows-1)*n],(k*n)+n*2,MPI_UNSIGNED_CHAR, i, 1, MPI_COMM_WORLD);}
      rows+=k-1;
    }
  }
  else{
    if(my_rank==1 || my_rank==(num_procs-1)){
      my_m+=1;
    }else{
      my_m+=2;
    }
    allocate_image(&u, my_m, my_n);
    allocate_image(&u_bar, my_m, my_n);
    my_image_chars = (unsigned char*) malloc((my_m*my_n)*sizeof(unsigned char));
    MPI_Recv(&my_image_chars[0], my_m*my_n ,MPI_UNSIGNED_CHAR, 0, 1, MPI_COMM_WORLD, &status);

    convert_jpeg_to_image(my_image_chars, &u);
    int p;
    // printf("Jobber med denoising!\n" );
    for (p = 0; p < iters; p ++) {
      iso_diffusion_denoising(&u, &u_bar, kappa);
      // printf("%d\n",my_rank );
      if(my_rank!=1){
        MPI_Sendrecv(&u.image_data[1][0], my_n, MPI_FLOAT, my_rank-1, 1, &u.image_data[0][0], my_n ,MPI_FLOAT, my_rank-1, 1, MPI_COMM_WORLD, &status);
      }
      if(my_rank!=(num_procs-1)){
        MPI_Sendrecv(&u.image_data[my_m-2][0], my_n, MPI_FLOAT, my_rank+1, 1, &u.image_data[my_m-1][0], my_n ,MPI_FLOAT, my_rank+1, 1, MPI_COMM_WORLD, &status);
      }
    }
    convert_image_to_jpeg(&u, my_image_chars);
    if(my_rank==1 || my_rank==(num_procs-1)){
      my_m-=1;
    }else{
      my_m-=2;
    }
    printf("%d Sender : %d til ROOT\n", my_rank, my_m*my_n );
    if(my_rank==1){
      MPI_Send(&my_image_chars[0], my_m*my_n, MPI_UNSIGNED_CHAR, 0, 1, MPI_COMM_WORLD);
    }else{
      MPI_Send(&my_image_chars[my_n], my_m*my_n, MPI_UNSIGNED_CHAR, 0, 1, MPI_COMM_WORLD);
    }
      free(my_image_chars);
      deallocate_image(&u);
      deallocate_image(&u_bar);
  }
  if(my_rank==0){
    int i, k=0, rows=0;
    for(i=1; i<num_procs; i++){
        k = getRank(m, num_procs, i);
        MPI_Recv(&image_chars[n*rows],k*n,MPI_UNSIGNED_CHAR, i, 1, MPI_COMM_WORLD, &status);
        rows+=k;
        printf("Root mottok : %dx%d aka %d piksler fra: %d\n",k, n, k*n,i );
    }
  }
  if (my_rank==0) {
    printf("Skriver til %s\n", output_jpeg_filename );
    export_JPEG_file(output_jpeg_filename, image_chars, m, n, c, 75);
    deallocate_image (&whole_image);
    free(image_chars);
  }
  MPI_Finalize();
  return EXIT_SUCCESS;
}
void allocate_image(image *u, int m, int n){
  //m er nedover n er bortover
  u->n=n;
  u->m=m;

  if ((u->image_data = (float**)malloc(m * sizeof(float*))) == NULL) {
    perror("malloc(): M");
    // free(u);
    exit(EXIT_FAILURE);
  }
  int i;
  for (i = 0; i < m; i++) {
    if ((u->image_data[i] = (float*)malloc(sizeof(float)*n)) == NULL) {
      perror("malloc(): N");
      // free(u->image_data);
      free(u);
      exit(EXIT_FAILURE);
    }
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
  int i, k;
  for (i = 0; i < u->m; i++) {
    for (k = 0; k < u->n; k++) {
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
void iso_diffusion_denoising(const image* u, const image *u_bar, float kappa){
  //smoothening function
  int i, j, k;
    for (i = 0; i < u->m; i ++) {
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
          if(j==0){
            u_bar->image_data[i][j] = u->image_data[i][j] + kappa * ( u->image_data[i - 1][j] + u->image_data[i][j + 1] - 2 * u->image_data[i][j]);
          }else if(j==u->n-1){
            u_bar->image_data[i][j] = u->image_data[i][j] + kappa * ( u->image_data[i - 1][j]  + u->image_data[i][j - 1]  - 2 * u->image_data[i][j]);
          }else{
            u_bar->image_data[i][j] = u->image_data[i][j] + kappa * ( u->image_data[i - 1][j]  + u->image_data[i][j - 1] + u->image_data[i][j + 1] - 3 * u->image_data[i][j]);
          }
        }else if(j==0){
          u_bar->image_data[i][j] = u->image_data[i][j] + kappa * (u->image_data[i - 1][j] + u->image_data[i + 1][j]  + u->image_data[i][j + 1] - 3 * u->image_data[i][j]);
        }else if(j==u->n-1){
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
