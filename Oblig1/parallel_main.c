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
  unsigned char *image_chars;
  char *input_jpeg_filename, *output_jpeg_filename;


  if(argc<5){
    printf("Wrong args! Usage: ./filename <kappa> <iters> <inutfile> <outputfile>\n");
  }
  //importing args
  printf("Importing args\n");
  kappa = atof(argv[1]);
  iters = atoi(argv[2]);
  input_jpeg_filename=argv[3];
  output_jpeg_filename=argv[4];

  //mpi initialization
  printf("Init MPI\n");
  MPI_Init (&argc, &argv); //oppretter childs
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank); //setter prosessid
  MPI_Comm_size (MPI_COMM_WORLD, &num_procs); // setter hvor mange prosessoer


  printf("My rank: %d\n",my_rank  );
  if (my_rank==0) {
    //importerer bilde
    printf("import\n" );
   import_JPEG_file(input_jpeg_filename, &image_chars, &m, &n, &c);
    //allokerer hele bildet
    printf("allocate\n");
    allocate_image (&whole_image, m, n);
    printf("allocate\n" );
    allocate_image (&u_bar, m, n);
    //image_chars er 1d arrayet
    printf("convert_jpeg_to_image\n"); //segfault
    convert_jpeg_to_image (image_chars, &whole_image);
  }
  printf("sending signal\n" );
  MPI_Bcast (&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
  printf("sending signal\n" );

  MPI_Bcast (&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
  //deler bildet
  printf("sending signal\n" );
  MPI_Bcast(whole_image.image_data, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

    /* divide the m x n pixels evenly among the MPI processes */

  //gi disse verdiene som startposisjoner for prosessene, & er value
 //dele inn i antall rader kis
  int start=my_rank*(m/num_procs);
  int slutt=start+(m/num_procs);
  //slutt verdi

  if(my_rank==(num_procs-1)){
    //siste delen av bilde
    slutt+=m%num_procs;
  }
  printf("m: %d, n: %d\n",m, n );

  printf("iso_diffusion_denoising\n");
  iso_diffusion_denoising(&whole_image, &u_bar, kappa, iters, start, slutt);



  // allocate_image (&u_bar, m, n);
    /* each process asks process 0 for a partitioned region */
    /* of image_chars and copy the values into u */
    /* ... */
    /* each process sends its resulting content of u_bar to process 0 */
    /* process 0 receives from each process incoming values and */
    /* copy them into the designated region of struct whole_image */
    /* ... */
  if (my_rank==0) {
    printf("something\n");
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
