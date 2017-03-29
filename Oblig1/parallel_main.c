#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

typedef struct{
  float** image_data; /* a 2D array of floats */
  int m; /* # pixels in y-direction */
  int n; /* # pixels in x-direction */
}image;

// make use of two functions from the simplejpeg library
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
    MPI_Init (&argc, &argv);
    3
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
    my_m = ...;
    my_n = ...;
    allocate_image (&u, my_m, my_n);
    allocate_image (&u_bar, my_m, my_n);
    /* each process asks process 0 for a partitioned region */
    /* of image_chars and copy the values into u */
    /* ... */
    convert_jpeg_to_image (my_image_chars, &u);
    iso_diffusion_denoising_parallel (&u, &u_bar, kappa, iters);
    /* each process sends its resulting content of u_bar to process 0 */
    /* process 0 receives from each process incoming values and */
    /* copy them into the designated region of struct whole_image */
    /* ... */
    if (my_rank==0) {
      convert_image_to_jpeg(&whole_image, image_chars);
      export_JPEG_file(output_jpeg_filename, image_chars, m, n, c, 75);
      deallocate_image (&whole_image);
}
