#include <stdio.h>
#include <stdlib.h>

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



void allocate_image(image *u, int m, int n){
  //m er horisontal n er vertikal
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


void iso_diffusion_denoising(const image* u, image *u_bar, float kappa, int iters){
  //smoothening function
  int i, p, j, k;

  for (p = 0; p < iters; p ++) {

  for (i = 1; i < (*u).m - 1; i ++) {
      for (j = 1; j < (*u).n - 1; j ++) {
          u_bar->image_data[i][j] = u->image_data[i][j] + kappa * (u->image_data[i - 1][j] + u->image_data[i + 1][j]  + u->image_data[i][j - 1] + u->image_data[i][j + 1] - 4 * u->image_data[i][j]);

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



int main(int argc, char *argv[]){
  printf("starting\n");
  int m, n, c, iters;
  float kappa;
  image u, u_bar;
  unsigned char *image_chars;
  char *input_jpeg_filename, *output_jpeg_filename;
  /* read from command line: kappa, iters, input_jpeg_filename, output_jpeg_filename */
  /* ... */
  if(argc<5){
    printf("Wrong args! Usage: ./filename <kappa> <iters> <inutfile> <outputfile>\n");
  }
  input_jpeg_filename=argv[3];
  output_jpeg_filename=argv[4];
  import_JPEG_file(input_jpeg_filename, &image_chars, &m, &n, &c);

  allocate_image (&u, m, n);

  allocate_image (&u_bar, m, n);

  convert_jpeg_to_image (image_chars, &u);

  iso_diffusion_denoising (&u, &u_bar, kappa, iters);

  convert_image_to_jpeg (&u_bar, image_chars);

  export_JPEG_file(output_jpeg_filename, image_chars, m, n, c, 75);

  deallocate_image (&u);
  deallocate_image (&u_bar);
  return 0;
}
