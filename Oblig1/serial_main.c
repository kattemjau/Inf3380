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


int main(int argc, char *argv[]){
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

void allocate_image(image *u, int m, int n){
  //m er horisontal n er vertikal
  u->n=n;
  u->m=m;

  u->image_data=(float**)malloc(m*sizeof(float**));

  for(int i=0; i<n; i++){
      u->image_data[i]=malloc(sizeof(float*)*n);
  }

}
void deallocate_image(image *u){
  for(int i=0; i<n; i++){
      free(u->image_data[i]);
  }
  free(u->image_data);

}
void convert_jpeg_to_image(const unsigned char* image_chars, image *u){
  //make 1d array(struct) to 2d array

}
void convert_image_to_jpeg(const image *u, unsignes char* image_chars){
  //make 2d array(struct) to 1d array


}
void iso_diffusion_denoising(image *u, image *u_bar, float kappa, int iters){
  //smoothening function
}
import_JPEG_file(char* input_jpeg_filename, int image_chars, int m,int n,int c){
  //import picture and change input_jpeg_filename pointer
}

export_JPEG_file(char* output_jpeg_filename, int image_chars, int m, int n, int c, int tall){
  //export picture to output_jpeg_filename
}
