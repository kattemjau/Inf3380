#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

void nummer1(){
  float sum=1;
  float n=2;
  for (;;) {
    sum=sum-(1/(pow(2,n)));
    n=n+2;
    sum=sum+(1/(pow(2,n)));
    n=n+2;
    if(sum>0.79 && sum<0.85){
      printf("forskjellen er: 4/5\n" );

    }else{
      printf("forskjellen er: %f\n", sum );

    }
  }


}

void nummer2(){
  int array[256]={0};
  int max=0;
  int min=0;
  for (size_t i = 0; i < 256; i++) {
    array[i]=rand();
    if(i==0){
      min=array[0];
    }
    if(array[i]>max){
      max=array[i];

    }
    if(min>array[i]){
      min=array[i];
    }

  }
  printf("min: %d\n", min );
  printf("max %d\n", max );

}

void darray(){
  int array[255][255];
  for (size_t i = 0; i < 255; i++) {

    memset(array[i], 0, sizeof(array[i]));

  }

  printf("array printes\n" );
}

void tredearray(){
  char ***array= (char***)malloc(12);
  for (size_t i = 0; i < 12; i++) {
    array[i]= (char**)malloc(12);
    for (size_t k = 0; k < 12; k++) {
      array[i][k]= (char*)malloc(12);

    }
  }
  array[1][1][1]='f';
  printf("schrar thingey: %c\n", array[1][1][1]);

  for (size_t i = 0; i < 12; i++) {
    for (size_t k = 0; k < 12; k++) {
      free(array[i][k]);
    }
    free(array[i]);
  }
  free(array);
}


int main() {
  // nummer1();
  // nummer2();
  // darray();
  tredearray();
  return 0;
}
