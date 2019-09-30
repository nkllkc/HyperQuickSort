#include "mpi.h"
#include <stdio.h>
#define MAX_SIZE 100000
#define MAX 99
#define PROCESS_ARRAY_SIZE 4
#define DIMENSION 4

int* generete_list(int myid, int n) {
  int* array = (int *) calloc(n, sizeof(int));
  int i;

  srand((unsigned) myid + 1);
  for (i = 0; i < n; i++) {
    array[i] = rand() % MAX;
  }

  return array;
}

void print_array(int* array, int n) {
  int i; 
  if (n == 0) {
    printf("\n");
    return;
  }

  for (i = 0; i < n - 1; i++) {
    printf("%d, ", array[i]);
  }

  printf("%d\n", array[n - 1]);
}

int compare(const void *a, const void *b){
   return  *(int*)a >= *(int*)b ? 1 : -1;
}

int partition(int *array, int n, const int pivot) {
  int i = 0, j = n - 1;
  while (true) {
    while ((array[j] > pivot) && (j >= 0)) j--;
    while ((array[i] <= pivot) && (i < n)) i++;
    if (i >= j) {
      return j;
    }
    int tmp = array[i];
    array[i] = array[j];
    array[j] = tmp;
  }
}

int main(int argc, char *argv[]) {
    MPI_Status status;
    int myid, partner;
    int l;
    int bitvalue, mask;
    int pivot;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    int* array = generete_list(myid, PROCESS_ARRAY_SIZE);
    printf("Starting array from process %d: ", myid);
    qsort(array, PROCESS_ARRAY_SIZE, sizeof(int), compare);
    print_array(array, PROCESS_ARRAY_SIZE);

    bitvalue = 1 << (DIMENSION - 1);
    mask = (1 << DIMENSION) - 1;
    for (l = DIMENSION; l >= 1; l--) {
      printf("myid %d, mask = %d\n", myid, mask);
      printf("myid & mask = %d\n", myid & mask);
      if ((myid & mask) == 0) {
        printf("TRALALA");
        pivot = array[PROCESS_ARRAY_SIZE/2];
      }
      MPI_Bcast(&pivot, 1, MPI_INT, 0, MPI_COMM_WORLD);

      j = partition(array, pivot);

      partner = myid ^ bitvalue;
      if ((myid & bitvalue) == 0) {
        MPI_Send(array + j + 1, n - j)
      } else {

      }

      mask = mask ^ bitvalue;
      bitvalue = bitvalue >> 1;
    }

    MPI_Finalize();
    return 0;
}
