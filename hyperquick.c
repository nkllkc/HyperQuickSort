#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>

#define DIMENSION 2
#define PROCESS_ARRAY_SIZE 4
#define MAX 99

int compare(const int *a, const int *b) {
    return *a - *b;
}

int* generate_array(int my_id) {
    int* array = (int*) calloc(PROCESS_ARRAY_SIZE, sizeof(int));
    int i;

    srand((unsigned) my_id + 1);
    for (i = 0; i < PROCESS_ARRAY_SIZE; i++) {
        array[i] = rand() % MAX;
    }

    return array;
}

void print_array(int* array, int size) {
    int i;

    if (size <= 0) {
        printf("Array is empty\n");
        return;
    }

    for (i = 0; i < size - 1; i++) {
        printf("%d, ", array[i]);
    }
    
    printf("%d\n", array[i]);
}

int partition(int* array, int lower_bound, int higher_bound, int pivot) {
    int i = lower_bound;
    int j = higher_bound;

    while(1) {
        while(array[i] <= pivot && i <= higher_bound) i++;
        while(array[j] > pivot && j >= lower_bound) j--;

        if (i >= j) {
            break;
        }

        int tmp = array[i];
        array[i] = array[j];
        array[j] = tmp;
    }

    return j;
}

void copy(int *source, int *destination, int size) {
    int i;
    for (i = 0; i < size; i++) {
        printf("copying element: %d\t", source[i]);
        destination[i] = source[i];
    }
    printf("\n");
}

void print_array_proc(int* array, int size, int id) {
    printf("%d: ", id);
    print_array(array, size);
}

int hyper_quick_sort(int* array, int my_id, int number_of_processes) {
    MPI_Status status;
    int mask, bitvalue, subcube_dimension;
    int pivot, number_of_elements;

    qsort(array, PROCESS_ARRAY_SIZE, sizeof(int), compare);
    print_array_proc(array, PROCESS_ARRAY_SIZE, my_id);

    number_of_elements = PROCESS_ARRAY_SIZE;
    
    bitvalue = 1 << (DIMENSION - 1);
    mask = (1 << DIMENSION) - 1;
    int initial_id = my_id;

    for (subcube_dimension = DIMENSION; subcube_dimension > 0; subcube_dimension--) {
        // Find the master of the subcube. Master has a leading 0.
        if ((my_id & mask) == 0) {
            pivot = rand() % 99;
        }

        int color = (1 << subcube_dimension) * initial_id / number_of_processes;
        MPI_Comm new_comm;
        MPI_Comm_split(MPI_COMM_WORLD, color, initial_id, &new_comm);
        int new_id, num_of_proc_in_group;
        MPI_Comm_size(new_comm, &num_of_proc_in_group);
        MPI_Comm_rank(new_comm, &new_id);

        MPI_Bcast(&pivot, 1, MPI_INT, 0, new_comm);
        printf("My_id: %d. Pivot: %d.\n", my_id, pivot);

        int j = partition(array, 0, number_of_elements - 1, pivot);

        print_array_proc(array, number_of_elements, my_id);
        int left_subarray_size = 
            j < 0 ? 0 : j + 1;
        int right_subarray_size = 
            j < 0 ? number_of_elements : number_of_elements - left_subarray_size; 

        printf("J: %d, l: %d, r: %d\n", j, left_subarray_size, right_subarray_size);

        // Just flip the highest bit of the id.
        int partner = my_id ^ bitvalue;
        
        // If this is left partner. 
        if ((new_id && bitvalue) == 0) {
            printf("Left\n");
            
            // Send the size of right sublist to the partner, 
            // so it can allocate the buffer.
            MPI_Send(&right_subarray_size, 1, MPI_INT, partner, 10, new_comm);
            
            // Receive the size of the partner's left sublist.
            int partners_sublist_size;
            MPI_Recv(&partners_sublist_size, 1, MPI_INT, partner, 10, new_comm, &status);
            printf("Partners sublist size: %d\n", partners_sublist_size);
            
            // Prepare the buffer for receiving new data.
            number_of_elements = left_subarray_size + partners_sublist_size;
            
            int *old_array = array;
            array = (int*) calloc(number_of_elements, sizeof(int));
            copy(old_array, array, left_subarray_size);

            print_array_proc(array, number_of_elements, my_id);

            // Send the right sublist list[j+1:nelement-1] to partner.
            MPI_Send(&old_array[left_subarray_size], right_subarray_size, MPI_INT, partner, 20, new_comm);
            // Receive the left sublist of partner.
            // Append the received list to my left list.
            MPI_Recv(&array[left_subarray_size], partners_sublist_size, MPI_INT, partner, 20, new_comm, &status);
            print_array_proc(array, number_of_elements, my_id);
            free(old_array);
        } else {
            // This is the right partner.
            printf("Right\n");

            // Send the size of left sublist to the partner, 
            // so it can allocate the buffer.
            MPI_Send(&left_subarray_size, 1, MPI_INT, partner, 10, new_comm);
            
            // Receive the size of the partner's right sublist.
            int partners_sublist_size;
            MPI_Recv(&partners_sublist_size, 1, MPI_INT, partner, 10, new_comm, &status);
            printf("Partners sublist size: %d\n", partners_sublist_size);
            
            // Prepare the buffer for receiving new data.
            number_of_elements = partners_sublist_size + right_subarray_size;

            int *old_array = array;
            array = (int*) calloc(number_of_elements, sizeof(int));
            copy(old_array + right_subarray_size, array + partners_sublist_size, right_subarray_size);

            print_array_proc(array, number_of_elements, my_id);
            
            // Send the left sublist list[0:j] to partner.
            MPI_Send(&old_array[0], left_subarray_size, MPI_INT, partner, 20, new_comm);
            // Receive the right sublist of partner.
            MPI_Recv(&array[0], partners_sublist_size, MPI_INT, partner, 20, new_comm, &status);
            // Append the received list to my right list.
            print_array_proc(array, number_of_elements, my_id);
            free(old_array);
        }

        mask >>= 1;
        bitvalue >>= 1;

        MPI_Comm_free(&new_comm);
    }

    qsort(array, number_of_elements, sizeof(int), compare);
    print_array_proc(array, number_of_elements, my_id);

    return number_of_elements;
}

int main(int argc, char* argv[]) {
    int number_of_processes;
    int my_id;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);

    printf("%d\n", my_id);
    int* array = generate_array(my_id);
    int number_of_elements = hyper_quick_sort(array, my_id, number_of_processes);

    MPI_Barrier(MPI_COMM_WORLD);

    // if (my_id != 0) {
    //     MPI_Send(array, number_of_elements, MPI_INT, 0, 1, MPI_COMM_WORLD);
    // }

    // int i;
    // print_array(array, number_of_elements);
    // printf(", ");
    // for (i = 1; i < number_of_processes; i++) {
        
    // }

    MPI_Finalize();

    return 0;
}
