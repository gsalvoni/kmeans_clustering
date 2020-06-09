#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <mpi.h>

#define nb_letters 26
#define nb_files 1000
#define increment_min 97 // ASCII value of a
#define increment_maj 65 // ASCII value of A
#define nb_iterations_max 100
#define threshold 0.00001

/*
 FUNCTIONS
 */

/*
 The function 'operations' measures the histogram of the occurrence and
 the average distance of each of the 26 English letters in the book *file
 and inserts the values in their corresponding arrays.
 */

void operations (FILE *file, float alphabet_distance[]){
    
    int character_nb = 0;
    int alphabet_occurrence[nb_letters] = {0};
    char str[150];
    
    // Matrix to store the first and the last times each letter appears in the book.
    int row;
    int first_last[nb_letters][2];
    // Initialisation of the matrix.
    for (row=0; row<nb_letters; row++){
        first_last[row][0] = 0;
        first_last[row][1] = 0;
    }
    
    // Book reading.
    int i;
    while(fgets(str, 150, file)!= NULL){
        for (i=0; i<strlen(str);  i++){
            // Each time a character is read, character_nb is incremented.
            character_nb++;
            
            /* Count the number of times each letter appears in the text
             and fill in the matrix first_last.
             */
            // If the letter is an uppercase.
            if(str[i]>(increment_maj-1) && str[i]<(increment_maj+26)){
                alphabet_occurrence[str[i]-increment_maj]++;
                
                if(first_last[str[i]-increment_maj][0] == 0)
                    first_last[str[i]-increment_maj][0] = character_nb;
                
                first_last[str[i]-increment_maj][1] = character_nb;
            }
            
            // If the letter is a lowercase.
            if(str[i]>(increment_min-1) && str[i]<(increment_min+26)){
                alphabet_occurrence[str[i]-increment_min]++;
                
                if(first_last[str[i]-increment_min][0] == 0)
                    first_last[str[i]-increment_min][0] = character_nb;
                
                first_last[str[i]-increment_min][1] = character_nb;
            }
        }
    }
    
    // Determination of the average distance for each letter.
    int k;
    for (k=0; k<nb_letters; k++){
        if (alphabet_occurrence[k] > 1)
            alphabet_distance[k] = (float)((first_last[k][1] - first_last[k][0]) - (alphabet_occurrence[k]-1))/(alphabet_occurrence[k]-1);
        // If a letter appears only once, its distance = the size of the whole text.
        else if (alphabet_occurrence[k] == 1)
            alphabet_distance[k] = character_nb;
        // If a letter is not in the text, its distance = 0.
        else if (alphabet_occurrence[k] == 0)
            alphabet_distance[k] = 0;
    }
}

/*
 The function 'create_file' creates a '.csv' file containing the book clustering
 among the different centroids.
 */

void create_file (int book_partition[], char book[]){
    FILE *file;
    
    file = fopen(book,"wb");
    
    int i;
    for (i=0; i<nb_files; i++){
        fprintf(file,"%d;%d",i+1,book_partition[i]);
        fprintf(file,"\n");
    }
    
    fclose(file);
}

/*
 The function 'read_book' reads a book and then computes the average distance for
 each letter.
 */

void read_book (char book[], float alphabet_distance[]){
    FILE *file;
    
    file = fopen(book,"r");
    
    // Check if the file is open.
    if (file == NULL){
        printf("File can't be read !\n");
    }
    
    operations(file, alphabet_distance);
    
    fclose(file);
    
}

/*
 The function 'euclidian_distance' measures the euclidian distance between a book and
 all the centroids. The output is the number of the centroid to which the book is the
 closest.
 */
int euclidian_distance (int nb_centroids, float centroids[nb_centroids][nb_letters], float alphabet_distance[nb_letters]){
    int i, j;
    float eucl_dist, tmp, min;
    int closest;
    
    for (j=0; j<nb_centroids; j++){
        eucl_dist = 0.0;
        tmp = 0.0;
        for (i=0; i<nb_letters; i++)
            tmp = tmp + (alphabet_distance[i] - centroids[j][i])*(alphabet_distance[i] - centroids[j][i]);
        eucl_dist = sqrtf(tmp);
        if (j==0){ // Initialisation of 'min' and 'closest' to the first centroid.
            min = eucl_dist;
            closest = 0;
        }
        else{
            if (eucl_dist <= min)
                closest = j;
        }
    }
    
    return closest;
}

/*
 The function 'new_centroids' adds the average distance of a book that have been
 previously associated to its nearest centroid. It also increments the number of
 books that are associated with this particular centroid.
 */
void new_centroids(float centroid_tmp[], float alphabet_distance[], int *nb_closest) {
    int i;
    (*nb_closest)++;
    for (i=0; i<nb_letters; i++)
        centroid_tmp[i] += alphabet_distance[i];
    
}

/*
 The function 'update_centroids' aims to move the data comprised in the new centroids
 back in the arrays of the initial centroids. This is to start the k-means algorithm
 with the new centroids that have been moved to the center of their associated books.
 */
void update_centroids(float centroid_init[], float centroid_new[]){
    int j;
    for (j=0; j<nb_letters; j++)
        centroid_init[j] = centroid_new[j];
}

/*
 MAIN
 */

int main(int argc, char *argv[]) {
    
    char book[100];
    int myrank, np, mytag = 21;
    int start_val, end_val;
    int i, j, k, l, check = 0;
    int nb_iterations, nb_centroids;
    
    /*
     Get the number of partitions from the command line. If nothing is
     specified, the number is set to 4 by default.
     */
    if (argc == 2){
        nb_centroids = atoi(argv[1]);
    }
    else{
        nb_centroids = 4;
    }
    
    /*
     Matrices initializations.
     */
    float alphabet_distance[nb_files][nb_letters];
    for (i=0; i<nb_files; i++){
        for (j=0; j<nb_letters; j++)
            alphabet_distance[i][j] = 0.0;
    }
    
    // Matrices that associate each book to its nearest centroid.
    int book_partition_tmp[nb_files] = {0};
    int book_partition[nb_files] = {0};
    
    // Matrices that store the centroids values.
    float centroids_init[nb_centroids][nb_letters];
    float centroids_new[nb_centroids][nb_letters];
    float centroids_tmp[nb_centroids][nb_letters];
    for (i=0; i<nb_centroids; i++){
        for (j=0; j<nb_letters; j++){
            centroids_init[i][j] = 0.0;
            centroids_new[i][j] = 0.0;
            centroids_tmp[i][j] = 0.0;
        }
    }
    
    // Matrices that contain the number of books associated to each centroid.
    int nb_closest[nb_centroids];
    int nb_closest_rank0[nb_centroids];
    for (i=0; i<nb_centroids; i++){
        nb_closest[i] = 0;
        nb_closest_rank0[i] = 0;
    }
    
    /*
     Initilization of MPI.
     */
    MPI_Status mystatus;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    
    /*
     Each processor reads the same number of books.
     */
    start_val = nb_files * myrank / np + 1;
    end_val = nb_files * (myrank+1) / np;
    
    for (i=start_val; i<=end_val; i++){
        sprintf(book,"book%01d.txt",i);
        read_book(book, alphabet_distance[i-1]);
    }
    
    /*
     On processor 0, random initialisation of the 4 centroids between its
     start value and its end value.
     */
    if (myrank == 0){
        for (i=0 ; i<nb_centroids; i++){
            for (j=0; j<nb_letters; j++)
                centroids_init[i][j] = alphabet_distance[(rand() % (end_val-1 + 1 - start_val-1)) + start_val-1][j];
        }
    }
    
    /*
     K-means algorithm.
     */
    for (k=0; k<=nb_iterations_max; k++){
       
        /*
         PROCESSOR 0.
         */
        if (myrank == 0) {
            
            for (i=1; i<np; i++){
                // Send the centroids and the nb_closest array to the other processors.
                MPI_Send(&centroids_init, nb_centroids*nb_letters, MPI_FLOAT, i, mytag, MPI_COMM_WORLD);
                MPI_Send(&nb_closest, nb_centroids, MPI_INT, i, mytag+1, MPI_COMM_WORLD);
            }
            
            // For each book : association to its nearest centroid.
            for (i=start_val-1; i<=end_val-1; i++){
                book_partition_tmp[i] = euclidian_distance(nb_centroids, centroids_init, alphabet_distance[i]);
                new_centroids(centroids_tmp[book_partition_tmp[i]], alphabet_distance[i], &nb_closest[book_partition_tmp[i]]);
            }
        
            // Compute the new centroids based on the positions of their respective associated books.
            MPI_Reduce(&centroids_tmp, &centroids_new, nb_centroids*nb_letters, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce(&nb_closest, &nb_closest_rank0, nb_centroids, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
            
            for (i = 0; i<nb_centroids; i++){
                for (j = 0; j<nb_letters; j++){
                    if (nb_closest_rank0[i] != 0)
                        centroids_new[i][j] /= nb_closest_rank0[i];
                }
            }
            
            // Check is the centroids have moved.
            for (i = 0; i<nb_centroids; i++){
                for (j = 0; j<nb_letters; j++){
                    if(fabs(centroids_new[i][j]-centroids_init[i][j]) > threshold){
                        // If the centroids are different : check = 1.
                        check = 1;
                        nb_iterations = k+1;
                        break;
                    }
                }
            }
            // Send the comparison result to the other processors to know what to do.
            if (check == 0){
                // If check = 0 : stop the algorithm.
                for (l=1; l<np; l++)
                    MPI_Send(&check, 1, MPI_INT, l, mytag-1, MPI_COMM_WORLD);
                break;
            }
            if (check == 1){
                // If check = 1 : new iteration, update the centroids.
                for (l=1; l<np; l++)
                    MPI_Send(&check, 1, MPI_INT, l, mytag-1, MPI_COMM_WORLD);
                
                check = 0;
                
                for (i = 0; i<nb_centroids; i++)
                    update_centroids(centroids_init[i], centroids_new[i]);
            
                for (i=0; i<nb_centroids; i++){
                    nb_closest[i] = 0;
                    nb_closest_rank0[i] = 0;
                }
                for (i=0; i<nb_centroids; i++){
                    for (j=0; j<nb_letters; j++){
                        centroids_tmp[i][j] = 0.0;
                        centroids_new[i][j] = 0.0;
                    }
                }
            }
        }
        /*
         OTHER PROCESSORS.
         */
        else {
            
            // Receive the centroids and the nb_closest array from processor 0.
            MPI_Recv(&centroids_init, nb_centroids*nb_letters, MPI_FLOAT, 0, mytag, MPI_COMM_WORLD, &mystatus);
            MPI_Recv(&nb_closest, nb_centroids, MPI_INT, 0, mytag+1, MPI_COMM_WORLD, &mystatus);
            
            // For each book : association to its nearest centroid.
            for (i=start_val-1; i<=end_val-1; i++){
                book_partition_tmp[i] = euclidian_distance(nb_centroids, centroids_init, alphabet_distance[i]);
                new_centroids(centroids_tmp[book_partition_tmp[i]], alphabet_distance[i], &nb_closest[book_partition_tmp[i]]);
            }
            
            // Send the computed centroids to processor 0 so that it computes the new ones.
            MPI_Reduce(&centroids_tmp, &centroids_new, nb_centroids*nb_letters, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce(&nb_closest, &nb_closest_rank0, nb_centroids, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
            
            // Receive the check variable for the processor 0.
            MPI_Recv(&check, 1, MPI_INT, 0, mytag-1, MPI_COMM_WORLD, &mystatus);
            
            // If check = 0 : stop the algorithm.
            if (check == 0){
                nb_iterations = k+1;
                break;
            }
            // If check = 1 : new iteration.
            if (check == 1){
                for (i=0; i<nb_centroids; i++){
                    for (j=0; j<nb_letters; j++)
                        centroids_tmp[i][j] = 0.0;
                }
            }
        }
    
    }
    
    // Regroup all the partition numbers in a vector on the processor 0.
    MPI_Reduce(&book_partition_tmp, &book_partition, nb_files, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    
    // Create the '.csv' file.
    if (myrank == 0){
        sprintf(book,"output.csv");
        create_file(book_partition,book);
    }
    if (myrank == 0)
        printf("Number of iterations : %d\n",nb_iterations);
    
    MPI_Finalize();
    return 0;
}



