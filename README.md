# Gutenberg Project - Book Clustering
A parallel program implementing a K-means algorithm for book clustering.

## Goal
Read an archive of books, clean it up using K-means algorithm, implement the code on multiple processors, and find those from the Gutenberg Project.

## Technologies
Project created with:
* C
* SLURM 
* MPI

## Method
Measure the Euclidian distance between each book and the centroids, associate books to the nearest centroid, move centroids, and repeat until convergence. Note: each processor reads the same number of books.
