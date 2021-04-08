# Lostruct_Array_Data
Creating some add on functions to run the lostruct package more efficiently on array data

Currently there are two functions that should make running the lostruct package a little bit easier on array data. The input files need to be in plink format and both the ped and the map files are needed.

Create_tped: Combines the ped and map files into a single dataframe and then transposes the dataframe to put it into a format that lostruct needs. Lostruct need a matrix or dataframe of [i,j] where i = markers and are located in rows and j = individual genotypes located in columns.

lostruct_run: Run the three main lostruct functions on the genotypic information provided by the tped. This function also calculated the mean position of each window determined using the window_size argument. 

I will also add some plotting functions to make plotting the data and further cleaning of the data a bit easier. This repository is actively being worked on. 
