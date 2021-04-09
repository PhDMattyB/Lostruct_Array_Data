# Lostruct_Array_Data
Created some handy add on functions to run the lostruct package more efficiently on array data

There are currently a few different functions that should make running the lostruct package a little bit easier on array data. The input files need to be in plink format and both the ped and the map files are needed.

Create_tped: Combines the ped and map files into a single dataframe and then transposes the dataframe to put it into a format that lostruct needs. Lostruct need a matrix or dataframe of [i,j] where i = markers and are located in rows and j = individual genotypes located in columns.

lostruct_run: Run the three main lostruct functions on the genotypic information provided by the tped. This function also calculated the mean position of each window determined using the window_size argument. 

MDS_survey: A basic MDS plot to survey and see what the data coming out of the lostruct package looks like. Nothing fancy. 

Outlier_hunter: Identified PC outliers according to the MDS. The outliers identified on MDS 1 and 2 are outside of +/- two standard deviations from the mean. The data set gets a new column called outlier_lab which identifies outliers on MDS1 and 2. This will come in clutch when plotting a new version of the MDS graphs which will identify the outlier PCs as well as when we plot this along the genome for each window. 

Outlier_plots: Plots three different plots in one (thanks patchwork). The first plot is just the MDS plot with the outlier windows coloured according to which MDS axis they belong two. The next two plots show the MDS score along the genome according to the mean window position of each of the windows. The two plots are similar but differ in which MDS axis the outlier lie on. 

This repository is actively being worked on and some new things should be out soon to make things even more easier. 
