##################################################################
## lostruc across the genome for labrador charr
##
## Matt Brachmann (PhDMattyB)
##
## 2021-03-25
##
##################################################################

## Load the data manipulation work horse
library(tidyverse)
library(janitor)
library(data.table)
library(lostruct)

## set the working directory for this script
setwd('~/Charr_Adaptive_Introgression/Charr_Project_1/GeneticData/')

## load in the map file for the SNP array data set
map = read_tsv('Charr_Lab_recode12_25.03.2021.map', 
               col_names = c('Chromosome', 
                             'MarkerID', 
                             'Genetic_dist', 
                             'Physical_dist'))

## load in the ped file for the SNP array data set
ped = read_table2('Charr_Lab_recode12_25.03.2021.ped', 
                  col_names = c('FamilyID', 
                                'IndividualID', 
                                'PaternalID', 
                                'MaternalID', 
                                'Sex', 
                                'Phenotype', 
                                map$MarkerID))

Create_tped = function(ped, map){
  ## Obtaining a vector of names for each individual
  message('Creating vector of individual names')
   indiv_names = ped %>% 
    filter(FamilyID != 'GDL') %>%
    select(IndividualID) %>% 
    t() %>% 
    as_tibble() %>% 
    row_to_names(row_number = 1) %>% 
    names()
   
   message('Merging ped and map files')
  ped %>% 
    filter(FamilyID != 'GDL') %>%
    select(contains('AX-')) %>% 
    t() %>% 
    as_tibble() %>% 
    rename_all(funs(c(indiv_names))) %>% 
    bind_cols(map) %>% 
    select(Chromosome, 
           MarkerID, 
           Genetic_dist, 
           Physical_dist, 
           everything())
}

tped = Create_tped(ped = ped, 
                   map = map) 

# Chr by Chr functions ----------------------------------------------------

## The map functions aren't necessarily working to 
## get the results up and running. We ran into a road block
## and had to get the results on a chromosome by chromosome basis
## we might as well create a function that runs well for a chr
## and outputs the results cleanly. 
## I feel like we learned a bit using the map functions
## we also need to include map data into this so we 
## actually know where we are in the genome
## the map functions didn't really get that done


## Need a matrix of genotypes for the lostruct functions
## individuals are in columns and markers are in rows


lostruct_run = function(data, 
                        chr, 
                        window_size, 
                        k_value){
  df = data %>% 
    dplyr::select(Chromosome, 
           5:length(tped)) %>% 
    filter(Chromosome == chr) %>% 
    dplyr::select(-Chromosome)
  
  message('calculating eigenvectors for windows')
  eigen = eigen_windows(df, 
                win = window_size, 
                k = k_value)
  message('calculating distance matrix')
  windist = pc_dist(eigen, 
          npc = k_value) %>% 
    as_tibble()
  
  message('calculating mean window size for the chr')
  window_data = data %>% 
    select(1:4) %>% 
    filter(Chromosome == 1) %>% 
    mutate(window = ceiling(row_number()/window_size)) %>% 
    group_by(window) %>% 
    mutate(mean_window = mean(Physical_dist)) %>% 
    distinct(mean_window, 
             .keep_all = T) %>% 
    filter(window %in% 1:nrow(windist))
  
  combo_data = bind_cols(window_data, 
                         windist)
  
  message('Scaling the data using MDS')
  MDS_data = cmdscale(combo_data[7:length(combo_data)], 
           eig = TRUE, 
           k = k_value)
  
  MDS_points = MDS_data$points %>% 
    as_tibble()
  
  combo_data = bind_cols(combo_data, 
                         MDS_points) %>% 
    rename(MDS_Points1 = V1...39,
           MDS_Points2 = V2...40,
           V1 = V1...7,
           V2 = V2...8) %>% 
    dplyr::select(-MarkerID, 
                  -Genetic_dist, 
                  -Physical_dist, 
                  Chromosome, 
                  mean_window, 
                  window:MDS_Points2)
  # output = list(combo_data, 
  #               MDS_data)
  # 
  # return(output)
}


lostruct_data = lostruct_run(data = tped, 
             chr = 1, 
             window_size = 20, 
             k_value = 2)

MDS_scatter = function(data){
  theme_set(theme_bw())
  
  MDS_points = data %>% 
    dplyr::select(MDS_Points1, 
                  MDS_Points2)
  
  window_distance = data %>%
    dplyr::select(contains('V'))
  
  MDS_points %>% 
    ggplot(aes(x = MDS_Points1, 
           y = MDS_Points2,
           col = rainbow(nrow(window_distance))))+
    geom_point(size = 3) +
    labs(x = 'MDS coordinate 1', 
         y = 'MDS coordinate 2')+
    theme(
      legend.position = 'none', 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      axis.title = element_text(size = 14), 
      axis.text = element_text(size = 12)
    )
  
}

MDS_scatter(lostruct_data)

MDS_outliers = function(data){
  MDS1_outliers = data %>% 
    ungroup() %>% 
    mutate(MDS1_cutoff = mean(MDS_Points1)+(2*sd(MDS_Points1))) %>% 
    filter(abs(MDS_Points1) > MDS1_cutoff)
  MDS2_outliers = data %>% 
    ungroup() %>% 
    mutate(MDS2_cutoff = mean(MDS_Points2)+(2*sd(MDS_Points2))) %>% 
    filter(abs(MDS_Points2) > MDS2_cutoff)
  
  
  
}

MDS1_outliers = lostruct_data %>% 
  ungroup() %>% 
  mutate(MDS1_cutoff = mean(MDS_Points1)+(2*sd(MDS_Points1))) %>% 
  filter(abs(MDS_Points1) > MDS1_cutoff) %>% 
  select(Chromosome, 
         window, 
         mean_window, 
         MDS_Points1, 
         MDS1_cutoff)

MDS2_outliers = lostruct_data %>% 
  ungroup() %>% 
  mutate(MDS2_cutoff = mean(MDS_Points2)+(2*sd(MDS_Points2))) %>% 
  filter(abs(MDS_Points2) > MDS2_cutoff) %>% 
  select(Chromosome, 
         window, 
         mean_window, 
         MDS_Points2, 
         MDS2_cutoff)

  
         



chrom_mds1_outliers=chrom_results[which(chrom_results$MDS1 > mean(chrom_results$MDS1)+2*sd(chrom_results$MDS1) | chrom_results$MDS1 < mean(chrom_results$MDS1)+-2*sd(chrom_results$MDS1)),]
chrom_mds2_outliers=chrom_results[which(chrom_results$MDS2 > mean(chrom_results$MDS2)+2*sd(chrom_results$MDS2) | chrom_results$MDS2 < mean(chrom_results$MDS2)+-2*sd(chrom_results$MDS2)),]

#If want values close to zero
chrom_zeros=chrom_results[which(chrom_results$MDS2 > -0.01  & chrom_results$MDS1 > -0.01 & chrom_results$MDS2 < 0.01  & chrom_results$MDS1 < 0.01 ),]


##

# Everything below this does not work fully and is essentially tri --------

# map function variation --------------------------------------------------

## code to split by chrs and  apply lostruc functions per chr
## are currently in flux

## map function methods

by_chr = tped %>% 
  select(Chromosome, 
         5:length(tped)) %>% 
  group_by(Chromosome) %>% 
  nest()

# by_chr$data[[24]]

cal_eigen_win = function(df){
  eigen_windows(df, 
                win = 20, 
                k = 2)
}

## HOLY SHIT IT'S ALIVE!!!!!!!
by_chr = by_chr %>%
  mutate(eigen_win = map(data, cal_eigen_win))

cal_windist = function(df){
  pc_dist(df, 
          npc = 2)
}

by_chr = by_chr %>% 
  mutate(windist = map(eigen_win, cal_windist))

cal_fit2d = function(df){
  cmdscale(df, 
           eig = TRUE, 
           k = 2)
}

by_chr = by_chr %>% 
  mutate(fit2d = map(windist, 
                     cal_fit2d))

by_chr = by_chr %>% 
  unnest_wider(fit2d)


lostruct_points = function(data, chr){
  data %>% 
    # chop(points)
    unnest(points) %>% 
    select(Chromosome, 
           points) %>%
    ungroup() %>%
    filter(Chromosome == chr)%>%
    select(-Chromosome) %>% 
    write_tsv('Lostruct_garbage_test.tsv')

} 

lostruct_windist = function(data, chr){
  data %>%
    select(Chromosome,
           windist) %>% 
    filter(Chromosome == chr) %>% 
    unnest(windist) %>% 
    ungroup() %>% 
    select(-Chromosome)
}


# lostruct processing post calculation ------------------------------------

# chrs = by_chr %>%
#   select(Chromosome) %>%
#   distinct() %>%
#   ungroup() %>%
#   as.list()


chr1_points = lostruct_points(by_chr, 1)

## Need to slice at the halfway point of the df
## the numbers in the slice will vary
## read the tables in without putting it into memory
points1 = read_table2('Lostruct_garbage_test.tsv', 
                      col_names = 'points1') %>% 
  slice(2:33)
## slice the data frame from halfway to the end
points2 = read_table2('Lostruct_garbage_test.tsv', 
                      col_names = 'points2') %>% 
  slice(34:65)
chr_points = bind_cols(points1, 
          points2) 
## write the final tsv file with the points for the chromosome
write_tsv(chr_points, 
          'chr1_CharrArray_win20_MSDS_points.tsv')

  
lostruct_windist(by_chr, 1) %>% 
## I usually hate base R but, damn did that just come in clutch
  write.table('chr1_CharrArray_win20_windist.txt')



##
# plotting -------------------------------------------------------------------

theme_set(theme_bw())

chr1_points = read_tsv('chr1_CharrArray_win20_MSDS_points.tsv')
## read in the matrix of the distances between windows
chr1_windist = read.table('chr1_CharrArray_win20_windist.txt') %>% 
  as_tibble()

full_df = bind_cols(chr1_points, 
                    chr1_windist) %>% 
  write_tsv('chr1_CharrArray_win20_fulldf.tsv')

# plot(chr1_points, 
#      xlab = 'Coordinate 1', 
#      ylab = 'Coordinate 2', 
#      col = rainbow(1.2*nrow(chr1_windist)), 
#      pch = 19)


ggplot(data = chr1_points, 
       aes(x = points1, 
           y = points2, 
           col = rainbow(nrow(chr1_points))))+
  geom_point(size = 3) +
  labs(x = 'Coordinate 1', 
       y = 'Coordinate 2', 
       title = 'MSDS plot for each PCA window Chr 1')+
  theme(
    legend.position = 'none', 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.title = element_text(size = 14), 
    axis.text = element_text(size = 12)
  )


# Need to define outliers -------------------------------------------------

chr1_fulldf = read_tsv('chr1_CharrArray_win20_fulldf.tsv')

chr1_mds1_outlier = chr1_fulldf[which(chr1_fulldf$points1 > mean(chr1_fulldf$points1)+ 2*sd(chr1_fulldf$points1) | chr1_fulldf$points1 < mean(chr1_fulldf$points1)+-2*sd(chr1_fulldf$points1)),]
chr1_mds2_outlier = chr1_fulldf[which(chr1_fulldf$points2 > mean(chr1_fulldf$points2)+ 2*sd(chr1_fulldf$points2) | chr1_fulldf$points2 < mean(chr1_fulldf$points2)+-2*sd(chr1_fulldf$points2)),]




