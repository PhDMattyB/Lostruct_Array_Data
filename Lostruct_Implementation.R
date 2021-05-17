##################################################################
## Lostruct functions implementation
##
## Matt Brachmann (PhDMattyB)
##
## 2021-05-12
##
##################################################################

## Load the packages needed to run the functions
library(tidyverse)
library(janitor)
library(data.table)
library(patchwork)
library(lostruct)
library(adegenet)
library(viridis)

# Lostruct run ------------------------------------------------------------
## set the working directory for this script
setwd('~/Charr_Adaptive_Introgression/Charr_Project_1/GeneticData/')

## Environmental metadata
env_data = read_csv('~/Charr_Adaptive_Introgression/Charr_Project_1/SampleSiteData/SampleSites_Coords_1June2020.csv')

## load in the map file for the SNP array data set
# map = read_tsv('Charr_Lab_recode12_25.03.2021.map', 
#                col_names = c('Chromosome', 
#                              'MarkerID', 
#                              'Genetic_dist', 
#                              'Physical_dist'))

map = read_tsv('~/Charr_Adaptive_Introgression/Charr_Project_1/GeneticData/Charr_All_Pops_recode12_ChrConvert.map') %>% 
  rename(Chromosome = CHR)

# load in the ped file for the SNP array data set
# OG_ped = read_table2('Charr_Lab_recode12_25.03.2021.ped',
#                      col_names = c('FamilyID',
#                                    'IndividualID',
#                                    'PaternalID',
#                                    'MaternalID',
#                                    'Sex',
#                                    'Phenotype',
#                                    map$MarkerID))

OG_ped = read_table2('~/Charr_Adaptive_Introgression/Charr_Project_1/GeneticData/Charr_All_Pops_recode12.ped', 
                     col_names = c('#FamilyID', 
                                   'IndividualID', 
                                   'ParentalID', 
                                   'MaternalID', 
                                   'Sex', 
                                   'Phenotype', 
                                   map$MarkerID)) %>% 
  filter(`#FamilyID` %in% c('ANA', 
                            'AVA', 
                            'BLD', 
                            'BRG', 
                            'ENG', 
                            'FRD', 
                            'FRN', 
                            'FRS', 
                            'GDL', 
                            'IGL', 
                            'IKA', 
                            'IKI', 
                            'IKL', 
                            'KAM', 
                            'KAN', 
                            'KIN', 
                            'KIY', 
                            'KOG', 
                            'KOM', 
                            'MBB', 
                            'MCC', 
                            'NAC', 
                            'NAK', 
                            'NOR', 
                            'PAL', 
                            'PAN', 
                            'PBP', 
                            'PUT', 
                            'R103', 
                            'R104', 
                            'R105', 
                            'R109', 
                            'R110', 
                            'R78', 
                            'REI', 
                            'STC', 
                            'SWA', 
                            'TOR', 
                            'UNH'))

tped = Create_tped(ped = OG_ped, 
                   map = map)
# Chr 1 -------------------------------------------------------------------
lostruct_data = lostruct_run(data = tped, 
                             chr = 1, 
                             window_size = 20, 
                             k_value = 2)

## If you get the error that the data frame isn't numeric
## use the two select_if arguments to hunt down the 
## column that isn't numeric
## Hunting down hidden non-numeric columns
# df %>% select_if(is.numeric)
# df %>% select_if(negate(is.numeric))

# MDS_survey(lostruct_data)

outliers = Outlier_hunter(data = lostruct_data,
                          sd_percentile = 2)

# Outlier_plots(normal_data = lostruct_data,
#               outlier_data = outliers)

outlier_full_data = Outlier_data(data = tped, 
                                 outlier_data = outliers, 
                                 chr = 1, 
                                 window_size = 20, 
                                 k_value = 2)


outlier_full_data

## CHanging the populations that are used in the lostruct
## analysis DRASTICALLY changes which regions are
## flagged as outliers and potential structural variants

## Identify which regions might actually show a SV
## Cycle through each outlier region through the four functions
## below. 
Chr1_map_win10 = map_maker(outlier_full_data$'10')
Chr1_ped_win10 = ped_maker(outlier_full_data$'10')


Chr1_win_data = Adegenet_PCA(outlier_ped = Chr1_ped_win, 
                              outlier_map = Chr1_map_win, 
                              OG_ped = OG_ped,
                              env = env_data)

# write_tsv(Chr1_win7_data, 
#           '~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/CHR1_REGION1_lostruct.txt')

Pop_that_pca(Chr1_win_data, 
                             pop_num = 39,
                             chr_num = 1, 
                             win_num = 11)


# CHR1 SV REGION HUNTING --------------------------------------------------
## Combo data for continuous regions
## CHR1 REGION1 WINDOWS 10 & 11
chr_map_win10= map_maker(outlier_full_data$'10')
chr_ped_win10 = ped_maker(outlier_full_data$'10')
chr_data_win10 = Adegenet_PCA(outlier_ped = chr_ped_win10, 
                              outlier_map = chr_map_win10, 
                              OG_ped = OG_ped,
                              env = env_data)

chr_map_win11 = map_maker(outlier_full_data$'11')
chr_ped_win11 = ped_maker(outlier_full_data$'11')
chr_data_win11 = Adegenet_PCA(outlier_ped = chr_ped_win11, 
                              outlier_map = chr_map_win11, 
                              OG_ped = OG_ped,
                              env = env_data)
chr_data_win11 = chr_data_win11 %>% 
  select(contains('AX-'))

chr1_combo_win1011 = bind_cols(chr_data_win10, 
                               chr_data_win11)

PCA_outlier_wins1011 = Pop_that_pca(chr1_combo_win1011, 
                                    pop_num = 39,
                                    chr_num = 1, 
                                    win_num = 1617)

write_tsv(chr1_combo_win1011, 
          '~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/Lostruct_chr1_region1.txt')
## write map files for the combo region

chr1_region1_map = bind_rows(chr_map_win10, 
                                     chr_map_win11) %>% 
  select(1:4) %>% 
  rename(`#Chromosome` = Chromosome, 
         `Marker ID` = MarkerID, 
         `Genetic distance` = Genetic_dist, 
         `Physical distance` = Physical_dist) %>% 
  write_tsv('~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/Chr1_region1_outlier_14.05.2021.map')

head(chr1_region1_map)
tail(chr1_region1_map)

## Calculate regions size
region_size = (23754321-21338526)/1000000
## CHR1 REGION1 is 2.42 Mb

## CHR 1 REGION 2 WINDOWS 15 & 16
chr_map_win15 = map_maker(outlier_full_data$'15')
chr_ped_win15 = ped_maker(outlier_full_data$'15')
chr_data_win15 = Adegenet_PCA(outlier_ped = chr_ped_win15, 
                              outlier_map = chr_map_win15, 
                              OG_ped = OG_ped,
                              env = env_data)

chr_map_win16 = map_maker(outlier_full_data$'16')
chr_ped_win16 = ped_maker(outlier_full_data$'16')
chr_data_win16 = Adegenet_PCA(outlier_ped = chr_ped_win16, 
                              outlier_map = chr_map_win16, 
                              OG_ped = OG_ped,
                              env = env_data)



chr_data_win16 = chr_data_win16 %>% 
  select(contains('AX-'))

chr1_combo_win1516 = bind_cols(chr_data_win15, 
                               chr_data_win16)

Pop_that_pca(chr1_combo_win1516, 
                                    pop_num = 39,
                                    chr_num = 1, 
                                    win_num = 1617)

write_tsv(chr1_combo_win1516, 
          '~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/Lostruct_chr1_region2.txt')
## write map files for the combo region

chr1_region2_map = bind_rows(chr_map_win15, 
                                     chr_map_win16) %>% 
  select(1:4) %>% 
  rename(`#Chromosome` = Chromosome, 
         `Marker ID` = MarkerID, 
         `Genetic distance` = Genetic_dist, 
         `Physical distance` = Physical_dist) %>% 
  write_tsv('~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/Chr1_region2_windows_14.05.2021.map')

tail(chr1_region2_map)
head(chr1_region2_map)

## Calculate regions size
chr1_region2_size = (31714898-28470108)/100000
## The region size of the outlier is
## 32.45Mb with 40 SNPs. 


## To make the combo map ped file
## Need to use the data object in the Data_explore script
## The data needed for adegenet is recoded as 12 format
## needs to be in AA CC GG TT format


## CHR 1 REGION 3 WINDOWS 22 & 23
chr_map_win22 = map_maker(outlier_full_data$'22')
chr_ped_win22 = ped_maker(outlier_full_data$'22')
chr_data_win22 = Adegenet_PCA(outlier_ped = chr_ped_win22, 
                              outlier_map = chr_map_win22, 
                              OG_ped = OG_ped,
                              env = env_data)

chr_map_win23 = map_maker(outlier_full_data$'23')
chr_ped_win23 = ped_maker(outlier_full_data$'23')
chr_data_win23 = Adegenet_PCA(outlier_ped = chr_ped_win23, 
                              outlier_map = chr_map_win23, 
                              OG_ped = OG_ped,
                              env = env_data)



chr_data_win23 = chr_data_win23 %>% 
  select(contains('AX-'))

chr1_combo_win2223 = bind_cols(chr_data_win22, 
                               chr_data_win23)

Pop_that_pca(chr1_combo_win2223, 
             pop_num = 39,
             chr_num = 1, 
             win_num = 1617)

write_tsv(chr1_combo_win2223, 
          '~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/Lostruct_chr1_region3.txt')
## write map files for the combo region

chr1_region3_map = bind_rows(chr_map_win22, 
                             chr_map_win23) %>% 
  select(1:4) %>% 
  rename(`#Chromosome` = Chromosome, 
         `Marker ID` = MarkerID, 
         `Genetic distance` = Genetic_dist, 
         `Physical distance` = Physical_dist) %>% 
  write_tsv('~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/Chr1_region3_windows_14.05.2021.map')

tail(chr1_region3_map)
head(chr1_region3_map)

## Calculate regions size
chr1_region3_size = (39740732-38247425)/100000
## The region size of the outlier is
## 14.93Mb with 40 SNPs. 


## To make the combo map ped file
## Need to use the data object in the Data_explore script
## The data needed for adegenet is recoded as 12 format
## needs to be in AA CC GG TT format

# CHR 2 -------------------------------------------------------------------
lostruct_data = lostruct_run(data = tped, 
                             chr = 2, 
                             window_size = 20, 
                             k_value = 2)
outliers = Outlier_hunter(data = lostruct_data,
                          sd_percentile = 2)

outlier_full_data = Outlier_data(data = tped, 
                                 outlier_data = outliers, 
                                 chr = 2, 
                                 window_size = 20, 
                                 k_value = 2)




Chr2_map = map_maker(outlier_full_data$'16')
Chr2_ped = ped_maker(outlier_full_data$'16')


Chr2_data = Adegenet_PCA(outlier_ped = Chr2_ped, 
                               outlier_map = Chr2_map, 
                               OG_ped = OG_ped,
                               env = env_data)

# write_tsv(Chr2_win23_data,
#           '~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/CHR2_REGION3_win23_lostruct.txt')

Pop_that_pca(Chr2_data, 
             pop_num = 39,
             chr_num = 2, 
             win_num = 00)



# CHR2 REGIONS ------------------------------------------------------------

## COMBINE CONTINUOUS REGIONS 9, 10, 11
Chr2_map_win9 = map_maker(outlier_full_data$'9')
Chr2_ped_win9 = ped_maker(outlier_full_data$'9')
Chr2_map_win10 = map_maker(outlier_full_data$'10')
Chr2_ped_win10 = ped_maker(outlier_full_data$'10')
Chr2_map_win11 = map_maker(outlier_full_data$'11')
Chr2_ped_win11 = ped_maker(outlier_full_data$'11')

Chr2_data_win9 = Adegenet_PCA(outlier_ped = Chr2_ped_win9, 
                              outlier_map = CHr2_map_win9, 
                              OG_ped = OG_ped,
                              env = env_data)

Chr2_data_win10 = Adegenet_PCA(outlier_ped = Chr2_ped_win10, 
                              outlier_map = Chr2_map_win10, 
                              OG_ped = OG_ped,
                              env = env_data)
Chr2_data_win11 = Adegenet_PCA(outlier_ped = Chr2_ped_win11, 
                              outlier_map = Chr2_map_win11, 
                              OG_ped = OG_ped,
                              env = env_data)


Chr2_data_win10 = Chr2_data_win10 %>% 
  select(contains('AX-'))
Chr2_data_win11 = Chr2_data_win11 %>% 
  select(contains('AX-'))


Chr2_win91011 = bind_cols(Chr2_data_win9, 
                         Chr2_data_win10, 
                         Chr2_data_win11)

Pop_that_pca(Chr2_win91011, 
             pop_num = 39)

write_tsv(Chr2_win91011, 
          '~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/Lostruct_CHR2_REGION1_win91011_14.05.2021.txt')

## map file for continuous regions
Chr2_region1 = bind_rows(Chr2_map_win9, 
          Chr2_map_win10, 
          Chr2_map_win11)%>% 
  select(1:4) %>% 
  rename(`#Chromosome` = Chromosome) %>% 
  write_tsv('~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/CHR2_REGION1_14.05.2021.map')

head(Chr2_region1)
tail(Chr2_region1)

## Calculate regions size
region_size = (16729756-13964830)/1000000
## The region size of the outlier is
## 2.76 Megabases with 40 SNPs.

## To make the combo map ped file
## Need to use the data object in the Data_explore script
## The data needed for adegenet is recoded as 12 format
## needs to be in AA CC GG TT format

## COMBINE CONTINUOUS REGIONS 14, 15, 16
Chr2_map_win14 = map_maker(outlier_full_data$'14')
Chr2_ped_win14 = ped_maker(outlier_full_data$'14')
Chr2_map_win15 = map_maker(outlier_full_data$'15')
Chr2_ped_win15 = ped_maker(outlier_full_data$'15')
Chr2_map_win16 = map_maker(outlier_full_data$'16')
Chr2_ped_win16 = ped_maker(outlier_full_data$'16')

Chr2_data_win14 = Adegenet_PCA(outlier_ped = Chr2_ped_win14, 
                              outlier_map = CHr2_map_win14, 
                              OG_ped = OG_ped,
                              env = env_data)

Chr2_data_win15 = Adegenet_PCA(outlier_ped = Chr2_ped_win15, 
                               outlier_map = Chr2_map_win15, 
                               OG_ped = OG_ped,
                               env = env_data)
Chr2_data_win16 = Adegenet_PCA(outlier_ped = Chr2_ped_win16, 
                               outlier_map = Chr2_map_win16, 
                               OG_ped = OG_ped,
                               env = env_data)


Chr2_data_win15 = Chr2_data_win15 %>% 
  select(contains('AX-'))
Chr2_data_win16 = Chr2_data_win16 %>% 
  select(contains('AX-'))


Chr2_win141516 = bind_cols(Chr2_data_win14, 
                          Chr2_data_win15, 
                          Chr2_data_win16)

Pop_that_pca(Chr2_win141516, 
             pop_num = 39)

write_tsv(Chr2_win141516, 
          '~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/Lostruct_CHR2_REGION2_win141516_14.05.2021.txt')

## map file for continuous regions
Chr2_region2 = bind_rows(Chr2_map_win14, 
                         Chr2_map_win15, 
                         Chr2_map_win16)%>% 
  select(1:4) %>% 
  rename(`#Chromosome` = Chromosome) %>% 
  write_tsv('~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/CHR2_REGION2_14.05.2021.map')

head(Chr2_region2)
tail(Chr2_region2)

## Calculate regions size
region_size = (24291345-19687700)/1000000
## The region size of the outlier is
## 4.60 Megabases with 40 SNPs.

## To make the combo map ped file
## Need to use the data object in the Data_explore script
## The data needed for adegenet is recoded as 12 format
## needs to be in AA CC GG TT format


# CHR 3 -------------------------------------------------------------------
lostruct_data = lostruct_run(data = tped, 
                             chr = 3, 
                             window_size = 20, 
                             k_value = 2)
outliers = Outlier_hunter(data = lostruct_data,
                          sd_percentile = 2)

outlier_full_data = Outlier_data(data = tped, 
                                 outlier_data = outliers, 
                                 chr = 3, 
                                 window_size = 20, 
                                 k_value = 2)




Chr3_map = map_maker(outlier_full_data$'17')
Chr3_ped = ped_maker(outlier_full_data$'17')


Chr3_data = Adegenet_PCA(outlier_ped = Chr3_ped, 
                         outlier_map = Chr3_map, 
                         OG_ped = OG_ped,
                         env = env_data)

Pop_that_pca(Chr3_data, 
             pop_num = 39,
             chr_num = 3, 
             win_num = 00)


# CHR3 REGIONS ------------------------------------------------------------
## COMBINE CONTINUOUS REGIONS 9 & 10
Chr3_map_win9 = map_maker(outlier_full_data$'9')
Chr3_ped_win9 = ped_maker(outlier_full_data$'9')
Chr3_map_win10 = map_maker(outlier_full_data$'10')
Chr3_ped_win10 = ped_maker(outlier_full_data$'10')

Chr3_data_win9 = Adegenet_PCA(outlier_ped = Chr3_ped_win9, 
                              outlier_map = CHr3_map_win9, 
                              OG_ped = OG_ped,
                              env = env_data)

Chr3_data_win10 = Adegenet_PCA(outlier_ped = Chr3_ped_win10, 
                               outlier_map = Chr3_map_win10, 
                               OG_ped = OG_ped,
                               env = env_data)

Chr3_data_win10 = Chr3_data_win10 %>% 
  select(contains('AX-'))

Chr3_win910 = bind_cols(Chr3_data_win9, 
                          Chr3_data_win10)

Pop_that_pca(Chr3_win910, 
             pop_num = 39)

write_tsv(Chr3_win910, 
          '~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/Lostruct_CHR3_REGION1_win910_14.05.2021.txt')

## map file for continuous regions
Chr3_region1 = bind_rows(Chr3_map_win9, 
                         Chr3_map_win10)%>% 
  select(1:4) %>% 
  rename(`#Chromosome` = Chromosome) %>% 
  write_tsv('~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/CHR3_REGION1_14.05.2021.map')

head(Chr3_region1)
tail(Chr3_region1)

## Calculate regions size
(23582676-21700050)/1000000
## The region size of the outlier is
## 1.88 Megabases with 40 SNPs.

## To make the combo map ped file
## Need to use the data object in the Data_explore script
## The data needed for adegenet is recoded as 12 format
## needs to be in AA CC GG TT format

## COMBINE CONTINUOUS REGIONS 16 & 17
Chr3_map_win16 = map_maker(outlier_full_data$'16')
Chr3_ped_win16 = ped_maker(outlier_full_data$'16')
Chr3_map_win17 = map_maker(outlier_full_data$'17')
Chr3_ped_win17 = ped_maker(outlier_full_data$'17')

Chr3_data_win16 = Adegenet_PCA(outlier_ped = Chr3_ped_win16, 
                              outlier_map = CHr3_map_win16, 
                              OG_ped = OG_ped,
                              env = env_data)

Chr3_data_win17 = Adegenet_PCA(outlier_ped = Chr3_ped_win17, 
                               outlier_map = Chr3_map_win17, 
                               OG_ped = OG_ped,
                               env = env_data)

Chr3_data_win17 = Chr3_data_win17 %>% 
  select(contains('AX-'))

Chr3_win1617 = bind_cols(Chr3_data_win16, 
                        Chr3_data_win17)

Pop_that_pca(Chr3_win1617, 
             pop_num = 39)

write_tsv(Chr3_win1617, 
          '~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/Lostruct_CHR3_REGION1_win1617_14.05.2021.txt')

## map file for continuous regions
Chr3_region2 = bind_rows(Chr3_map_win16, 
                         Chr3_map_win17)%>% 
  select(1:4) %>% 
  rename(`#Chromosome` = Chromosome) %>% 
  write_tsv('~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/CHR3_REGION2_14.05.2021.map')

head(Chr3_region2)
tail(Chr3_region2)

## Calculate regions size
(33706958-30788823)/1000000
## The region size of the outlier is
## 2.92 Megabases with 40 SNPs.

## To make the combo map ped file
## Need to use the data object in the Data_explore script
## The data needed for adegenet is recoded as 12 format
## needs to be in AA CC GG TT format


# CHR 5 -------------------------------------------------------------------
lostruct_data = lostruct_run(data = tped, 
                             chr = 5, 
                             window_size = 20, 
                             k_value = 2)
outliers = Outlier_hunter(data = lostruct_data,
                          sd_percentile = 2)

outlier_full_data = Outlier_data(data = tped, 
                                 outlier_data = outliers, 
                                 chr = 5, 
                                 window_size = 20, 
                                 k_value = 2)




Chr5_map = map_maker(outlier_full_data$'23')
Chr5_ped = ped_maker(outlier_full_data$'23')


Chr5_data = Adegenet_PCA(outlier_ped = Chr5_ped, 
                         outlier_map = Chr5_map, 
                         OG_ped = OG_ped,
                         env = env_data)

Pop_that_pca(Chr5_data, 
             pop_num = 39,
             chr_num = 3, 
             win_num = 00)


# CHR 5 REGIONS -----------------------------------------------------------
## COMBINE CONTINUOUS REGIONS 15, 16, 17, & 18
Chr5_map_win15 = map_maker(outlier_full_data$'15')
Chr5_ped_win15 = ped_maker(outlier_full_data$'15')
Chr5_map_win16 = map_maker(outlier_full_data$'16')
Chr5_ped_win16 = ped_maker(outlier_full_data$'16')
Chr5_map_win17 = map_maker(outlier_full_data$'17')
Chr5_ped_win17 = ped_maker(outlier_full_data$'17')
Chr5_map_win18 = map_maker(outlier_full_data$'18')
Chr5_ped_win18 = ped_maker(outlier_full_data$'18')


Chr5_data_win15 = Adegenet_PCA(outlier_ped = Chr5_ped_win15, 
                              outlier_map = Chr5_map_win15, 
                              OG_ped = OG_ped,
                              env = env_data)

Chr5_data_win16 = Adegenet_PCA(outlier_ped = Chr5_ped_win16, 
                               outlier_map = Chr5_map_win16, 
                               OG_ped = OG_ped,
                               env = env_data)
Chr5_data_win17 = Adegenet_PCA(outlier_ped = Chr5_ped_win17, 
                               outlier_map = Chr5_map_win17, 
                               OG_ped = OG_ped,
                               env = env_data)
Chr5_data_win18 = Adegenet_PCA(outlier_ped = Chr5_ped_win18, 
                               outlier_map = Chr5_map_win18, 
                               OG_ped = OG_ped,
                               env = env_data)

Chr5_data_win16 = Chr5_data_win16 %>% 
  select(contains('AX-'))
Chr5_data_win17 = Chr5_data_win17 %>% 
  select(contains('AX-'))
Chr5_data_win18 = Chr5_data_win18 %>% 
  select(contains('AX-'))

Chr5_win15161718 = bind_cols(Chr5_data_win15, 
                        Chr5_data_win16, 
                        Chr5_data_win17, 
                        Chr5_data_win18)

Pop_that_pca(Chr5_win15161718, 
             pop_num = 39)

write_tsv(Chr5_win15161718, 
          '~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/Lostruct_Chr5_REGION1_win15161718_14.05.2021.txt')

## map file for continuous regions
Chr5_region1 = bind_rows(Chr5_map_win15, 
                         Chr5_map_win16, 
                         Chr5_map_win17, 
                         Chr5_map_win18)%>% 
  select(1:4) %>% 
  rename(`#Chromosome` = Chromosome) %>% 
  write_tsv('~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/Chr5_REGION1_14.05.2021.map')

head(Chr5_region1)
tail(Chr5_region1)

## Calculate regions size
(26593414-22808372)/1000000
## The region size of the outlier is
## 3.79 Megabases with 40 SNPs.

## To make the combo map ped file
## Need to use the data object in the Data_explore script
## The data needed for adegenet is recoded as 12 format
## needs to be in AA CC GG TT format

## CHR 5 REGION 2 win 23
Chr5_map_win23 = map_maker(outlier_full_data$'23')
Chr5_ped_win23 = ped_maker(outlier_full_data$'23')


Chr5_data_win23 = Adegenet_PCA(outlier_ped = Chr5_ped_win23, 
                               outlier_map = Chr5_map_win23, 
                               OG_ped = OG_ped,
                               env = env_data)

Pop_that_pca(Chr5_data_win23, 
             pop_num = 39)

write_tsv(Chr5_data_win23, 
          '~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/Lostruct_Chr5_REGION2_win23_14.05.2021.txt')

## map file for continuous regions
Chr5_region2 = Chr5_map_win23 %>% 
  select(1:4) %>% 
  rename(`#Chromosome` = Chromosome) %>% 
  write_tsv('~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/Chr5_REGION2_14.05.2021.map')

head(Chr5_region2)
tail(Chr5_region2)

## Calculate regions size
(32207878-31284810)/1000000
## The region size of the outlier is
## 0.92 Megabases with 40 SNPs.

## To make the combo map ped file
## Need to use the data object in the Data_explore script
## The data needed for adegenet is recoded as 12 format
## needs to be in AA CC GG TT format


# CHR6 --------------------------------------------------------------------
lostruct_data = lostruct_run(data = tped, 
                             chr = 6, 
                             window_size = 20, 
                             k_value = 2)
outliers = Outlier_hunter(data = lostruct_data,
                          sd_percentile = 2)

outlier_full_data = Outlier_data(data = tped, 
                                 outlier_data = outliers, 
                                 chr = 6, 
                                 window_size = 20, 
                                 k_value = 2)




Chr6_map = map_maker(outlier_full_data$'10')
Chr6_ped = ped_maker(outlier_full_data$'10')


Chr6_data = Adegenet_PCA(outlier_ped = Chr6_ped, 
                         outlier_map = Chr6_map, 
                         OG_ped = OG_ped,
                         env = env_data)

Pop_that_pca(Chr6_data, 
             pop_num = 39,
             chr_num = 3, 
             win_num = 00)


# CHR6 REGIONS ------------------------------------------------------------
## COMBINE CONTINUOUS REGIONS 1 & 2
Chr6_map_win1 = map_maker(outlier_full_data$'1')
Chr6_ped_win1 = ped_maker(outlier_full_data$'1')
Chr6_map_win2 = map_maker(outlier_full_data$'2')
Chr6_ped_win2 = ped_maker(outlier_full_data$'2')


Chr6_data_win1 = Adegenet_PCA(outlier_ped = Chr6_ped_win1, 
                               outlier_map = Chr6_map_win1, 
                               OG_ped = OG_ped,
                               env = env_data)

Chr6_data_win2 = Adegenet_PCA(outlier_ped = Chr6_ped_win2, 
                               outlier_map = Chr6_map_win2, 
                               OG_ped = OG_ped,
                               env = env_data)

Chr6_data_win2 = Chr6_data_win2 %>% 
  select(contains('AX-'))

Chr6_win12 = bind_cols(Chr6_data_win1, 
                             Chr6_data_win2)

Pop_that_pca(Chr6_win12, 
             pop_num = 39)

write_tsv(Chr6_win12, 
          '~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/Lostruct_Chr6_REGION1_win12_14.05.2021.txt')

## map file for continuous regions
Chr6_region1 = bind_rows(Chr6_map_win1, 
                         Chr6_map_win2)%>% 
  select(1:4) %>% 
  rename(`#Chromosome` = Chromosome) %>% 
  write_tsv('~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/Chr6_REGION1_14.05.2021.map')

head(Chr6_region1)
tail(Chr6_region1)

## REGION 1 SHOWS NO PATTERN OF AN SV


## CHR6 REGION2 WIN 9 & 10
Chr6_map_win9 = map_maker(outlier_full_data$'9')
Chr6_ped_win9 = ped_maker(outlier_full_data$'9')
Chr6_map_win10 = map_maker(outlier_full_data$'10')
Chr6_ped_win10 = ped_maker(outlier_full_data$'10')


Chr6_data_win9 = Adegenet_PCA(outlier_ped = Chr6_ped_win9, 
                              outlier_map = Chr6_map_win9, 
                              OG_ped = OG_ped,
                              env = env_data)

Chr6_data_win10 = Adegenet_PCA(outlier_ped = Chr6_ped_win10, 
                              outlier_map = Chr6_map_win10, 
                              OG_ped = OG_ped,
                              env = env_data)

Chr6_data_win10 = Chr6_data_win10 %>% 
  select(contains('AX-'))

Chr6_win910 = bind_cols(Chr6_data_win9, 
                       Chr6_data_win10)

Pop_that_pca(Chr6_win910, 
             pop_num = 39)

write_tsv(Chr6_win910, 
          '~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/Lostruct_Chr6_REGION2_win910_14.05.2021.txt')

## map file for continuous regions
Chr6_region2 = bind_rows(Chr6_map_win9, 
                         Chr6_map_win10)%>% 
  select(1:4) %>% 
  rename(`#Chromosome` = Chromosome) %>% 
  write_tsv('~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/Chr6_REGION2_14.05.2021.map')

head(Chr6_region2)
tail(Chr6_region2)

(16112100-13579128)/1000000
## The region size of the outlier is
## 2.53 Megabases with 40 SNPs.

## To make the combo map ped file
## Need to use the data object in the Data_explore script
## The data needed for adegenet is recoded as 12 format
## needs to be in AA CC GG TT format


# CHR 7 -------------------------------------------------------------------
lostruct_data = lostruct_run(data = tped, 
                             chr = 7, 
                             window_size = 20, 
                             k_value = 2)
outliers = Outlier_hunter(data = lostruct_data,
                          sd_percentile = 2)

outlier_full_data = Outlier_data(data = tped, 
                                 outlier_data = outliers, 
                                 chr = 7, 
                                 window_size = 20, 
                                 k_value = 2)




Chr7_map = map_maker(outlier_full_data$'5')
Chr7_ped = ped_maker(outlier_full_data$'5')


Chr7_data = Adegenet_PCA(outlier_ped = Chr7_ped, 
                         outlier_map = Chr7_map, 
                         OG_ped = OG_ped,
                         env = env_data)

Pop_that_pca(Chr7_data, 
             pop_num = 39,
             chr_num = 7, 
             win_num = 00)


# CHR 7 regions -----------------------------------------------------------
Chr7_map_win4 = map_maker(outlier_full_data$'4')
Chr7_ped_win4 = ped_maker(outlier_full_data$'4')
Chr7_map_win5 = map_maker(outlier_full_data$'5')
Chr7_ped_win5 = ped_maker(outlier_full_data$'5')


Chr7_data_win4 = Adegenet_PCA(outlier_ped = Chr7_ped_win4, 
                              outlier_map = Chr7_map_win4, 
                              OG_ped = OG_ped,
                              env = env_data)

Chr7_data_win5 = Adegenet_PCA(outlier_ped = Chr7_ped_win5, 
                              outlier_map = Chr7_map_win5, 
                              OG_ped = OG_ped,
                              env = env_data)

Chr7_data_win5 = Chr7_data_win5 %>% 
  select(contains('AX-'))

Chr7_win45 = bind_cols(Chr7_data_win4, 
                       Chr7_data_win5)

Pop_that_pca(Chr7_win45, 
             pop_num = 39)

write_tsv(Chr7_win45, 
          '~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/Lostruct_Chr7_REGION1_win42_14.05.2021.txt')

## map file for continuous regions
Chr7_region1 = bind_rows(Chr7_map_win4, 
                         Chr7_map_win5)%>% 
  select(1:4) %>% 
  rename(`#Chromosome` = Chromosome) %>% 
  write_tsv('~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/Chr7_REGION1_14.05.2021.map')

head(Chr7_region1)
tail(Chr7_region1)

## region size
(9888946-5783985)/1000000


# CHR 8 -------------------------------------------------------------------
lostruct_data = lostruct_run(data = tped, 
                             chr = 8, 
                             window_size = 20, 
                             k_value = 2)
outliers = Outlier_hunter(data = lostruct_data,
                          sd_percentile = 2)

outlier_full_data = Outlier_data(data = tped, 
                                 outlier_data = outliers, 
                                 chr = 8, 
                                 window_size = 20, 
                                 k_value = 2)




Chr8_map = map_maker(outlier_full_data$'6')
Chr8_ped = ped_maker(outlier_full_data$'6')


Chr8_data = Adegenet_PCA(outlier_ped = Chr8_ped, 
                         outlier_map = Chr8_map, 
                         OG_ped = OG_ped,
                         env = env_data)

Pop_that_pca(Chr8_data, 
             pop_num = 39,
             chr_num = 7, 
             win_num = 00)


ggsave('~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/Chr8_outlier_window6.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm')


# CHR 9 -------------------------------------------------------------------
lostruct_data = lostruct_run(data = tped, 
                             chr = 9, 
                             window_size = 20, 
                             k_value = 2)
outliers = Outlier_hunter(data = lostruct_data,
                          sd_percentile = 2)

outlier_full_data = Outlier_data(data = tped, 
                                 outlier_data = outliers, 
                                 chr = 9, 
                                 window_size = 20, 
                                 k_value = 2)




Chr9_map = map_maker(outlier_full_data$'4')
Chr9_ped = ped_maker(outlier_full_data$'4')


Chr9_data = Adegenet_PCA(outlier_ped = Chr9_ped, 
                         outlier_map = Chr9_map, 
                         OG_ped = OG_ped,
                         env = env_data)

Pop_that_pca(Chr9_data, 
             pop_num = 39,
             chr_num = 7, 
             win_num = 00)


ggsave('~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/Chr9_outlier_window12.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm')


# CHR 9 REGION ------------------------------------------------------------
Chr9_map_win4 = map_maker(outlier_full_data$'4')
Chr9_ped_win4 = ped_maker(outlier_full_data$'4')


Chr9_data_win4 = Adegenet_PCA(outlier_ped = Chr9_ped_win4, 
                              outlier_map = Chr9_map_win4, 
                              OG_ped = OG_ped,
                              env = env_data)

Pop_that_pca(Chr9_data_win4, 
             pop_num = 39)

write_tsv(Chr9_data_win4, 
          '~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/Lostruct_Chr9_REGION1_win4_14.05.2021.txt')

## map file for continuous regions
Chr9_region1 = Chr9_map_win4%>% 
  select(1:4) %>% 
  rename(`#Chromosome` = Chromosome) %>% 
  write_tsv('~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/Chr9_REGION1_14.05.2021.map')

head(Chr9_region1)
tail(Chr9_region1)

## region size
(3499661-2633420)/1000000

##
# Chr 11 ------------------------------------------------------------------

lostruct_data = lostruct_run(data = tped, 
                             chr = 11, 
                             window_size = 20, 
                             k_value = 2)

# MDS_survey(lostruct_data)

outliers = Outlier_hunter(data = lostruct_data,
                          sd_percentile = 2)

# Outlier_plots(normal_data = lostruct_data,
#               outlier_data = outliers)

outlier_full_data = Outlier_data(data = tped, 
                                 outlier_data = outliers, 
                                 chr = 11, 
                                 window_size = 20, 
                                 k_value = 2)


outlier_full_data


# chr 11 Combine consequtive outlier windows -------------------------------------

## combining windows 7 and 8
chr_map_win7 = map_maker(outlier_full_data$'7')
chr_ped_win7 = ped_maker(outlier_full_data$'7')
chr_data_win7 = Adegenet_PCA(outlier_ped = chr_ped_win7, 
                             outlier_map = chr_map_win7, 
                             OG_ped = OG_ped,
                             env = env_data)

chr_map_win8 = map_maker(outlier_full_data$'8')
chr_ped_win8 = ped_maker(outlier_full_data$'8')
chr_data_win8 = Adegenet_PCA(outlier_ped = chr_ped_win8, 
                             outlier_map = chr_map_win8, 
                             OG_ped = OG_ped,
                             env = env_data)



chr_data_win8 = chr_data_win8 %>% 
  select(contains('AX-'))

chr_combo_win78 = bind_cols(chr_data_win7, 
                            chr_data_win8)

PCA_outlier_wins7and8 = Pop_that_pca(chr_combo_win78, 
                                     pop_num = 38,
                                     chr_num = 11, 
                                     win_num = 78)

## IMPORTANT NOTE
## When GDL in not inluded window 15 comes up as an outlier
## when GDL is included in the analysis GDL doesn't come up
## as an outlier


# Make combo map and ped files for continuous regions ---------------------

## Making the map files for continouous regions
Chr_map_win7 = map_maker(outlier_full_data$'7')
Chr_map_win8 = map_maker(outlier_full_data$'8')

chr11_region1_outlier_map = bind_rows(Chr_map_win7, 
                                      Chr_map_win8) %>% 
  select(1:4) %>% 
  rename(`#Chromosome` = Chromosome, 
         `Marker ID` = MarkerID, 
         `Genetic distance` = Genetic_dist, 
         `Physical distance` = Physical_dist)

write_tsv(chr11_region1_outlier_map, 
          '~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/Chr11_region1_outlier_windows.map')

## Calculate regions size
region_size = 14801093-11095966
region_size/1000000
## The region size of the outlier is
## 3.7 Megabases with 40 SNPs. 
## SNP density is DEFINITELY and issue in identifying 
## these structural variants. 

## Making the ped files for continuous regions
chr_ped_win7 = ped_maker(outlier_full_data$'7')
chr_ped_win8 = ped_maker(outlier_full_data$'8')

combo_ped = bind_cols(chr_ped_win7, 
                      chr_ped_win8)
combo_ped = OG_ped %>% 
  select(1:6) %>% 
  bind_cols(combo_ped) %>%
  rename(`#FamilyID` = 1) 


OG_ped %>% 
  select(1:6, 
         'AX-181987266':'AX-181928438') %>% 
  rename(`#FamilyID` = 1) %>% 
  View()
write_tsv('~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/Chr11_region1_outlier_windows.ped')
# test = names(combo_ped) %>% 
#   as_tibble() %>% 
#   slice(7:46) %>% 
#   rename(MarkerID = value)

# %>% 
#   write_tsv('~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/Chr11_region1_outlier_windows.ped')
