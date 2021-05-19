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


# Chr 10 ------------------------------------------------------------------
lostruct_data = lostruct_run(data = tped, 
                             chr = 10, 
                             window_size = 20, 
                             k_value = 2)
outliers = Outlier_hunter(data = lostruct_data,
                          sd_percentile = 2)

outlier_full_data = Outlier_data(data = tped, 
                                 outlier_data = outliers, 
                                 chr = 10, 
                                 window_size = 20, 
                                 k_value = 2)




Chr10_map = map_maker(outlier_full_data$'5')
Chr10_ped = ped_maker(outlier_full_data$'5')


Chr10_data = Adegenet_PCA(outlier_ped = Chr10_ped, 
                         outlier_map = Chr10_map, 
                         OG_ped = OG_ped,
                         env = env_data)

Pop_that_pca(Chr10_data, 
             pop_num = 39,
             chr_num = 7, 
             win_num = 00)


# ggsave('~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/Chr10_outlier_window12.tiff', 
#        plot = last_plot(), 
#        dpi = 'retina', 
#        units = 'cm')


# Chr 10 REGION -----------------------------------------------------------

Chr10_map_win4 = map_maker(outlier_full_data$'4')
Chr10_ped_win4 = ped_maker(outlier_full_data$'4')
Chr10_map_win5 = map_maker(outlier_full_data$'5')
Chr10_ped_win5 = ped_maker(outlier_full_data$'5')


Chr10_data_win4 = Adegenet_PCA(outlier_ped = Chr10_ped_win4, 
                              outlier_map = Chr10_map_win4, 
                              OG_ped = OG_ped,
                              env = env_data)

Chr10_data_win5 = Adegenet_PCA(outlier_ped = Chr10_ped_win5, 
                              outlier_map = Chr10_map_win5, 
                              OG_ped = OG_ped,
                              env = env_data)

Chr10_data_win5 = Chr10_data_win5 %>% 
  select(contains('AX-'))

Chr10_win45 = bind_cols(Chr10_data_win4, 
                       Chr10_data_win5)

Pop_that_pca(Chr10_win45, 
             pop_num = 39)

write_tsv(Chr10_win45, 
          '~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/Lostruct_Chr10_REGION1_win45_14.05.2021.txt')

## map file for continuous regions
Chr10_region1 = bind_rows(Chr10_map_win4, 
                         Chr10_map_win5)%>% 
  select(1:4) %>% 
  rename(`#Chromosome` = Chromosome) %>% 
  write_tsv('~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/CHR10_REGION1_14.05.2021.map')

head(Chr10_region1)
tail(Chr10_region1)

## region size
(7773641-4844048)/1000000


# CHR 11 ------------------------------------------------------------------

lostruct_data = lostruct_run(data = tped, 
                             chr = 11, 
                             window_size = 20, 
                             k_value = 2)
outliers = Outlier_hunter(data = lostruct_data,
                          sd_percentile = 2)

outlier_full_data = Outlier_data(data = tped, 
                                 outlier_data = outliers, 
                                 chr = 11, 
                                 window_size = 20, 
                                 k_value = 2)

Chr11_map = map_maker(outlier_full_data$'5')
Chr11_ped = ped_maker(outlier_full_data$'5')


Chr11_data = Adegenet_PCA(outlier_ped = Chr11_ped, 
                          outlier_map = Chr11_map, 
                          OG_ped = OG_ped,
                          env = env_data)

Pop_that_pca(Chr11_data, 
             pop_num = 39,
             chr_num = 7, 
             win_num = 00)

# CHR 11 REGIONS ----------------------------------------------------------
## CHR11 region1 COmbining windows 2 and 3
Chr11_map_win2 = map_maker(outlier_full_data$'2')
Chr11_ped_win2 = ped_maker(outlier_full_data$'2')
Chr11_map_win3 = map_maker(outlier_full_data$'3')
Chr11_ped_win3 = ped_maker(outlier_full_data$'3')


Chr11_data_win2 = Adegenet_PCA(outlier_ped = Chr11_ped_win2, 
                               outlier_map = Chr11_map_win2, 
                               OG_ped = OG_ped,
                               env = env_data)

Chr11_data_win3 = Adegenet_PCA(outlier_ped = Chr11_ped_win3, 
                               outlier_map = Chr11_map_win3, 
                               OG_ped = OG_ped,
                               env = env_data)

Chr11_data_win3 = Chr11_data_win3 %>% 
  select(contains('AX-'))

Chr11_win23 = bind_cols(Chr11_data_win2, 
                        Chr11_data_win3)

Pop_that_pca(Chr11_win23, 
             pop_num = 39)

write_tsv(Chr11_win23, 
          '~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/Lostruct_Chr11_REGION1_win23_14.05.2021.txt')

## map file for continuous regions
Chr11_region1 = bind_rows(Chr11_map_win2, 
                          Chr11_map_win3)%>% 
  select(1:4) %>% 
  rename(`#Chromosome` = Chromosome) %>% 
  write_tsv('~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/Chr11_REGION1_14.05.2021.map')

head(Chr11_region1)
tail(Chr11_region1)

## region size
(4994889-2611798)/1000000
## 2.38 Mb


## CHR 11 REGION 2 window 5
Chr11_map_win5 = map_maker(outlier_full_data$'5')
Chr11_ped_win5 = ped_maker(outlier_full_data$'5')

Chr11_data_win5 = Adegenet_PCA(outlier_ped = Chr11_ped_win5, 
                               outlier_map = Chr11_map_win5, 
                               OG_ped = OG_ped,
                               env = env_data)

Pop_that_pca(Chr11_data_win5, 
             pop_num = 39)

write_tsv(Chr11_data_win5, 
          '~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/Lostruct_Chr11_REGION2_win5_14.05.2021.txt')

## map file for continuous regions
Chr11_region2 = Chr11_map_win5 %>% 
  select(1:4) %>% 
  rename(`#Chromosome` = Chromosome) %>% 
  write_tsv('~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/Chr11_REGION2_14.05.2021.map')

head(Chr11_region2)
tail(Chr11_region2)

## region size
(10393035-7226516)/1000000
## 3.17Mb


## combing all three regions. THey overlap based on the 
## pcadmix results


Chr11_map_win2 = map_maker(outlier_full_data$'2')
Chr11_ped_win2 = ped_maker(outlier_full_data$'2')
Chr11_map_win3 = map_maker(outlier_full_data$'3')
Chr11_ped_win3 = ped_maker(outlier_full_data$'3')
Chr11_map_win5 = map_maker(outlier_full_data$'5')
Chr11_ped_win5 = ped_maker(outlier_full_data$'5')


Chr11_data_win2 = Adegenet_PCA(outlier_ped = Chr11_ped_win2, 
                               outlier_map = Chr11_map_win2, 
                               OG_ped = OG_ped,
                               env = env_data)

Chr11_data_win3 = Adegenet_PCA(outlier_ped = Chr11_ped_win3, 
                               outlier_map = Chr11_map_win3, 
                               OG_ped = OG_ped,
                               env = env_data)
Chr11_data_win5 = Adegenet_PCA(outlier_ped = Chr11_ped_win5, 
                               outlier_map = Chr11_map_win5, 
                               OG_ped = OG_ped,
                               env = env_data)

Chr11_data_win3 = Chr11_data_win3 %>% 
  select(contains('AX-'))
Chr11_data_win5 = Chr11_data_win5 %>% 
  select(contains('AX-'))

Chr11_win235 = bind_cols(Chr11_data_win2, 
                        Chr11_data_win3, 
                        Chr11_data_win5)

Pop_that_pca(Chr11_win235, 
             pop_num = 39)

write_tsv(Chr11_win235, 
          '~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/Lostruct_Chr11_REGION3_win235_14.05.2021.txt')

## map file for continuous regions
Chr11_region3 = bind_rows(Chr11_map_win2, 
                          Chr11_map_win3, 
                          Chr11_map_win5)%>% 
  select(1:4) %>% 
  rename(`#Chromosome` = Chromosome) %>% 
  write_tsv('~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/Chr11_REGION3_14.05.2021.map')

head(Chr11_region3)
tail(Chr11_region3)

## region size
(10393035-2611798)/1000000
## 7.78 Mb

# CHR 12 ------------------------------------------------------------------

lostruct_data = lostruct_run(data = tped, 
                             chr = 12, 
                             window_size = 20, 
                             k_value = 2)
outliers = Outlier_hunter(data = lostruct_data,
                          sd_percentile = 2)

outlier_full_data = Outlier_data(data = tped, 
                                 outlier_data = outliers, 
                                 chr = 12, 
                                 window_size = 20, 
                                 k_value = 2)

Chr12_map = map_maker(outlier_full_data$'16')
Chr12_ped = ped_maker(outlier_full_data$'16')


Chr12_data = Adegenet_PCA(outlier_ped = Chr12_ped, 
                          outlier_map = Chr12_map, 
                          OG_ped = OG_ped,
                          env = env_data)

Pop_that_pca(Chr12_data, 
             pop_num = 39,
             chr_num = 7, 
             win_num = 00)

# CHR 12 REGIONS ----------------------------------------------------------

Chr12_map_win8 = map_maker(outlier_full_data$'8')
Chr12_ped_win8 = ped_maker(outlier_full_data$'8')
Chr12_map_win9 = map_maker(outlier_full_data$'9')
Chr12_ped_win9 = ped_maker(outlier_full_data$'9')


Chr12_data_win8 = Adegenet_PCA(outlier_ped = Chr12_ped_win8, 
                               outlier_map = Chr12_map_win8, 
                               OG_ped = OG_ped,
                               env = env_data)

Chr12_data_win9 = Adegenet_PCA(outlier_ped = Chr12_ped_win9, 
                               outlier_map = Chr12_map_win9, 
                               OG_ped = OG_ped,
                               env = env_data)

Chr12_data_win9 = Chr12_data_win9 %>% 
  select(contains('AX-'))

Chr12_win89 = bind_cols(Chr12_data_win8, 
                        Chr12_data_win9)

Pop_that_pca(Chr12_win89, 
             pop_num = 39)

ggsave('~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/CHR12_REGION1_win89_PCA.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm')

## Doesn't look like  a structural variant

## Combined windows 14, 15, 16
Chr12_map_win14 = map_maker(outlier_full_data$'14')
Chr12_ped_win14 = ped_maker(outlier_full_data$'14')
Chr12_map_win15 = map_maker(outlier_full_data$'15')
Chr12_ped_win15 = ped_maker(outlier_full_data$'15')
Chr12_map_win16 = map_maker(outlier_full_data$'16')
Chr12_ped_win16 = ped_maker(outlier_full_data$'16')


Chr12_data_win14 = Adegenet_PCA(outlier_ped = Chr12_ped_win14, 
                               outlier_map = Chr12_map_win14, 
                               OG_ped = OG_ped,
                               env = env_data)

Chr12_data_win15 = Adegenet_PCA(outlier_ped = Chr12_ped_win15, 
                               outlier_map = Chr12_map_win15, 
                               OG_ped = OG_ped,
                               env = env_data)
Chr12_data_win16 = Adegenet_PCA(outlier_ped = Chr12_ped_win16, 
                                outlier_map = Chr12_map_win16, 
                                OG_ped = OG_ped,
                                env = env_data)

Chr12_data_win15 = Chr12_data_win15 %>% 
  select(contains('AX-'))
Chr12_data_win16 = Chr12_data_win16 %>% 
  select(contains('AX-'))

Chr12_win141516 = bind_cols(Chr12_data_win14, 
                        Chr12_data_win15, 
                        Chr12_data_win16)

Pop_that_pca(Chr12_win141516, 
             pop_num = 39)

## Doesn't look like a structural variant. 

ggsave('~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/CHR12_REGION2_win141516_PCA.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm')


# CHR 13 ------------------------------------------------------------------


lostruct_data = lostruct_run(data = tped, 
                             chr = 13, 
                             window_size = 20, 
                             k_value = 2)
outliers = Outlier_hunter(data = lostruct_data,
                          sd_percentile = 2)

outlier_full_data = Outlier_data(data = tped, 
                                 outlier_data = outliers, 
                                 chr = 13, 
                                 window_size = 20, 
                                 k_value = 2)

Chr13_map = map_maker(outlier_full_data$'1')
Chr13_ped = ped_maker(outlier_full_data$'1')


Chr13_data = Adegenet_PCA(outlier_ped = Chr13_ped, 
                          outlier_map = Chr13_map, 
                          OG_ped = OG_ped,
                          env = env_data)

Pop_that_pca(Chr13_data, 
             pop_num = 39,
             chr_num = 7, 
             win_num = 00)

## The outliers at window 1 and 6 don't look like structural 
## variants

ggsave('~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/CHR13_REGION1_win1_PCA.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm')

# Chr 14 ------------------------------------------------------------------
lostruct_data = lostruct_run(data = tped, 
                             chr = 14, 
                             window_size = 20, 
                             k_value = 2)
outliers = Outlier_hunter(data = lostruct_data,
                          sd_percentile = 2)

outlier_full_data = Outlier_data(data = tped, 
                                 outlier_data = outliers, 
                                 chr = 14, 
                                 window_size = 20, 
                                 k_value = 2)

Chr14_map = map_maker(outlier_full_data$'11')
Chr14_ped = ped_maker(outlier_full_data$'11')


Chr14_data = Adegenet_PCA(outlier_ped = Chr14_ped, 
                          outlier_map = Chr14_map, 
                          OG_ped = OG_ped,
                          env = env_data)

Pop_that_pca(Chr14_data, 
             pop_num = 39,
             chr_num = 7, 
             win_num = 00)


# CHR 14 regions ----------------------------------------------------------

Chr14_map_win7 = map_maker(outlier_full_data$'7')
Chr14_ped_win7 = ped_maker(outlier_full_data$'7')
Chr14_map_win8 = map_maker(outlier_full_data$'8')
Chr14_ped_win8 = ped_maker(outlier_full_data$'8')


Chr14_data_win7 = Adegenet_PCA(outlier_ped = Chr14_ped_win7, 
                               outlier_map = Chr14_map_win7, 
                               OG_ped = OG_ped,
                               env = env_data)

Chr14_data_win8 = Adegenet_PCA(outlier_ped = Chr14_ped_win8, 
                               outlier_map = Chr14_map_win8, 
                               OG_ped = OG_ped,
                               env = env_data)

Chr14_data_win8 = Chr14_data_win8 %>% 
  select(contains('AX-'))

Chr14_win78 = bind_cols(Chr14_data_win7, 
                         Chr14_data_win8)

Pop_that_pca(Chr14_win78, 
             pop_num = 39)

## Definitely not a structural variant
ggsave('~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/CHR14_windows78_PCA.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm')

## combining windows 10 and 11 to identify structural variants
Chr14_map_win10 = map_maker(outlier_full_data$'10')
Chr14_ped_win10 = ped_maker(outlier_full_data$'10')
Chr14_map_win11 = map_maker(outlier_full_data$'11')
Chr14_ped_win11 = ped_maker(outlier_full_data$'11')


Chr14_data_win10 = Adegenet_PCA(outlier_ped = Chr14_ped_win10, 
                               outlier_map = Chr14_map_win10, 
                               OG_ped = OG_ped,
                               env = env_data)

Chr14_data_win11 = Adegenet_PCA(outlier_ped = Chr14_ped_win11, 
                               outlier_map = Chr14_map_win11, 
                               OG_ped = OG_ped,
                               env = env_data)

Chr14_data_win11 = Chr14_data_win11 %>% 
  select(contains('AX-'))

Chr14_win1011 = bind_cols(Chr14_data_win10, 
                        Chr14_data_win11)

Pop_that_pca(Chr14_win1011, 
             pop_num = 39)
## doesn't really look like a stuctural variant
ggsave('~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/CHR14_windows1011_PCA.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm')

# CHR 15 ------------------------------------------------------------------
lostruct_data = lostruct_run(data = tped, 
                             chr = 15, 
                             window_size = 20, 
                             k_value = 2)
outliers = Outlier_hunter(data = lostruct_data,
                          sd_percentile = 2)

outlier_full_data = Outlier_data(data = tped, 
                                 outlier_data = outliers, 
                                 chr = 15, 
                                 window_size = 20, 
                                 k_value = 2)

## No outliers

# CHR16 -------------------------------------------------------------------

lostruct_data = lostruct_run(data = tped, 
                             chr = 16, 
                             window_size = 20, 
                             k_value = 2)
outliers = Outlier_hunter(data = lostruct_data,
                          sd_percentile = 2)

outlier_full_data = Outlier_data(data = tped, 
                                 outlier_data = outliers, 
                                 chr = 16, 
                                 window_size = 20, 
                                 k_value = 2)
Chr16_map = map_maker(outlier_full_data$'23')
Chr16_ped = ped_maker(outlier_full_data$'23')


Chr16_data = Adegenet_PCA(outlier_ped = Chr16_ped, 
                          outlier_map = Chr16_map, 
                          OG_ped = OG_ped,
                          env = env_data)

Pop_that_pca(Chr16_data, 
             pop_num = 39,
             chr_num = 7, 
             win_num = 00)

# CHR 16 regions------------------------------------------------------------------
Chr16_map_win18 = map_maker(outlier_full_data$'18')
Chr16_ped_win18 = ped_maker(outlier_full_data$'18')
Chr16_map_win19 = map_maker(outlier_full_data$'19')
Chr16_ped_win19 = ped_maker(outlier_full_data$'19')


Chr16_data_win18 = Adegenet_PCA(outlier_ped = Chr16_ped_win18, 
                               outlier_map = Chr16_map_win18, 
                               OG_ped = OG_ped,
                               env = env_data)

Chr16_data_win19 = Adegenet_PCA(outlier_ped = Chr16_ped_win19, 
                               outlier_map = Chr16_map_win19, 
                               OG_ped = OG_ped,
                               env = env_data)

Chr16_data_win19 = Chr16_data_win19 %>% 
  select(contains('AX-'))

Chr16_win1819 = bind_cols(Chr16_data_win18, 
                         Chr16_data_win19)

Pop_that_pca(Chr16_win1819, 
             pop_num = 39)

## there's no structural variant here

## combine regions 22 & 23
Chr16_map_win22 = map_maker(outlier_full_data$'22')
Chr16_ped_win22 = ped_maker(outlier_full_data$'22')
Chr16_map_win23 = map_maker(outlier_full_data$'23')
Chr16_ped_win23 = ped_maker(outlier_full_data$'23')


Chr16_data_win22 = Adegenet_PCA(outlier_ped = Chr16_ped_win22, 
                                outlier_map = Chr16_map_win22, 
                                OG_ped = OG_ped,
                                env = env_data)

Chr16_data_win23 = Adegenet_PCA(outlier_ped = Chr16_ped_win23, 
                                outlier_map = Chr16_map_win23, 
                                OG_ped = OG_ped,
                                env = env_data)

Chr16_data_win23 = Chr16_data_win23 %>% 
  select(contains('AX-'))

Chr16_win2223 = bind_cols(Chr16_data_win22, 
                          Chr16_data_win23)

Pop_that_pca(Chr16_win2223, 
             pop_num = 39)

## Doesn't really look like a structural variant


# CHR 17 ------------------------------------------------------------------
lostruct_data = lostruct_run(data = tped, 
                             chr = 17, 
                             window_size = 20, 
                             k_value = 2)
outliers = Outlier_hunter(data = lostruct_data,
                          sd_percentile = 2)

outlier_full_data = Outlier_data(data = tped, 
                                 outlier_data = outliers, 
                                 chr = 17, 
                                 window_size = 20, 
                                 k_value = 2)

Chr17_map = map_maker(outlier_full_data$'11')
Chr17_ped = ped_maker(outlier_full_data$'11')


Chr17_data = Adegenet_PCA(outlier_ped = Chr17_ped, 
                          outlier_map = Chr17_map, 
                          OG_ped = OG_ped,
                          env = env_data)

Pop_that_pca(Chr17_data, 
             pop_num = 39,
             chr_num = 7, 
             win_num = 00)

# CHR 17 Regions ----------------------------------------------------------

Chr17_map_win2 = map_maker(outlier_full_data$'2')
Chr17_ped_win2 = ped_maker(outlier_full_data$'2')
Chr17_map_win3 = map_maker(outlier_full_data$'3')
Chr17_ped_win3 = ped_maker(outlier_full_data$'3')
Chr17_map_win4 = map_maker(outlier_full_data$'4')
Chr17_ped_win4 = ped_maker(outlier_full_data$'4')
Chr17_map_win5 = map_maker(outlier_full_data$'5')
Chr17_ped_win5 = ped_maker(outlier_full_data$'5')


Chr17_data_win2 = Adegenet_PCA(outlier_ped = Chr17_ped_win2, 
                               outlier_map = Chr17_map_win2, 
                               OG_ped = OG_ped,
                               env = env_data)

Chr17_data_win3 = Adegenet_PCA(outlier_ped = Chr17_ped_win3, 
                               outlier_map = Chr17_map_win3, 
                               OG_ped = OG_ped,
                               env = env_data)
Chr17_data_win4 = Adegenet_PCA(outlier_ped = Chr17_ped_win4, 
                               outlier_map = Chr17_map_win4, 
                               OG_ped = OG_ped,
                               env = env_data)

Chr17_data_win5 = Adegenet_PCA(outlier_ped = Chr17_ped_win5, 
                               outlier_map = Chr17_map_win5, 
                               OG_ped = OG_ped,
                               env = env_data)

Chr17_data_win3 = Chr17_data_win3 %>% 
  select(contains('AX-'))
Chr17_data_win4 = Chr17_data_win4 %>% 
  select(contains('AX-'))
Chr17_data_win5 = Chr17_data_win5 %>% 
  select(contains('AX-'))

Chr17_win2345 = bind_cols(Chr17_data_win2, 
                         Chr17_data_win3, 
                         Chr17_data_win4, 
                         Chr17_data_win5)

Pop_that_pca(Chr17_win2345, 
             pop_num = 39)

write_tsv(Chr17_win2345, 
          '~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/Lostruct_Chr17_REGION1_win2345_14.05.2021.txt')

## map file for continuous regions
Chr17_region3 = bind_rows(Chr17_map_win2, 
                          Chr17_map_win3,
                          Chr17_map_win4,
                          Chr17_map_win5)%>% 
  select(1:4) %>% 
  rename(`#Chromosome` = Chromosome) %>% 
  write_tsv('~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/Chr17_REGION1_14.05.2021.map')

head(Chr17_region3)
tail(Chr17_region3)

## region size
(5836567-1864217)/1000000
## 3.97 Mb

## REGIONS 8-11 were flagged as outliers, 
## just want to see if it's actualy an SV
## Test out regions 8-11
Chr17_map_win8 = map_maker(outlier_full_data$'8')
Chr17_ped_win8 = ped_maker(outlier_full_data$'8')
Chr17_map_win9 = map_maker(outlier_full_data$'9')
Chr17_ped_win9 = ped_maker(outlier_full_data$'9')
Chr17_map_win10 = map_maker(outlier_full_data$'10')
Chr17_ped_win10 = ped_maker(outlier_full_data$'10')
Chr17_map_win11 = map_maker(outlier_full_data$'11')
Chr17_ped_win11 = ped_maker(outlier_full_data$'11')


Chr17_data_win8 = Adegenet_PCA(outlier_ped = Chr17_ped_win8, 
                               outlier_map = Chr17_map_win8, 
                               OG_ped = OG_ped,
                               env = env_data)

Chr17_data_win9 = Adegenet_PCA(outlier_ped = Chr17_ped_win9, 
                               outlier_map = Chr17_map_win9, 
                               OG_ped = OG_ped,
                               env = env_data)
Chr17_data_win10 = Adegenet_PCA(outlier_ped = Chr17_ped_win10, 
                               outlier_map = Chr17_map_win10, 
                               OG_ped = OG_ped,
                               env = env_data)

Chr17_data_win11 = Adegenet_PCA(outlier_ped = Chr17_ped_win11, 
                               outlier_map = Chr17_map_win11, 
                               OG_ped = OG_ped,
                               env = env_data)

Chr17_data_win9 = Chr17_data_win9 %>% 
  select(contains('AX-'))
Chr17_data_win10 = Chr17_data_win10 %>% 
  select(contains('AX-'))
Chr17_data_win11 = Chr17_data_win11 %>% 
  select(contains('AX-'))

Chr17_win891011 = bind_cols(Chr17_data_win8, 
                          Chr17_data_win9, 
                          Chr17_data_win10, 
                          Chr17_data_win11)

Pop_that_pca(Chr17_win891011, 
             pop_num = 39)

## doesn't really look like an SV
ggsave('~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/CHR17_win891011_PCA.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm')


# CHR18 -------------------------------------------------------------------

lostruct_data = lostruct_run(data = tped, 
                             chr = 18, 
                             window_size = 20, 
                             k_value = 2)
outliers = Outlier_hunter(data = lostruct_data,
                          sd_percentile = 2)

outlier_full_data = Outlier_data(data = tped, 
                                 outlier_data = outliers, 
                                 chr = 18, 
                                 window_size = 20, 
                                 k_value = 2)

Chr18_map = map_maker(outlier_full_data$'18')
Chr18_ped = ped_maker(outlier_full_data$'18')


Chr18_data = Adegenet_PCA(outlier_ped = Chr18_ped, 
                          outlier_map = Chr18_map, 
                          OG_ped = OG_ped,
                          env = env_data)

Pop_that_pca(Chr18_data, 
             pop_num = 39,
             chr_num = 7, 
             win_num = 00)

# CHR 18 REGIONS ----------------------------------------------------------
Chr18_map_win6 = map_maker(outlier_full_data$'6')
Chr18_ped_win6 = ped_maker(outlier_full_data$'6')
Chr18_map_win7 = map_maker(outlier_full_data$'7')
Chr18_ped_win7 = ped_maker(outlier_full_data$'7')
Chr18_map_win8 = map_maker(outlier_full_data$'8')
Chr18_ped_win8 = ped_maker(outlier_full_data$'8')


Chr18_data_win6 = Adegenet_PCA(outlier_ped = Chr18_ped_win6, 
                               outlier_map = Chr18_map_win6, 
                               OG_ped = OG_ped,
                               env = env_data)

Chr18_data_win7 = Adegenet_PCA(outlier_ped = Chr18_ped_win7, 
                               outlier_map = Chr18_map_win7, 
                               OG_ped = OG_ped,
                               env = env_data)
Chr18_data_win8 = Adegenet_PCA(outlier_ped = Chr18_ped_win8, 
                               outlier_map = Chr18_map_win8, 
                               OG_ped = OG_ped,
                               env = env_data)


Chr18_data_win7 = Chr18_data_win7 %>% 
  select(contains('AX-'))
Chr18_data_win8 = Chr18_data_win8 %>% 
  select(contains('AX-'))

Chr18_win678 = bind_cols(Chr18_data_win6,
                         Chr18_data_win7, 
                          Chr18_data_win8)

Pop_that_pca(Chr18_win678, 
             pop_num = 39)

## yep that's fucked. Theres structure but it's weird
ggsave('~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/CHR18_win678_PCA.tiff', 
       plot = last_plot(), 
       dpi = 'retina', 
       units = 'cm')


## Need to combine 16, 17, 18
Chr18_map_win16 = map_maker(outlier_full_data$'16')
Chr18_ped_win16 = ped_maker(outlier_full_data$'16')
Chr18_map_win17 = map_maker(outlier_full_data$'17')
Chr18_ped_win17 = ped_maker(outlier_full_data$'17')
Chr18_map_win18 = map_maker(outlier_full_data$'18')
Chr18_ped_win18 = ped_maker(outlier_full_data$'18')


Chr18_data_win16 = Adegenet_PCA(outlier_ped = Chr18_ped_win16, 
                               outlier_map = Chr18_map_win16, 
                               OG_ped = OG_ped,
                               env = env_data)

Chr18_data_win17 = Adegenet_PCA(outlier_ped = Chr18_ped_win17, 
                               outlier_map = Chr18_map_win17, 
                               OG_ped = OG_ped,
                               env = env_data)
Chr18_data_win18 = Adegenet_PCA(outlier_ped = Chr18_ped_win18, 
                               outlier_map = Chr18_map_win18, 
                               OG_ped = OG_ped,
                               env = env_data)


Chr18_data_win17 = Chr18_data_win17 %>% 
  select(contains('AX-'))
Chr18_data_win18 = Chr18_data_win18 %>% 
  select(contains('AX-'))

Chr18_win161718 = bind_cols(Chr18_data_win16,
                         Chr18_data_win17, 
                         Chr18_data_win18)

Pop_that_pca(Chr18_win161718, 
             pop_num = 39)

write_tsv(Chr18_win161718, 
          '~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/Lostruct_Chr18_REGION2_win161718_14.05.2021.txt')

## map file for continuous regions
Chr18_region2 = bind_rows(Chr18_map_win16, 
                          Chr18_map_win17,
                          Chr18_map_win18)%>% 
  select(1:4) %>% 
  rename(`#Chromosome` = Chromosome) %>% 
  write_tsv('~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/Chr18_REGION2_14.05.2021.map')

head(Chr18_region2)
tail(Chr18_region2)

## region size
(33322062-27640399)/1000000
## 5.68 Mb


# CHR19 -------------------------------------------------------------------

lostruct_data = lostruct_run(data = tped, 
                             chr = 19, 
                             window_size = 20, 
                             k_value = 2)
outliers = Outlier_hunter(data = lostruct_data,
                          sd_percentile = 2)

outlier_full_data = Outlier_data(data = tped, 
                                 outlier_data = outliers, 
                                 chr = 19, 
                                 window_size = 20, 
                                 k_value = 2)

Chr19_map = map_maker(outlier_full_data$'11')
Chr19_ped = ped_maker(outlier_full_data$'11')


Chr19_data = Adegenet_PCA(outlier_ped = Chr19_ped, 
                          outlier_map = Chr19_map, 
                          OG_ped = OG_ped,
                          env = env_data)

Pop_that_pca(Chr19_data, 
             pop_num = 39,
             chr_num = 7, 
             win_num = 00)


# CHR 19 region -----------------------------------------------------------

Chr19_map_win10 = map_maker(outlier_full_data$'10')
Chr19_ped_win10 = ped_maker(outlier_full_data$'10')
Chr19_map_win11 = map_maker(outlier_full_data$'11')
Chr19_ped_win11 = ped_maker(outlier_full_data$'11')


Chr19_data_win10 = Adegenet_PCA(outlier_ped = Chr19_ped_win10, 
                                outlier_map = Chr19_map_win10, 
                                OG_ped = OG_ped,
                                env = env_data)

Chr19_data_win11 = Adegenet_PCA(outlier_ped = Chr19_ped_win11, 
                                outlier_map = Chr19_map_win11, 
                                OG_ped = OG_ped,
                                env = env_data)


Chr19_data_win11 = Chr19_data_win11 %>% 
  select(contains('AX-'))

Chr19_win1011 = bind_cols(Chr19_data_win10,
                            Chr19_data_win11)

Pop_that_pca(Chr19_win1011, 
             pop_num = 39)

write_tsv(Chr19_win1011, 
          '~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/Lostruct_Chr19_REGION1_win1011_14.05.2021.txt')

## map file for continuous regions
Chr19_region1 = bind_rows(Chr19_map_win10, 
                          Chr19_map_win11)%>% 
  select(1:4) %>% 
  rename(`#Chromosome` = Chromosome) %>% 
  write_tsv('~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/Chr19_REGION1_14.05.2021.map')

head(Chr19_region1)
tail(Chr19_region1)

## region size
(20969278-17697514)/1000000


# CHR20 -------------------------------------------------------------------

lostruct_data = lostruct_run(data = tped, 
                             chr = 20, 
                             window_size = 20, 
                             k_value = 2)
outliers = Outlier_hunter(data = lostruct_data,
                          sd_percentile = 2)

outlier_full_data = Outlier_data(data = tped, 
                                 outlier_data = outliers, 
                                 chr = 20, 
                                 window_size = 20, 
                                 k_value = 2)

Chr20_map = map_maker(outlier_full_data$'13')
Chr20_ped = ped_maker(outlier_full_data$'13')


Chr20_data = Adegenet_PCA(outlier_ped = Chr20_ped, 
                          outlier_map = Chr20_map, 
                          OG_ped = OG_ped,
                          env = env_data)

Pop_that_pca(Chr20_data, 
             pop_num = 39,
             chr_num = 7, 
             win_num = 00)

# CHR20 REGIONS -----------------------------------------------------------

Chr20_map_win10 = map_maker(outlier_full_data$'10')
Chr20_ped_win10 = ped_maker(outlier_full_data$'10')
Chr20_map_win11 = map_maker(outlier_full_data$'11')
Chr20_ped_win11 = ped_maker(outlier_full_data$'11')
Chr20_map_win12 = map_maker(outlier_full_data$'12')
Chr20_ped_win12 = ped_maker(outlier_full_data$'12')
Chr20_map_win13 = map_maker(outlier_full_data$'13')
Chr20_ped_win13 = ped_maker(outlier_full_data$'13')


Chr20_data_win10 = Adegenet_PCA(outlier_ped = Chr20_ped_win10, 
                                outlier_map = Chr20_map_win10, 
                                OG_ped = OG_ped,
                                env = env_data)

Chr20_data_win11 = Adegenet_PCA(outlier_ped = Chr20_ped_win11, 
                                outlier_map = Chr20_map_win11, 
                                OG_ped = OG_ped,
                                env = env_data)
Chr20_data_win12 = Adegenet_PCA(outlier_ped = Chr20_ped_win12, 
                                outlier_map = Chr20_map_win12, 
                                OG_ped = OG_ped,
                                env = env_data)
Chr20_data_win13 = Adegenet_PCA(outlier_ped = Chr20_ped_win13, 
                                outlier_map = Chr20_map_win13, 
                                OG_ped = OG_ped,
                                env = env_data)


Chr20_data_win11 = Chr20_data_win11 %>% 
  select(contains('AX-'))
Chr20_data_win12 = Chr20_data_win12 %>% 
  select(contains('AX-'))
Chr20_data_win13 = Chr20_data_win13 %>% 
  select(contains('AX-'))


Chr20_win10111213 = bind_cols(Chr20_data_win10,
                          Chr20_data_win11, 
                          Chr20_data_win12, 
                          Chr20_data_win13)

Pop_that_pca(Chr20_win10111213, 
             pop_num = 39)

write_tsv(Chr20_win10111213, 
          '~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/Lostruct_Chr20_REGION1_win10111213_14.05.2021.txt')

## map file for continuous regions
Chr20_region1 = bind_rows(Chr20_map_win10, 
                          Chr20_map_win11, 
                          Chr20_map_win12, 
                          Chr20_map_win13)%>% 
  select(1:4) %>% 
  rename(`#Chromosome` = Chromosome) %>% 
  write_tsv('~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/CHR20_REGION1_14.05.2021.map')

head(Chr20_region1)
tail(Chr20_region1)

## region size
(27645135-22501825)/1000000


# CHR21 -------------------------------------------------------------------

lostruct_data = lostruct_run(data = tped, 
                             chr = 21, 
                             window_size = 20, 
                             k_value = 2)
outliers = Outlier_hunter(data = lostruct_data,
                          sd_percentile = 2)

outlier_full_data = Outlier_data(data = tped, 
                                 outlier_data = outliers, 
                                 chr = 21, 
                                 window_size = 20, 
                                 k_value = 2)

Chr21_map_win23 = map_maker(outlier_full_data$'23')
Chr21_ped_win23 = ped_maker(outlier_full_data$'23')
Chr21_map_win24 = map_maker(outlier_full_data$'24')
Chr21_ped_win24 = ped_maker(outlier_full_data$'24')
Chr21_map_win25 = map_maker(outlier_full_data$'25')
Chr21_ped_win25 = ped_maker(outlier_full_data$'25')

Chr21_data_win23 = Adegenet_PCA(outlier_ped = Chr21_ped_win23, 
                                outlier_map = Chr21_map_win23, 
                                OG_ped = OG_ped,
                                env = env_data)

Chr21_data_win24 = Adegenet_PCA(outlier_ped = Chr21_ped_win24, 
                                outlier_map = Chr21_map_win24, 
                                OG_ped = OG_ped,
                                env = env_data)
Chr21_data_win25 = Adegenet_PCA(outlier_ped = Chr21_ped_win25, 
                                outlier_map = Chr21_map_win25, 
                                OG_ped = OG_ped,
                                env = env_data)


Chr21_data_win24 = Chr21_data_win24 %>% 
  select(contains('AX-'))
Chr21_data_win25 = Chr21_data_win25 %>% 
  select(contains('AX-'))


Chr21_win232425 = bind_cols(Chr21_data_win23,
                              Chr21_data_win24, 
                              Chr21_data_win25)

Pop_that_pca(Chr21_win232425, 
             pop_num = 39)
## Definitely not an SV

# CHR22 -------------------------------------------------------------------
lostruct_data = lostruct_run(data = tped, 
                             chr = 22, 
                             window_size = 20, 
                             k_value = 2)
outliers = Outlier_hunter(data = lostruct_data,
                          sd_percentile = 2)

outlier_full_data = Outlier_data(data = tped, 
                                 outlier_data = outliers, 
                                 chr = 22, 
                                 window_size = 20, 
                                 k_value = 2)

Chr22_map_win10 = map_maker(outlier_full_data$'10')
Chr22_ped_win10 = ped_maker(outlier_full_data$'10')
Chr22_map_win11 = map_maker(outlier_full_data$'11')
Chr22_ped_win11 = ped_maker(outlier_full_data$'11')

Chr22_data_win10 = Adegenet_PCA(outlier_ped = Chr22_ped_win10, 
                                outlier_map = Chr22_map_win10, 
                                OG_ped = OG_ped,
                                env = env_data)

Chr22_data_win11 = Adegenet_PCA(outlier_ped = Chr22_ped_win11, 
                                outlier_map = Chr22_map_win11, 
                                OG_ped = OG_ped,
                                env = env_data)


Chr22_data_win11 = Chr22_data_win11 %>% 
  select(contains('AX-'))

Chr22_win1011 = bind_cols(Chr22_data_win10,
                            Chr22_data_win11)

Pop_that_pca(Chr22_win1011, 
             pop_num = 39)

write_tsv(Chr22_win1011, 
          '~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/Lostruct_Chr22_REGION1_win1011_14.05.2021.txt')

## map file for continuous regions
Chr22_region1 = bind_rows(Chr22_map_win10, 
                          Chr22_map_win11)%>% 
  select(1:4) %>% 
  rename(`#Chromosome` = Chromosome) %>% 
  write_tsv('~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/CHR22_REGION1_14.05.2021.map')

head(Chr22_region1)
tail(Chr22_region1)

## region size
(20025543-17288288)/1000000

# CHR23 -------------------------------------------------------------------
lostruct_data = lostruct_run(data = tped, 
                             chr = 23, 
                             window_size = 20, 
                             k_value = 2)
outliers = Outlier_hunter(data = lostruct_data,
                          sd_percentile = 2)

outlier_full_data = Outlier_data(data = tped, 
                                 outlier_data = outliers, 
                                 chr = 23, 
                                 window_size = 20, 
                                 k_value = 2)

# CHR23 REGIONS -----------------------------------------------------------
## combining windows 9, 10, 11
Chr23_map_win9 = map_maker(outlier_full_data$'9')
Chr23_ped_win9 = ped_maker(outlier_full_data$'9')
Chr23_map_win10 = map_maker(outlier_full_data$'10')
Chr23_ped_win10 = ped_maker(outlier_full_data$'10')
Chr23_map_win11 = map_maker(outlier_full_data$'11')
Chr23_ped_win11 = ped_maker(outlier_full_data$'11')

Chr23_data_win9 = Adegenet_PCA(outlier_ped = Chr23_ped_win9, 
                                outlier_map = Chr23_map_win9, 
                                OG_ped = OG_ped,
                                env = env_data)

Chr23_data_win10 = Adegenet_PCA(outlier_ped = Chr23_ped_win10, 
                                outlier_map = Chr23_map_win10, 
                                OG_ped = OG_ped,
                                env = env_data)

Chr23_data_win11 = Adegenet_PCA(outlier_ped = Chr23_ped_win11, 
                                outlier_map = Chr23_map_win11, 
                                OG_ped = OG_ped,
                                env = env_data)

Chr23_data_win10 = Chr23_data_win10 %>% 
  select(contains('AX-'))

Chr23_data_win11 = Chr23_data_win11 %>% 
  select(contains('AX-'))

Chr23_win91011 = bind_cols(Chr23_data_win9, 
                          Chr23_data_win10,
                          Chr23_data_win11)

Pop_that_pca(Chr23_win91011, 
             pop_num = 39)

write_tsv(Chr23_win91011, 
          '~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/Lostruct_CHR23_REGION1_win91011_14.05.2021.txt')

## map file for continuous regions
Chr23_region1 = bind_rows(Chr23_map_win9, 
                          Chr23_map_win10, 
                          Chr23_map_win11)%>% 
  select(1:4) %>% 
  rename(`#Chromosome` = Chromosome) %>% 
  write_tsv('~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/CHR23_REGION1_14.05.2021.map')

head(Chr23_region1)
tail(Chr23_region1)

## region size
(24975155-17149434)/1000000


## combine windows 29 and 30
Chr23_map_win29 = map_maker(outlier_full_data$'29')
Chr23_ped_win29 = ped_maker(outlier_full_data$'29')
Chr23_map_win30 = map_maker(outlier_full_data$'30')
Chr23_ped_win30 = ped_maker(outlier_full_data$'30')

Chr23_data_win29 = Adegenet_PCA(outlier_ped = Chr23_ped_win29, 
                               outlier_map = Chr23_map_win29, 
                               OG_ped = OG_ped,
                               env = env_data)

Chr23_data_win30 = Adegenet_PCA(outlier_ped = Chr23_ped_win30, 
                                outlier_map = Chr23_map_win30, 
                                OG_ped = OG_ped,
                                env = env_data)

Chr23_data_win30 = Chr23_data_win30 %>% 
  select(contains('AX-'))

Chr23_win2930 = bind_cols(Chr23_data_win29, 
                           Chr23_data_win30)

Pop_that_pca(Chr23_win2930, 
             pop_num = 39)

write_tsv(Chr23_win2930, 
          '~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/Lostruct_CHR23_REGION2_win2930_14.05.2021.txt')

## map file for continuous regions
Chr23_region2 = bind_rows(Chr23_map_win29, 
                          Chr23_map_win30)%>% 
  select(1:4) %>% 
  rename(`#Chromosome` = Chromosome) %>% 
  write_tsv('~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/CHR23_REGION2_14.05.2021.map')

head(Chr23_region2)
tail(Chr23_region2)

## region size
(67132585-63255761)/1000000


# CHR 24 ------------------------------------------------------------------
lostruct_data = lostruct_run(data = tped, 
                             chr = 24, 
                             window_size = 20, 
                             k_value = 2)
outliers = Outlier_hunter(data = lostruct_data,
                          sd_percentile = 2)

outlier_full_data = Outlier_data(data = tped, 
                                 outlier_data = outliers, 
                                 chr = 24, 
                                 window_size = 20, 
                                 k_value = 2)


# CHR25 -------------------------------------------------------------------

lostruct_data = lostruct_run(data = tped, 
                             chr = 25, 
                             window_size = 20, 
                             k_value = 2)
outliers = Outlier_hunter(data = lostruct_data,
                          sd_percentile = 2)

outlier_full_data = Outlier_data(data = tped, 
                                 outlier_data = outliers, 
                                 chr = 25, 
                                 window_size = 20, 
                                 k_value = 2)

Chr25_map_win13 = map_maker(outlier_full_data$'13')
Chr25_ped_win13 = ped_maker(outlier_full_data$'13')
Chr25_map_win14 = map_maker(outlier_full_data$'14')
Chr25_ped_win14 = ped_maker(outlier_full_data$'14')

Chr25_data_win13 = Adegenet_PCA(outlier_ped = Chr25_ped_win13, 
                                outlier_map = Chr25_map_win13, 
                                OG_ped = OG_ped,
                                env = env_data)

Chr25_data_win14 = Adegenet_PCA(outlier_ped = Chr25_ped_win14, 
                                outlier_map = Chr25_map_win14, 
                                OG_ped = OG_ped,
                                env = env_data)


Chr25_data_win14 = Chr25_data_win14 %>% 
  select(contains('AX-'))

Chr25_win1314 = bind_cols(Chr25_data_win13,
                          Chr25_data_win14)

Pop_that_pca(Chr25_win1314, 
             pop_num = 39)

## No dice Jim Rice

# CHR26 -------------------------------------------------------------------

lostruct_data = lostruct_run(data = tped, 
                             chr = 26, 
                             window_size = 20, 
                             k_value = 2)
outliers = Outlier_hunter(data = lostruct_data,
                          sd_percentile = 2)

outlier_full_data = Outlier_data(data = tped, 
                                 outlier_data = outliers, 
                                 chr = 26, 
                                 window_size = 20, 
                                 k_value = 2)


# CHR26 REGIONS -----------------------------------------------------------

Chr26_map_win1 = map_maker(outlier_full_data$'1')
Chr26_ped_win1 = ped_maker(outlier_full_data$'1')
Chr26_map_win2 = map_maker(outlier_full_data$'2')
Chr26_ped_win2 = ped_maker(outlier_full_data$'2')

Chr26_data_win1 = Adegenet_PCA(outlier_ped = Chr26_ped_win1, 
                                outlier_map = Chr26_map_win1, 
                                OG_ped = OG_ped,
                                env = env_data)

Chr26_data_win2 = Adegenet_PCA(outlier_ped = Chr26_ped_win2, 
                                outlier_map = Chr26_map_win2, 
                                OG_ped = OG_ped,
                                env = env_data)


Chr26_data_win2 = Chr26_data_win2 %>% 
  select(contains('AX-'))

Chr26_win12 = bind_cols(Chr26_data_win1,
                          Chr26_data_win2)

Pop_that_pca(Chr26_win12, 
             pop_num = 39)

## No dice Jim Rice


## combining windows 11 and 12
Chr26_map_win11 = map_maker(outlier_full_data$'11')
Chr26_ped_win11 = ped_maker(outlier_full_data$'11')
Chr26_map_win12 = map_maker(outlier_full_data$'12')
Chr26_ped_win12 = ped_maker(outlier_full_data$'12')

Chr26_data_win11 = Adegenet_PCA(outlier_ped = Chr26_ped_win11, 
                               outlier_map = Chr26_map_win11, 
                               OG_ped = OG_ped,
                               env = env_data)

Chr26_data_win12 = Adegenet_PCA(outlier_ped = Chr26_ped_win12, 
                               outlier_map = Chr26_map_win12, 
                               OG_ped = OG_ped,
                               env = env_data)


Chr26_data_win12 = Chr26_data_win12 %>% 
  select(contains('AX-'))

Chr26_win1112 = bind_cols(Chr26_data_win11,
                        Chr26_data_win12)

Pop_that_pca(Chr26_win1112, 
             pop_num = 39)
## Nah Bruh
# CHR27 -------------------------------------------------------------------

lostruct_data = lostruct_run(data = tped, 
                             chr = 27, 
                             window_size = 20, 
                             k_value = 2)
outliers = Outlier_hunter(data = lostruct_data,
                          sd_percentile = 2)

outlier_full_data = Outlier_data(data = tped, 
                                 outlier_data = outliers, 
                                 chr = 27, 
                                 window_size = 20, 
                                 k_value = 2)

# CHR28 -------------------------------------------------------------------

lostruct_data = lostruct_run(data = tped, 
                             chr = 28, 
                             window_size = 20, 
                             k_value = 2)
outliers = Outlier_hunter(data = lostruct_data,
                          sd_percentile = 2)

outlier_full_data = Outlier_data(data = tped, 
                                 outlier_data = outliers, 
                                 chr = 28, 
                                 window_size = 20, 
                                 k_value = 2)


# CHR 29 ------------------------------------------------------------------

lostruct_data = lostruct_run(data = tped, 
                             chr = 29, 
                             window_size = 20, 
                             k_value = 2)
outliers = Outlier_hunter(data = lostruct_data,
                          sd_percentile = 2)

outlier_full_data = Outlier_data(data = tped, 
                                 outlier_data = outliers, 
                                 chr = 29, 
                                 window_size = 20, 
                                 k_value = 2)

Chr29_map = map_maker(outlier_full_data$'16')
Chr29_ped = ped_maker(outlier_full_data$'16')


Chr29_data = Adegenet_PCA(outlier_ped = Chr29_ped, 
                         outlier_map = Chr29_map, 
                         OG_ped = OG_ped,
                         env = env_data)

Pop_that_pca(Chr29_data, 
             pop_num = 39,
             chr_num = 7, 
             win_num = 00)
## NAH

# CHR30 -------------------------------------------------------------------

lostruct_data = lostruct_run(data = tped, 
                             chr = 30, 
                             window_size = 20, 
                             k_value = 2)
outliers = Outlier_hunter(data = lostruct_data,
                          sd_percentile = 2)

outlier_full_data = Outlier_data(data = tped, 
                                 outlier_data = outliers, 
                                 chr = 30, 
                                 window_size = 20, 
                                 k_value = 2)

Chr30_map_win17 = map_maker(outlier_full_data$'17')
Chr30_ped_win17 = ped_maker(outlier_full_data$'17')
Chr30_map_win18 = map_maker(outlier_full_data$'18')
Chr30_ped_win18 = ped_maker(outlier_full_data$'18')

Chr30_data_win17 = Adegenet_PCA(outlier_ped = Chr30_ped_win17, 
                               outlier_map = Chr30_map_win17, 
                               OG_ped = OG_ped,
                               env = env_data)

Chr30_data_win18 = Adegenet_PCA(outlier_ped = Chr30_ped_win18, 
                               outlier_map = Chr30_map_win18, 
                               OG_ped = OG_ped,
                               env = env_data)


Chr30_data_win18 = Chr30_data_win18 %>% 
  select(contains('AX-'))

Chr30_win1718 = bind_cols(Chr30_data_win17,
                        Chr30_data_win18)

Pop_that_pca(Chr30_win1718, 
             pop_num = 39)
## NAH

# CHR 31 ------------------------------------------------------------------
lostruct_data = lostruct_run(data = tped, 
                             chr = 31, 
                             window_size = 20, 
                             k_value = 2)
outliers = Outlier_hunter(data = lostruct_data,
                          sd_percentile = 2)

outlier_full_data = Outlier_data(data = tped, 
                                 outlier_data = outliers, 
                                 chr = 31, 
                                 window_size = 20, 
                                 k_value = 2)
Chr31_map_win19 = map_maker(outlier_full_data$'19')
Chr31_ped_win19 = ped_maker(outlier_full_data$'19')
Chr31_map_win18 = map_maker(outlier_full_data$'18')
Chr31_ped_win18 = ped_maker(outlier_full_data$'18')
Chr31_map_win20 = map_maker(outlier_full_data$'20')
Chr31_ped_win20 = ped_maker(outlier_full_data$'20')


Chr31_data_win19 = Adegenet_PCA(outlier_ped = Chr31_ped_win19, 
                                outlier_map = Chr31_map_win19, 
                                OG_ped = OG_ped,
                                env = env_data)

Chr31_data_win18 = Adegenet_PCA(outlier_ped = Chr31_ped_win18, 
                                outlier_map = Chr31_map_win18, 
                                OG_ped = OG_ped,
                                env = env_data)

Chr31_data_win20 = Adegenet_PCA(outlier_ped = Chr31_ped_win20, 
                                outlier_map = Chr31_map_win20, 
                                OG_ped = OG_ped,
                                env = env_data)


Chr31_data_win19 = Chr31_data_win19 %>% 
  select(contains('AX-'))
Chr31_data_win20 = Chr31_data_win20 %>% 
  select(contains('AX-'))


Chr31_win181920 = bind_cols(Chr31_data_win18, 
                          Chr31_data_win19, 
                          Chr31_data_win20)

Pop_that_pca(Chr31_win181920, 
             pop_num = 39)
## NAH


# CHR 32 ------------------------------------------------------------------
lostruct_data = lostruct_run(data = tped, 
                             chr = 32, 
                             window_size = 20, 
                             k_value = 2)
outliers = Outlier_hunter(data = lostruct_data,
                          sd_percentile = 2)

outlier_full_data = Outlier_data(data = tped, 
                                 outlier_data = outliers, 
                                 chr = 32, 
                                 window_size = 20, 
                                 k_value = 2)


Chr32_map_win2 = map_maker(outlier_full_data$'2')
Chr32_ped_win2 = ped_maker(outlier_full_data$'2')
Chr32_data_win2 = Adegenet_PCA(outlier_ped = Chr32_ped_win2, 
                               outlier_map = Chr32_map_win2, 
                               OG_ped = OG_ped,
                               env = env_data)

Pop_that_pca(Chr32_data_win2,
             pop_num = 39)

write_tsv(Chr32_data_win2, 
          '~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/Lostruct_Chr32_REGION1_win2_14.05.2021.txt')

## map file for continuous regions
Chr32_region1 = Chr32_map_win2 %>% 
  select(1:4) %>% 
  rename(`#Chromosome` = Chromosome) %>% 
  write_tsv('~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/Chr32_REGION1_14.05.2021.map')

head(Chr32_region1)
tail(Chr32_region1)

## region size
(4504500-2008042)/1000000


# CHR33 -------------------------------------------------------------------

lostruct_data = lostruct_run(data = tped, 
                             chr = 33, 
                             window_size = 20, 
                             k_value = 2)
outliers = Outlier_hunter(data = lostruct_data,
                          sd_percentile = 2)

outlier_full_data = Outlier_data(data = tped, 
                                 outlier_data = outliers, 
                                 chr = 33, 
                                 window_size = 20, 
                                 k_value = 2)



# CHR34 -------------------------------------------------------------------

lostruct_data = lostruct_run(data = tped, 
                             chr = 34, 
                             window_size = 20, 
                             k_value = 2)
outliers = Outlier_hunter(data = lostruct_data,
                          sd_percentile = 2)

outlier_full_data = Outlier_data(data = tped, 
                                 outlier_data = outliers, 
                                 chr = 34, 
                                 window_size = 20, 
                                 k_value = 2)

## COmbine windows 5 & 6
Chr34_map_win22 = map_maker(outlier_full_data$'22')
Chr34_ped_win22 = ped_maker(outlier_full_data$'22')
Chr34_map_win23 = map_maker(outlier_full_data$'23')
Chr34_ped_win23 = ped_maker(outlier_full_data$'23')
Chr34_map_win24 = map_maker(outlier_full_data$'24')
Chr34_ped_win24 = ped_maker(outlier_full_data$'24')

Chr34_data_win22 = Adegenet_PCA(outlier_ped = Chr34_ped_win22, 
                               outlier_map = Chr34_map_win22, 
                               OG_ped = OG_ped,
                               env = env_data)

Chr34_data_win23 = Adegenet_PCA(outlier_ped = Chr34_ped_win23, 
                                outlier_map = Chr34_map_win23, 
                                OG_ped = OG_ped,
                                env = env_data)
Chr34_data_win24 = Adegenet_PCA(outlier_ped = Chr34_ped_win24, 
                                outlier_map = Chr34_map_win24, 
                                OG_ped = OG_ped,
                                env = env_data)


Chr34_data_win23 = Chr34_data_win23 %>% 
  select(contains('AX-'))
Chr34_data_win24 = Chr34_data_win24 %>% 
  select(contains('AX-'))


Chr34_win222324 = bind_cols(Chr34_data_win22, 
                           Chr34_data_win23, 
                         Chr34_data_win24)

Pop_that_pca(Chr34_win222324, 
             pop_num = 39)

## not super clear


# CHR35 -------------------------------------------------------------------
lostruct_data = lostruct_run(data = tped, 
                             chr = 35, 
                             window_size = 20, 
                             k_value = 2)
outliers = Outlier_hunter(data = lostruct_data,
                          sd_percentile = 2)

outlier_full_data = Outlier_data(data = tped, 
                                 outlier_data = outliers, 
                                 chr = 35, 
                                 window_size = 20, 
                                 k_value = 2)

## potential SV on window 2

## combine data on windows 5 and 6