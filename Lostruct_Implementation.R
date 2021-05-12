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
OG_ped = read_table2('Charr_Lab_recode12_25.03.2021.ped',
                     col_names = c('FamilyID',
                                   'IndividualID',
                                   'PaternalID',
                                   'MaternalID',
                                   'Sex',
                                   'Phenotype',
                                   map$MarkerID))

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
                            'UNH')) %>% 
  rename(Chromosome = 1)

tped = Create_tped(ped = OG_ped, 
                   map = map) 
# Chr 1 -------------------------------------------------------------------
lostruct_data = lostruct_run(data = tped, 
                             chr = 1, 
                             window_size = 20, 
                             k_value = 2)
chr = 1
window_size = 20
k_value = 2
df = tped %>% 
  dplyr::select(Chromosome, 
                5:length(tped)) %>% 
  filter(Chromosome == chr) %>% 
  dplyr::select(-Chromosome) 

## WHY WONT THIS FUCKING WORK!!!!!!
## IT worked before!!!

eigen = eigen_windows(df, 
                      win = window_size, 
                      k = k_value)

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


## getting the data for chr1 region 1
Chr1_map_win7 = map_maker(outlier_full_data$'7')
Chr1_ped_win7 = ped_maker(outlier_full_data$'7')


Chr1_win7_data = Adegenet_PCA(outlier_ped = Chr1_ped_win7, 
                              outlier_map = Chr1_map_win7, 
                              OG_ped = OG_ped,
                              env = env_data)

write_tsv(Chr1_win7_data, 
          '~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/CHR1_REGION1_lostruct.txt')

PCA_chr1_win7 = Pop_that_pca(Chr1_win7_data, 
                             pop_num = 38,
                             chr_num = 11, 
                             win_num = 7)

## combine data for continuous regions
chr1_region1_outlier_map = Chr1_map_win7 %>% 
  select(1:4) %>% 
  rename(`#Chromosome` = Chromosome) %>% 
  write_tsv('~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/Chr1_region1_win7_outlier_windows.map')
head(chr1_region1_outlier_map)
tail(chr1_region1_outlier_map)

## region size 
region_size = 17710378 - 15748822
region_size/1000000
## CHR1 REGION 1 is 1.96 Mb

## CHR 1 REGION 2
## Combining continuous outlier regions
## combining windows 7 and 8
chr_map_win16 = map_maker(outlier_full_data$'16')
chr_ped_win16 = ped_maker(outlier_full_data$'16')
chr_data_win16 = Adegenet_PCA(outlier_ped = chr_ped_win16, 
                              outlier_map = chr_map_win16, 
                              OG_ped = OG_ped,
                              env = env_data)

chr_map_win17 = map_maker(outlier_full_data$'17')
chr_ped_win17 = ped_maker(outlier_full_data$'17')
chr_data_win17 = Adegenet_PCA(outlier_ped = chr_ped_win17, 
                              outlier_map = chr_map_win17, 
                              OG_ped = OG_ped,
                              env = env_data)



chr_data_win17 = chr_data_win17 %>% 
  select(contains('AX-'))

chr1_combo_win1617 = bind_cols(chr_data_win16, 
                               chr_data_win17)

PCA_outlier_wins1617 = Pop_that_pca(chr1_combo_win1617, 
                                    pop_num = 38,
                                    chr_num = 1, 
                                    win_num = 1617)

write_tsv(chr1_combo_win1617, 
          '~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/Lostruct_chr1_region2.txt')
## write map files for the combo region

chr1_region2_outlier_map = bind_rows(chr_map_win16, 
                                     chr_map_win17) %>% 
  select(1:4) %>% 
  rename(`#Chromosome` = Chromosome, 
         `Marker ID` = MarkerID, 
         `Genetic distance` = Genetic_dist, 
         `Physical distance` = Physical_dist) %>% 
  write_tsv('~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/Chr1_region2_outlier_windows.map')

head(chr1_region2_outlier_map)
tail(chr1_region2_outlier_map)

## Calculate regions size
region_size = 32789668-30591434
region_size/1000000
## The region size of the outlier is
## 2.2 Megabases with 40 SNPs. 
## SNP density is DEFINITELY and issue in identifying 
## these structural variants. 

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


outlier_full_data$'5' %>% 
  names() %>% 
  view()

Chr2_map_win23 = map_maker(outlier_full_data$'23')
Chr2_ped_win23 = ped_maker(outlier_full_data$'23')


Chr2_win23_data = Adegenet_PCA(outlier_ped = Chr2_ped_win23, 
                               outlier_map = Chr2_map_win23, 
                               OG_ped = OG_ped,
                               env = env_data)

write_tsv(Chr2_win23_data,
          '~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/CHR2_REGION3_win23_lostruct.txt')

Pop_that_pca(Chr2_win23_data, 
             pop_num = 38,
             chr_num = 2, 
             win_num = 23)

## Regions with only 20 snps. Not continuous regions
write_tsv(Chr2_map_win20, 
          '~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/CHR2_REGION3_Win20.map')
head(Chr2_map_win20)
tail(Chr2_map_win20)
## make ped file from data_explore R script


## COMBO CONTINUOUS REGIONS!!!
## Need to combine regions 22 and 23. They may make up 
## a larger region along the chromosome

## CHR 1 REGION 2
## Combining continuous outlier regions

Chr2_win22_data
Chr2_win23_data = Chr2_win23_data %>% 
  select(contains('AX-'))

Chr2_win2223 = bind_cols(Chr2_win22_data, 
                         Chr2_win23_data)

Chr2_win2223 %>% 
  select(FamilyID) %>% 
  distinct() %>% 
  arrange(FamilyID) %>% 
  View()

Pop_that_pca(Chr2_win2223, 
             pop_num = 38)

write_tsv(Chr2_win2223, 
          '~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/Lostruct_CHR2_REGION4_win2223_lostruct.txt')

## map file for continuous regions
Chr2_region4 = bind_rows(Chr2_map_win22, 
          Chr2_map_win23)%>% 
  select(1:4) %>% 
  rename(`#Chromosome` = Chromosome) %>% 
  write_tsv('~/Charr_Adaptive_Introgression/Charr_Project_1/Lostruc/CHR2_REGION4_win2223.map')

head(Chr2_region4)
tail(Chr2_region4)

## Calculate regions size
region_size = 32789668-30591434
region_size/1000000
## The region size of the outlier is
## 2.2 Megabases with 40 SNPs. 
## SNP density is DEFINITELY and issue in identifying 
## these structural variants. 

## To make the combo map ped file
## Need to use the data object in the Data_explore script
## The data needed for adegenet is recoded as 12 format
## needs to be in AA CC GG TT format

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
