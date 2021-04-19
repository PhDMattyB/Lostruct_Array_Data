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
library(patchwork)
library(lostruct)
library(adegenet)
library(viridis)

## set the working directory for this script
setwd('~/Charr_Adaptive_Introgression/Charr_Project_1/GeneticData/')

## load in the map file for the SNP array data set
map = read_tsv('Charr_Lab_recode12_25.03.2021.map', 
               col_names = c('Chromosome', 
                             'MarkerID', 
                             'Genetic_dist', 
                             'Physical_dist'))

## load in the ped file for the SNP array data set
OG_ped = read_table2('Charr_Lab_recode12_25.03.2021.ped', 
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

tped = Create_tped(ped = OG_ped, 
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
    filter(Chromosome == chr) %>% 
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
    as_tibble() %>% 
    rename(MDS_Points1 = V1, 
           MDS_Points2 = V2)
  
  combo_data = bind_cols(combo_data, 
                         MDS_points) %>% 
    # rename(MDS_Points1 = V1...25,
    #        MDS_Points2 = V2...26,
    #        V1 = V1...7,
    #        V2 = V2...8) %>% 
    dplyr::select(-MarkerID, 
                  -Genetic_dist, 
                  -Physical_dist, 
                  Chromosome, 
                  mean_window, 
                  window:length(.))
  # output = list(combo_data, 
  #               MDS_data)
  # 
  # return(output)
}


lostruct_data = lostruct_run(data = tped, 
             chr = 3, 
             window_size = 20, 
             k_value = 2)

lostruct_data$window

MDS_survey = function(data){
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

MDS_survey(lostruct_data)

Outlier_hunter = function(data){
  
  MDS1_outliers = data %>% 
    ungroup() %>% 
    mutate(MDS_cutoff = mean(MDS_Points1)+(2*sd(MDS_Points1))) %>% 
    filter(abs(MDS_Points1) > MDS_cutoff)
  MDS2_outliers = data %>% 
    ungroup() %>% 
    mutate(MDS_cutoff = mean(MDS_Points2)+(2*sd(MDS_Points2))) %>% 
    filter(abs(MDS_Points2) > MDS_cutoff)
  label = rep('MDS1 outlier', 
              nrow(MDS1_outliers)) %>% 
    as_tibble()
  MDS1_outliers = bind_cols(MDS1_outliers, 
                            label) %>%
    rename(outlier_lab = value)
  label = rep('MDS2 outlier', 
              nrow(MDS2_outliers)) %>% 
    as_tibble()
  MDS2_outliers = bind_cols(MDS2_outliers, 
                            label)%>% 
    rename(outlier_lab = value)
  Outliers = bind_rows(MDS1_outliers, 
                       MDS2_outliers)
  
  # 
  # MDS1_normy = data %>%
  #   ungroup() %>%
  #   mutate(MDS_cutoff = mean(MDS_Points1)+(2*sd(MDS_Points1))) %>%
  #   filter(abs(MDS_Points1) < MDS_cutoff)
  # label = rep('MDS1 non-outlier',
  #             nrow(MDS1_normy)) %>%
  #   as_tibble() %>%
  #   rename(outlier_lab = value)
  # MDS1_normy = bind_cols(MDS1_normy,
  #                        label)
  # MDS2_normy = data %>%
  #   ungroup() %>%
  #   mutate(MDS_cutoff = mean(MDS_Points2)+(2*sd(MDS_Points2))) %>%
  #   filter(abs(MDS_Points2) < MDS_cutoff)
  # label = rep('MDS2 non-outlier',
  #             nrow(MDS2_normy)) %>%
  #   as_tibble() %>%
  #   rename(outlier_lab = value)
  # MDS2_normy = bind_cols(MDS2_normy,
  #                        label)
  # 
  # # 
  # non_outliers = bind_rows(MDS1_normy,
  #                          MDS2_normy)%>%
  #   arrange(window) %>%
  #   distinct(window,
  #            .keep_all = T)
  # # 
  # outlier_df = bind_rows(non_outliers,
  #                        Outliers) %>%
  #   arrange(window)

  
  # outlier_df %>% 
  # group_by(window, 
  #          outlier_lab) %>% 
  #   filter(!duplicated(outlier_lab == 'MDS1 outlier') | outlier_lab != 'MDS1 outlier') %>%
  #   filter(!duplicated(outlier_lab == 'MDS2 outlier') | outlier_lab != 'MDS2 outlier') %>% 
  #   ungroup() %>% 
  #   distinct(window, 
  #            .keep_all = T)
  
}

outliers = Outlier_hunter(lostruct_data)


Outlier_plots = function(normal_data, 
                              outlier_data){
  
  theme_set(theme_bw())
  
  colour_pal = c('#F2055C',
                 '#1F26A6')
  
  plot1 = normal_data %>% 
    ggplot(aes(x = MDS_Points1, 
               y = MDS_Points2))+
    geom_point(size = 3, 
               col = 'Grey', 
               fill = 'Grey')+
    geom_point(data = outlier_data, 
               aes(x = MDS_Points1, 
                   y = MDS_Points2, 
                   col = outlier_lab, 
                   fill = outlier_lab), 
               size = 3)+
    scale_colour_manual(values = colour_pal)+
    labs(x = 'MDS coordinate 1', 
         y = 'MDS coordinate 2', 
         col = 'Outlier label', 
         fill = 'Outlier label')+
    theme(
      legend.title = element_text(size = 14), 
      legend.text = element_text(size = 12),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      axis.title = element_text(size = 14), 
      axis.text = element_text(size = 12)
    )
  
  
  outliers_mds1 = outlier_data %>% 
    filter(outlier_lab == 'MDS1 outlier')
  outliers_mds2 = outlier_data %>% 
    filter(outlier_lab == 'MDS2 outlier')
  
  plot2 = normal_data %>% 
    ggplot(aes(x = mean_window, 
               y = abs(MDS_Points1)))+
    geom_point(col = 'Grey', 
               size = 3)+
    geom_point(data = outliers_mds1, 
               aes(x = mean_window, 
                   y = abs(MDS_Points1)),
               col = '#F2055C',
               size = 3)+
    labs(x = 'Mean window position (bp)', 
         y = 'MDS score (axis 1)')+
    theme(
      legend.position = 'none',
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      axis.title = element_text(size = 14), 
      axis.text = element_text(size = 12)
    )
  
  
 plot3 = normal_data %>% 
    ggplot(aes(x = mean_window, 
               y = abs(MDS_Points2)))+
    geom_point(col = 'Grey', 
               size = 3)+
    geom_point(data = outliers_mds2, 
               aes(x = mean_window, 
                   y = abs(MDS_Points2)),
               col = '#1F26A6',
               size = 3)+
    labs(x = 'Mean window position (bp)', 
         y = 'MDS score (axis 2)')+
    theme(
      legend.position = 'none',
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      axis.title = element_text(size = 14), 
      axis.text = element_text(size = 12)
    )
  
 combo = plot1/(plot2 + plot3)
    
    return(combo)
}

Outlier_plots(normal_data = lostruct_data, 
              outlier_data = outliers)

##

## Next up. Need to peforma PCA on those regions to see
## if they're actually structural variants. Going to need
## the FULL genomic information for the outlier regions. 

Outlier_data = function(data,
                            outlier_data,
                            chr,
                            window_size, 
                            k_value){
  df = data %>% 
    dplyr::select(Chromosome, 
                  5:length(tped)) %>% 
    filter(Chromosome == chr) %>% 
    dplyr::select(-Chromosome)
  eigen = eigen_windows(df, 
                        win = window_size, 
                        k = k_value)
  windist = pc_dist(eigen, 
                    npc = k_value) %>% 
    as_tibble()
  
  tped_data = data %>% 
    # select(1:4) %>% 
    filter(Chromosome == chr) %>% 
    mutate(window = ceiling(row_number()/window_size)) %>% 
    group_by(window) %>% 
    mutate(mean_window = mean(Physical_dist)) %>% 
    filter(window %in% outliers$window) %>% 
    dplyr::select(Chromosome, 
           MarkerID,
           Genetic_dist, 
           Physical_dist, 
           window, 
           mean_window, 
           contains('_')) 
  
  ## Splitting by window to get a dataframe for each outlier window
  by_window = split(tped_data, tped_data$window) 
  
}

 outlier_full_data = Outlier_data(data = tped, 
                 outlier_data = outliers, 
                 chr = 2, 
                 window_size = 20, 
                 k_value = 2)

 map_maker = function(data){
   make_map = data %>% 
     ungroup() %>% 
     select(Chromosome,
            MarkerID,
            Genetic_dist,
            Physical_dist,
            window,
            mean_window)
   
   
}

 ped_maker = function(data){
   marker_names = data %>% 
     ungroup() %>% 
     select(MarkerID) %>% 
     t() %>% 
     as_tibble() %>% 
     row_to_names(row_number = 1) %>% 
     names()
   
   make_ped = data %>% 
     ungroup() %>% 
     select(contains('_'), 
            -Genetic_dist, 
            -Physical_dist, 
            -mean_window) %>% 
     t() %>% 
     as_tibble() %>% 
     rename_all(funs(c(marker_names)))
   
 }

 outlier_full_data
 
 env_data = read_csv('~/Charr_Adaptive_Introgression/Charr_Project_1/SampleSiteData/SampleSites_Coords_1June2020.csv')
 
 Chr1_win17_map = map_maker(outlier_full_data$'17')
 Chr1_win17_ped = ped_maker(outlier_full_data$'17')
 
 Chr2_map = map_maker(outlier_full_data$'23')
 Chr2_ped = ped_maker(outlier_full_data$'23')
 
 
 ## Now we need to run a PCA in adegenet!!
 
Adegenet_PCA = function(outlier_ped, 
                        outlier_map,
                        OG_ped,
                        env){
  
   env_data = env %>% 
     rename(FamilyID = Population)
   
   ped_data = OG_ped %>% 
     dplyr::select(FamilyID, 
                   IndividualID)%>% 
     filter(FamilyID != 'GDL') 
   
   metadata = left_join(ped_data, 
                        env_data, 
                        by = 'FamilyID')
   
   ped = bind_cols(metadata, 
                   outlier_ped)
   
   data = as.data.frame(ped)
   
   genind = df2genind(data[,8:length(data)], 
                          ploidy = 2, 
                          ind.names = data[,2],
                          sep = '\t', 
                          pop = data[,1]
                      )
   pca_data = tab(genind, 
       freq=TRUE, 
       NA.method="mean")
   
   pca = dudi.pca(df=pca_data,
            center = T, 
            scale = F,
            nf = 2, 
            scannf = F)
   
   pca_data = pca$li %>% 
     as_tibble()
   
   ped = ped %>% 
     as_tibble() %>% 
     bind_cols(., 
               pca_data)
   
   ggplot(data = ped, 
          aes(x = Axis1, 
              y = Axis2))+
     geom_point(aes(col = FamilyID))
   
   return(ped)
 }

# chr1_win17_data = Adegenet_PCA(outlier_ped = Chr1_win17_ped, 
#                               outlier_map = Chr1_win17_map, 
#                               OG_ped = ped,
#                               env = env_data)


chr2_data = Adegenet_PCA(outlier_ped = Chr2_ped, 
                               outlier_map = Chr2_map, 
                               OG_ped = ped,
                               env = env_data)

# theme_set(theme_bw())
Pop_that_pca = function(data, 
                        pop_num,
                        chr_num,
                        win_num){
  
  data %>% 
    ggplot(aes(x = Axis1, 
               y = Axis2))+
    geom_point(aes(col = FamilyID), 
               size = 2)+
    labs(x = 'PCA axis 1', 
         y = 'PCA axis 2', 
         title = paste0('Chr', 
                        chr_num,
                        " ",
                        'Win',
                        win_num,
                        " ",
                        'PCA plot'))+
    scale_color_manual(values = magma(n = pop_num))+
    theme(
      legend.position = 'none', 
      panel.grid = element_blank(), 
      axis.title = element_text(size = 14), 
      axis.text = element_text(size = 12)
    )
    
}

Pop_that_pca(chr2_data, 
             pop_num = 37,
             chr_num = 2, 
             win_num = 23)
