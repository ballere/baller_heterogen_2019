require('visreg')
require('mgcv')
require('tableone')
require('dplyr')
require('plm')
require('MatchIt')
require('tidyr')
require('ggplot2')
require('reshape')
require('emmeans')
require('cowplot')
require('stringr')
require('visreg')
require('rasterVis')
require('lattice')
require('circlize')
## set directories ##
#local_wkdir <- '~/Google Drive/TDSlab/SCZ_gene_imaging/'
#remote_wkdir <- '~/Desktop/BBL/data/joy/BBL/studies/pnc/'

#############################
###Make Demographics Table###
#############################
#############################
####### Demographics ########
#############################

#######Matched group all clusters ##############

make_demographics_table<- function(data_frame, hydra_cluster) {
  #subset demographics
  cluster <- paste0("Hydra_k", hydra_cluster)
  listVars <- c("Race", "Sex", "Maternal Ed", "Age", "Depression", "Cluster") #Race 1 = caucasian, Maternal Ed = years, age = years, dep 1 = dep, 0 = non_dep
  for_parse <- paste0("data.frame(data_frame$race_binarized, data_frame$sex, data_frame$medu1, data_frame$age_in_years, data_frame$dep_binarized, data_frame$", cluster, ")")
  demo <- eval(parse(text = for_parse)) 
  names(demo) <- c(listVars)

  #Change categorical values to have names
  demo$Depression <- ifelse(demo$Depression == 1, "Depressed", "Non-depressed")
  demo$Race <- ifelse(demo$Race == 1, "Caucasian", "Non-caucasian")
  demo$Sex <- ifelse(demo$Sex == 1, "Male", "Female")
 # demo$Cluster <- ifthen(demo$Cluster == -1, "TD")
#  demo$Cluster <- ifthen(demo$Cluster == 1, "Cluster 1")
#  demo$Cluster <- ifthen(demo$Cluster == 2, "Cluster 2")
#  demo$Cluster <- ifthen(demo$Cluster == 1, "Cluster 3")
  
  #Define Categorical Variables
  cat_variables <- c("Race", "Depression", "Sex", "Cluster")
  title <- c(paste0("Hydra_k", hydra_cluster, "demographics"))

  #create demographics table
  demo_table <- CreateTableOne(vars = listVars, data = demo, factorVars = cat_variables, strata = c("Cluster"))
  print(demo_table, showAllLevels = TRUE)
}

total_people_per_cluster <- function(data_frame, hydra_cluster)
{
  #returns vector with a number for each cluster
  if (hydra_cluster > 1 & hydra_cluster < 10) {
    length_controls <- eval(parse(text = paste0("length(which(data_frame$Hydra_k", hydra_cluster, " == -1))")))
    cluster_num_vector <- c(length_controls)
    for (cluster_counter in 1:hydra_cluster){
      length_cluster <- eval(parse(text = paste0("length(which(data_frame$Hydra_k", hydra_cluster, " == ", cluster_counter, "))")))
      cluster_num_vector <- c(cluster_num_vector, length_cluster)
    }
    return(cluster_num_vector)
  }
  else {
    print("Error: Cluster number must be between 2 and 10")
  }
}

total_people_per_cluster_by_group <- function(data_frame, variable, hydra_cluster, group_val)
{
  #returns vector with a number for each cluster per group
  #make sure variable is in quotes
  #group_val is the group you'd like. For exampke, if variable= sex and group_val = 1, it will send back number of men, 
       #if variable = race_binarized and group_val = 1, it will send back #caucasians
  if (hydra_cluster > 1 & hydra_cluster < 10) {
    length_controls <- eval(parse(text = paste0("length(which(data_frame$Hydra_k", hydra_cluster, " == -1 & data_frame$", variable, " == ", group_val, "))")))
    cluster_num_vector <- c(length_controls)
    for (cluster_counter in 1:hydra_cluster){
      length_cluster <- eval(parse(text = paste0("length(which(data_frame$Hydra_k", hydra_cluster, " == ", cluster_counter, " & data_frame$", variable, " == ", group_val, "))")))
      cluster_num_vector <- c(cluster_num_vector, length_cluster)
    }
    return(cluster_num_vector)
  }
  else {
    print("Error: Cluster number must be between 2 and 10")
  }
}

chi_sq <- function(data_frame, variable, hydra_cluster) {
  #make sure variable is in quotes
  chisq <- eval(parse(text = paste0("chisq.test(data_frame$", variable, ", data_frame$Hydra_k", hydra_cluster, ")")))
  return(chisq)
}

chi_sq_no_controls <- function(data_frame, variable, hydra_cluster)
{
  chisq <- eval(parse(text = paste0("chisq.test(data_frame$", variable, ", (data_frame$Hydra_k", hydra_cluster, " != -1))")))
  return(chisq)
}

chi_sq_p <- function(data_frame, variable, hydra_cluster) {
 #make sure variable is in quotes
  chisq <- eval(parse(text = paste0("chisq.test(data_frame$", variable, ", data_frame$Hydra_k", hydra_cluster, ")")))
  print(chisq)
  chisq_pvalue <- chisq$p.value
  return(chisq_pvalue)
}

get_cluster_titles <- function(hydra_cluster){
  # build vector of titles
  cluster_titles <- c("TD")
  for (cluster_counter in 1:hydra_cluster){
    title_to_add <- paste0("Subtype ", cluster_counter)
    cluster_titles <- c(cluster_titles, title_to_add)
  }
  return(cluster_titles)
}
#get_cluster_titles <- function(hydra_cluster){
  # build vector of titles
 # cluster_titles <- c("TD")
  #for (cluster_counter in 1:hydra_cluster){
   # title_to_add <- paste0("Cluster ", cluster_counter)
    #cluster_titles <- c(cluster_titles, title_to_add)
#  }
 # return(cluster_titles)
#}
get_cluster_titles_no_TD <- function(hydra_cluster){
  # build vector of titles
  cluster_titles <- NULL
  for (cluster_counter in 1:hydra_cluster){
    title_to_add <- paste0("Cluster ", cluster_counter)
    cluster_titles <- c(cluster_titles, title_to_add)
  }
  return(cluster_titles)
}

get_cluster_numerical_vector <- function(hydra_cluster){
  # build vector of titles
  cluster_vector <- c("-1")
  for (cluster_counter in 1:hydra_cluster){
    title_to_add <- paste0(cluster_counter)
    cluster_vector <- c(cluster_vector, title_to_add)
  }
  return(cluster_vector)
}
 
get_cluster_numerical_vector_no_TD <- function(hydra_cluster){
  # build vector of titles
  cluster_vector <- NULL
  for (cluster_counter in 1:hydra_cluster){
    title_to_add <- paste0(cluster_counter)
    cluster_vector <- c(cluster_vector, title_to_add)
  }
  return(cluster_vector)
}

get_variable_mean_vector <- function(data_frame, variable, hydra_cluster){
  #returns a vector of means, depends on # hydra clusters
  means <- eval(parse(text = paste0("c(mean(data_frame$", variable, "[which(data_frame$Hydra_k", hydra_cluster, " == -1)]))")))
  for (cluster_counter in 1:hydra_cluster){
    mean_to_add <- eval(parse(text = paste0("c(mean(data_frame$", variable, "[which(data_frame$Hydra_k", hydra_cluster, " == ", cluster_counter, ")]))")))
    means <- c(means, mean_to_add)
  }
  return(means)
}

get_variable_mean_vector_no_TD <- function(data_frame, variable, hydra_cluster){
  #returns a vector of means, depends on # hydra clusters
  #means <- eval(parse(text = paste0("c(mean(data_frame$", variable, "[which(data_frame$Hydra_k", hydra_cluster, " == -1)]))")))
  means <- NULL
  for (cluster_counter in 1:hydra_cluster){
    mean_to_add <- eval(parse(text = paste0("c(mean(data_frame$", variable, "[which(data_frame$Hydra_k", hydra_cluster, " == ", cluster_counter, ")]))")))
    means <- c(means, mean_to_add)
  }
  return(means)
}

get_variable_sd_vector <- function(data_frame, variable, hydra_cluster){
  #returns vector of sds, depends on # hydra clusters
  sds <- eval(parse(text = paste0("c(sd(data_frame$", variable, "[which(data_frame$Hydra_k", hydra_cluster, " == -1)]))")))
  for (cluster_counter in 1:hydra_cluster){
    sd_to_add <- eval(parse(text = paste0("c(sd(data_frame$", variable, "[which(data_frame$Hydra_k", hydra_cluster, " == ", cluster_counter, ")]))")))
    sds <- c(sds, sd_to_add)
  }
  return(sds)
}

get_variable_sd_vector_no_TD <- function(data_frame, variable, hydra_cluster){
  #returns vector of sds, depends on # hydra clusters
  #sds <- eval(parse(text = paste0("c(sd(data_frame$", variable, "[which(data_frame$Hydra_k", hydra_cluster, " == -1)]))")))
  sds <- NULL
  for (cluster_counter in 1:hydra_cluster){
    sd_to_add <- eval(parse(text = paste0("c(sd(data_frame$", variable, "[which(data_frame$Hydra_k", hydra_cluster, " == ", cluster_counter, ")]))")))
    sds <- c(sds, sd_to_add)
  }
  return(sds)
}

get_num_subj_per_cluster <- function(data_frame, hydra_cluster)
{
  nums <- eval(parse(text = paste0("c(length(which(data_frame$Hydra_k", hydra_cluster, " == -1)))")))
  for (cluster_counter in 1:hydra_cluster){
    num_to_add <- eval(parse(text = paste0("c(length(which(data_frame$Hydra_k", hydra_cluster, " == ", cluster_counter, ")))")))
    nums <- c(nums, num_to_add)
  }
  return(nums)
}

get_num_subj_per_cluster_no_TD <- function(data_frame, hydra_cluster)
{
 # nums <- eval(parse(text = paste0("c(length(which(data_frame$Hydra_k", hydra_cluster, " == -1)))")))
  nums <- NULL
  for (cluster_counter in 1:hydra_cluster){
    num_to_add <- eval(parse(text = paste0("c(length(which(data_frame$Hydra_k", hydra_cluster, " == ", cluster_counter, ")))")))
    nums <- c(nums, num_to_add)
  }
  return(nums)
}

data_frame_mean_sd_sem <- function(data_frame, variable, hydra_cluster){
  cluster_titles <- get_cluster_titles(hydra_cluster = hydra_cluster)
  variable_mean <- get_variable_mean_vector(data_frame = data_frame, variable = variable, hydra_cluster = hydra_cluster)
  variable_sd <- get_variable_sd_vector(data_frame = data_frame, variable = variable, hydra_cluster = hydra_cluster)
  num_per_cluster <- get_num_subj_per_cluster(data_frame = data_frame, hydra_cluster = hydra_cluster)
  variable_sem <- variable_sd/sqrt(num_per_cluster)
  
  #put all together in one data frame
  df_mean_sd_sem <- data.frame(cl = cluster_titles, mean = variable_mean, sd = variable_sd, sem = variable_sem)
  return(df_mean_sd_sem)
}

data_frame_mean_sd_sem_no_TD <- function(data_frame, variable, hydra_cluster){
  cluster_titles <- get_cluster_titles_no_TD(hydra_cluster = hydra_cluster)
  variable_mean <- get_variable_mean_vector_no_TD(data_frame = data_frame, variable = variable, hydra_cluster = hydra_cluster)
  variable_sd <- get_variable_sd_vector_no_TD(data_frame = data_frame, variable = variable, hydra_cluster = hydra_cluster)
  num_per_cluster <- get_num_subj_per_cluster_no_TD(data_frame = data_frame, hydra_cluster = hydra_cluster)
  variable_sem <- variable_sd/sqrt(num_per_cluster)
  
  #put all together in one data frame
  df_mean_sd_sem <- data.frame(cl = cluster_titles, mean = variable_mean, sd = variable_sd, sem = variable_sem)
  return(df_mean_sd_sem)
}



###### Plotting #####
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

plot_continuous_variables <- function(data_frame, var1, var2, hydra_cluster, optional_variable_name_string){
#makes bar plots and line plot with error bars, specifically for plotting numerical things like age and medu. 
  #NOT for categorical variables like sex/race 
  
  num_total_groups <- hydra_cluster + 1
  cluster_titles <- get_cluster_titles(hydra_cluster = hydra_cluster)
  dat_var1_sd_sem <- data_frame_mean_sd_sem(data_frame = data_frame, variable = var1, hydra_cluster = hydra_cluster)
  dat_var2_sd_sem <- data_frame_mean_sd_sem(data_frame = data_frame, variable = var2, hydra_cluster = hydra_cluster)
  dat_var1_var2_sem_melted <- melt(c(dat_var1_sd_sem$sem,dat_var2_sd_sem$sem), id.vars = "cl")


  #get everything into the right format
  var1_and_var2 <- data.frame(cl=cluster_titles, var1 = dat_var1_sd_sem$mean, var2 = dat_var2_sd_sem$mean)
  var1_and_var2_for_plot <- melt(var1_and_var2, id.vars = "cl")
  var1_and_var2_for_plot$sem <- dat_var1_var2_sem_melted$value
  names(var1_and_var2_for_plot) <- c("cluster", "group", "years", "sem")
  
  #change names from var1 and var2 to their actual values
  if (!missing(optional_variable_name_string))
  {
    var1 = optional_variable_name_string[1]
    var2 = optional_variable_name_string[2]
  }
  
  replace_group_names <- c(rep(var1, num_total_groups), rep(var2, num_total_groups))
  #replace_group_names <- c(rep(var1, 4), rep(var2, 4))
  var1_and_var2_for_plot$group <- replace_group_names

  #plot
  title_of_plot <- paste0("Hydra_k", hydra_cluster, " ", var1, " and ", var2)
  p1 <- ggplot(data = var1_and_var2_for_plot, aes(x = group, y = years, group = cluster)) + 
    geom_line(aes(color=cluster)) +
    geom_point(aes(color=cluster)) + 
    geom_errorbar(aes(ymin=years-sem, ymax=years+sem), width=.1) +
    ggtitle(title_of_plot)

  p2 <- ggplot(dat_var1_sd_sem, aes(x = cl, y = mean, fill = cl)) + geom_col() + 
    geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem),width=.2,position=position_dodge(.9)) + 
    scale_x_discrete(limits=cluster_titles) + ylim(0, 18) + xlab("Clusters") + ylab(paste0(var1, " in Years")) + 
    ggtitle(paste0(var1, " by Cluster")) + scale_fill_discrete(breaks=cluster_titles) +
    guides(fill=guide_legend(title=NULL)) 
  
  p3 <- ggplot(dat_var2_sd_sem, aes(x = cl, y = mean, fill = cl)) + geom_col() + 
    geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem),width=.2,position=position_dodge(.9)) + 
    scale_x_discrete(limits=cluster_titles) + ylim(0, 18) + xlab("Clusters") + ylab(paste0(var2, " in Years")) + 
    ggtitle(paste0(var2, " by Cluster")) + scale_fill_discrete(breaks=cluster_titles) +
    guides(fill=guide_legend(title=NULL)) 
    
  #send back all three to the rmarkdown document that I called
  list_to_return <- list(p1, p2, p3)
   # multiplot(p1,p2,p3, cols=3)
  return(list_to_return)
  
}

better_levelplot <- function(adj, node_names, title) {
  adj_norm <- adj/max(abs(adj))
  limit = max(abs(adj_norm))
  #keycol=c('#FFFDE7','#FFFF00', '#F57F17', '#D50000',"#212121","#311B92","#2979FF","#B2EBF2")
  #keycol=c("d11225", "", "", "", "", "", "", "")
  plot<-levelplot(adj_norm, par.settings = RdBuTheme(region = brewer.pal(9, 'Reds')), #RdBuTheme(), 
                  colorkey=list(labels=list(cex=2,font=1,col="black")),
                  at = seq(limit,0,length.out = 12),xlab="",ylab = "", strip = F, contour = F, region= T,main=list(label=title,cex=3),
                  scales=list(x=list(at = 1:length(node_names), labels=node_names,rot=90, tck = 0, cex = 2),
                              y=list(at = 1:length(node_names),labels=node_names, tck = 0, cex = 2)))
  return(plot)
}

## better_level_plot ##
better_levelplot_edge <- function(adj, node_names_x, node_names_y, title) {
  adj_norm <- adj/max(abs(adj))
  limit = max(abs(adj_norm))
  keycol=c('#FFFDE7','#FFFF00', '#F57F17', '#D50000',"#212121","#311B92","#2979FF","#B2EBF2")
  plot<-levelplot(adj_norm, par.settings = BuRdTheme(), 
                  at = seq(limit,-limit,length.out = 12),xlab="",ylab = "", strip = F, contour = F, region= T,main=title,
                  scales=list(x=list(at = 1:length(node_names_x), labels=node_names_x,rot=90, tck = 0),
                              y=list(at = 1:length(node_names_y),labels=node_names_y, tck = 0)))
  return(plot)
}

better_chorDiagram_norm <- function(adj, node_names) {
  circos.clear()
  model_color=c('#FFFDE7','#FFFF00', '#F57F17', '#D50000',"#212121","#311B92","#2979FF","#B2EBF2")
  adj_norm <- adj/max(abs(adj))
  rownames(adj_norm) = node_names
  colnames(adj_norm) = node_names
  chordDiagram(adj_norm, directional = TRUE, transparency = 0.5, self.link = TRUE, grid.col = model_color) 
}

better_chorDiagram <- function(adj, node_names) {
  circos.clear()
  model_color=c('#CC2626','#CC7926', '#CCC326', '#26CCB8',"#2689CC","#3126CC","#A526CC","#CC2665")
  rownames(adj) = node_names
  colnames(adj) = node_names
  par(cex=2)
  chordDiagram(adj, directional = TRUE, transparency = 0.5, self.link = TRUE)#, grid.col = model_color) 
}

better_chorDiagram_13 <- function(adj, node_names) {
  circos.clear()
  model_color=c('#CC2626','#CC6D26', '#CCA526', '#C9CC26', '#8FCC26', '#36CC26', '#26CC8F', '#268CCC', '#262EEC', '#8426CC', '#CC26CC', '#CC2676', '#75696F')
  rownames(adj) = node_names
  colnames(adj) = node_names
  par(cex=1.5)
  chordDiagram(adj, directional = TRUE, transparency = 0.5, self.link = TRUE, grid.col = model_color) 
}

#### Stats#####
fdr_anova <- function(data_frame) {
  models_anova <- lapply(data_frame, summary)
  
  #Pull p-values
  p_anova <- sapply(data_frame, function(v) v$"Pr(>F)"[1]) #$coef[,"Pr(>F)"][2]) #get the p value for dep binarized
  
  #Convert to data frame
  p_anova <- as.data.frame(p_anova)
  
  #print BEFORE FDR correction 
  #print("Anova scores, BEFORE FDR correction: i.e., uncorrected")
  print(p_anova)
  
  #Print original p-values to three decimal places
  p_round_anova <- round(p_anova,3)
  
  #FDR correct p-values
  pfdr_anova <- p.adjust(p_anova[,1],method="fdr")
  
  #Convert to data frame
  pfdr_anova <- as.data.frame(pfdr_anova)
  row.names(pfdr_anova) <- names(data_frame)
  
  #To print fdr-corrected p-values to three decimal places
  pfdr_round_anova <- round(pfdr_anova,3)
  
  #List the components that survive FDR correction
  components_fdr_anova <- row.names(pfdr_anova)[pfdr_anova<0.05]
  
  #make a data frame with names and fdr values (rounded to 3 decimals)
  names_and_fdr_values_anova <- data.frame(cbind(components_fdr_anova, round(pfdr_anova[pfdr_anova<0.05],3)))
  
  #add titles to names_and_fdr tables
  names(names_and_fdr_values_anova) <- c("component", "p_FDR_corr")
  
 # print("Mean centered age that was then squared, FDR corrected")
  #print(names_and_fdr_values_anova)
  return(names_and_fdr_values_anova)
  
}

pairwise_contrasts_3clusters <- function(data_frame_lm, fdr_anova) {
 #for 3 clusters
  print(fdr_anova)
  emmodel_df <- lapply(data_frame_lm, function(x) {as.list(ref_grid(x))})
  emgrid_df <- lapply(emmodel_df, function(x) {as.emmGrid(x)})
  
  #run emmeans
  emmeans_df <- lapply(emgrid_df, function(x) {emmeans(x, "Hydra_k3")})
  
  #run pairwise contrasts
  empairs_df <- lapply(emmeans_df, function(x) {pairs(x)})
  

  #Only include stuff that was fdr corrected (i.e., only keep parts of the model (or only display) ones that are corrected),this will be null if nothing was corrected
  empairs_FDR_corrected <- empairs_df[fdr_anova[,1]]
  
  print(empairs_FDR_corrected)
  #contrast names, -1 = controls, 1-3 are clusters
  contrast_names <- c("-1 - 1", "-1 - 2", "-1 - 3", "1 - 2", "1 - 3", "2 - 3")
  #go through each fdr corrected brain region, and extract p values
  #contrast_table <- lapply(fdr_anova[,1], function(x) {round(summary(x)$p.value,3)})
  contrast_table <- lapply(empairs_FDR_corrected, function(x) {round(summary(x)$p.value,3)})
  #get the names of the brain regions that were fdr corrected
  brain_regions <- names(contrast_table)
  #build table that will hold the name of the brain region and the p values
  pairwise_table <- data.frame(matrix(nrow = length(brain_regions), ncol = 6))
  #give the appropriate names
  rownames(pairwise_table) <- brain_regions
  colnames(pairwise_table) <- contrast_names
  
  #loop through each brain region, and manually assign the columns to be the p values
  for (region in brain_regions)
  {
    pair_pval <- contrast_table[[region]]
    pairwise_table[region,] <- pair_pval
  }
  pairwise_table_with_fdr <- pairwise_table
  pairwise_table_with_fdr$p_FDR_corr <- fdr_anova$p_FDR_corr
  
  pairwise <- list(empairs_df, empairs_FDR_corrected, pairwise_table, pairwise_table_with_fdr)
  return(pairwise)
  
}

make_pairwise_contrast_names <- function(num_clusters) {
  #get the numerical vector per number of clusters, will loop through them, and generate vector with all pairs
  clusters_vector1 <- get_cluster_numerical_vector(hydra_cluster = num_clusters)
  clusters_vector2 <- clusters_vector1
  pairs_vector <- NULL
  for(pair1 in clusters_vector1){
    for(pair2 in clusters_vector2)
      if (pair1 < pair2) {
        pair_text <- paste0(pair1, " - ", pair2)
        pairs_vector <- c(pairs_vector, pair_text)
      }
  }
  return(pairs_vector)
}


pairwise_contrasts_generic_num_clusters <- function(data_frame_lm, fdr_anova, num_clusters) {
  #for 3 clusters
  print(fdr_anova)
  emmodel_df <- lapply(data_frame_lm, function(x) {as.list(ref_grid(x))})
  emgrid_df <- lapply(emmodel_df, function(x) {as.emmGrid(x)})
  
  #run emmeans
  hydra_var <- paste0("Hydra_k", num_clusters)
  emmeans_df <- lapply(emgrid_df, function(x) {emmeans(x, hydra_var)})
  
  #run pairwise contrasts
  empairs_df <- lapply(emmeans_df, function(x) {pairs(x)})
  
  
  #Only include stuff that was fdr corrected (i.e., only keep parts of the model (or only display) ones that are corrected),this will be null if nothing was corrected
  empairs_FDR_corrected <- empairs_df[fdr_anova[,1]]
  
 # print(empairs_FDR_corrected)
  #contrast names, -1 = controls, other numbers are clusters
  contrast_names <- make_pairwise_contrast_names(num_clusters = num_clusters)
  #go through each fdr corrected brain region, and extract p values
  #contrast_table <- lapply(fdr_anova[,1], function(x) {round(summary(x)$p.value,3)})
  contrast_table <- lapply(empairs_FDR_corrected, function(x) {round(summary(x)$p.value,3)})
  #get the names of the brain regions that were fdr corrected
  brain_regions <- names(contrast_table)
  #build table that will hold the name of the brain region and the p values
  #number of columns (or number of unique pairs is)
  pairwise_table <- data.frame(matrix(nrow = length(brain_regions), ncol = length(contrast_names)))
  #give the appropriate names
  rownames(pairwise_table) <- brain_regions
  colnames(pairwise_table) <- contrast_names
  
  #loop through each brain region, and manually assign the columns to be the p values
  for (region in brain_regions)
  {
    pair_pval <- contrast_table[[region]]
    pairwise_table[region,] <- pair_pval
  }
  pairwise_table_with_fdr <- pairwise_table
  pairwise_table_with_fdr$p_FDR_corr <- fdr_anova$p_FDR_corr
  
  pairwise <- list(empairs_df, empairs_FDR_corrected, pairwise_table, pairwise_table_with_fdr)
  return(pairwise)
  
}



############## Get community names #############
get_community_names <- function(num_communities) {
  #takes a number of communities, returns a list of each community, within, and between network labels
  if (num_communities == 7) {
    com_names <- c("visual", "somatomotor", "dorsalAttention", "salienceVentralAttention", "limbic", "frontoparietalControl", "default")
    within <- c("visual_visual", "somatomotor_somatomotor", "dorsalAttention_dorsalAttention", "salienceVentralAttention_salienceVentralAttention", "limbic_limbic", "frontoparietalControl_frontoparietalControl", "default_default")
    between <- c("visual_somatomotor", "visual_dorsalAttention", "somatomotor_dorsalAttention", "visual_salienceVentralAttention", "somatomotor_salienceVentralAttention", "dorsalAttention_salienceVentralAttention", "visual_limbic", "somatomotor_limbic", "dorsalAttention_limbic",                       
    "salienceVentralAttention_limbic", "visual_frontoparietalControl",                  
    "somatomotor_frontoparietalControl", "dorsalAttention_frontoparietalControl",         
    "salienceVentralAttention_frontoparietalControl", "limbic_frontoparietalControl",                  
    "visual_default", "somatomotor_default",                          
    "dorsalAttention_default", "salienceVentralAttention_default",              
    "limbic_default", "frontoparietalControl_default")     
     } else if (num_communities == 17) {
       networks <- read.csv("/Users/eballer/BBL/from_chead/ballerDepHeterogen/results/csvs/community_names.csv")
       com_names <- networks[2]
       com_names <- c(sapply(com_names, as.character))
       within <- networks[2]
       within <- c(sapply(within, as.character))
      
      #make list of between
       num_between <- (num_communities * (num_communities -1))/2
       between <- array(data = "NA", dim = num_between)
       nrow <- 1
       num_netA <- 1
       num_netB <- 1
       for (network_a in within) {
         for (network_b in within) {
           #make sure you aren't adding betworks that are either WITHIN, or have already been added to the matrix
           if ((network_a != network_b) & (num_netA < num_netB)) {
              between[nrow] <- paste0(network_a, "_", network_b)
              nrow <- nrow + 1
           }
           num_netB <- num_netB + 1
         }
         num_netA <- num_netA + 1
         num_netB <- 1
       }
       between <- c(sapply(between, as.character))
       #make the character lists of within networks, aka visualCentral_visualCentral
       for (row in 1:length(within)) {
         within[row] <- paste0(within[row], "_", within[row])
       }
        
  } else {
    com_names <- c("error")
    within <- c("error")
    between <- c("error")
  }
  all_names<- list(com_names, within, between)
  return(all_names)
}
