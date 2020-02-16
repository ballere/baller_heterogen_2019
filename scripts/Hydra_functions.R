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

  #Define Categorical Variables
  cat_variables <- c("Race", "Depression", "Sex", "Cluster")
  title <- c("Hydra_k3 demographics")

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
    title_to_add <- paste0("Cluster", cluster_counter)
    cluster_titles <- c(cluster_titles, title_to_add)
  }
  return(cluster_titles)
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

get_variable_sd_vector <- function(data_frame, variable, hydra_cluster){
  #returns vector of sds, depends on # hydra clusters
  sds <- eval(parse(text = paste0("c(sd(data_frame$", variable, "[which(data_frame$Hydra_k", hydra_cluster, " == -1)]))")))
  for (cluster_counter in 1:hydra_cluster){
    sd_to_add <- eval(parse(text = paste0("c(sd(data_frame$", variable, "[which(data_frame$Hydra_k", hydra_cluster, " == ", cluster_counter, ")]))")))
    sds <- c(sds, sd_to_add)
  }
  return(sds)
}

data_frame_mean_sd_sem <- function(data_frame, variable, hydra_cluster){
  cluster_titles <- get_cluster_titles(hydra_cluster = hydra_cluster)
  variable_mean <- get_variable_mean_vector(data_frame = data_frame, variable = variable, hydra_cluster = hydra_cluster)
  variable_sd <- get_variable_sd_vector(data_frame = data_frame, variable = variable, hydra_cluster = hydra_cluster)
  variable_sem <- variable_sd/sqrt(nrow(data_frame))
  
  #put all together in one data frame
  df_mean_sd_sem <- data.frame(cl = cluster_titles, mean = variable_mean, sd = variable_sd, sem = variable_sem)
  return(df_mean_sd_sem)
}

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
  replace_group_names <- c(rep(var1, 4), rep(var2, 4))
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


