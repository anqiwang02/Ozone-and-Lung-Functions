#fisher interval for secondary genes

#clear the environment
rm(list=ls()) 

#load library

library(dplyr)
library(tidyr)
library(parallel)
library(ggplot2)

n_cores <- suppressWarnings(as.integer(Sys.getenv("FISHER_N_CORES", unset = NA)))
if (is.na(n_cores) || n_cores < 1L) {
  n_cores <- max(1L, parallel::detectCores() - 1L)
}

# Save all PDF graphs under ozone/fisher_int_sec_graph_<today's date>/
graph_out_dir <- paste0("fisher_int_sec_graph_", format(Sys.Date(), "%Y-%m-%d"))
dir.create(graph_out_dir, recursive = TRUE, showWarnings = FALSE)

#Section 1: cg12662193 of PLSCR1. effect sizes=0.00891

load("../Exp_data_secondary_CpG.RData")
load("../../Data/w0true.Rdata")
 


mn <- which(colnames(final_secondary_CpG)=="cg12662193")   #provide the column number for the CpG site 

#all possible randomizations
n <- 2^17

#create a data frame that contains the CpG site of interest
final1<-
  final_secondary_CpG %>%
  select(c(1:4),mn) %>%
  rename(Y_obs = cg12662193)

#data wrangling
#ozone, 1=O3
final2 <-
  final1%>%
  filter(exp==1) %>%
  rename(Yi_wiO3 = Y_obs)%>%
  select(3,4,5)

#clean air, -1=CA
final3 <-
  final1 %>%
  filter(exp==-1)%>%
  rename(Yi_wiCA = Y_obs)%>%
  select(3,4,5)


#full join table final4 with final2 or final 3, respectively
final4 <-
  full_join(final1, final2, by = c("id2" = "id2", "exp" = "exp")) %>%
  select(1,2,3,4,6,5)

final4<-
  full_join(final4, final3, by = c("id2" = "id2", "exp" = "exp")) %>%
  select(1,2,3,4,5,7,6)
 


#for lower boundry. This CpG site has a positive effect size -->calculating p value should use >=
#should have no absolute value when calculating p value. -0.01
 
a=200
vectora_low = as.matrix(seq(-0.005, 0.017, length.out = a)) #if changing the number of a values here length.out=# of a values

x<-list()

for(i in 1:a) { # of a value in the loop
  #substitute the NA value in Yi_1 with the correponding value in Yi_0 but "a" more
  final4 <-
    full_join(final1, final2, by = c("id2" = "id2", "exp" = "exp")) %>%
    select(1,2,3,4,6,5)
  
  final4<-
    full_join(final4, final3, by = c("id2" = "id2", "exp" = "exp")) %>%
    select(1,2,3,4,5,7,6)
  
  #substitute the NA value in Yi_1 with the correponding value in Yi_0 but a more
  final4$Yi_wiO3[is.na(final4$Yi_wiO3)] <- final4$Yi_wiCA[final4$exp==-1]+vectora_low[i,1]
  final4$Yi_wiCA[is.na(final4$Yi_wiCA)] <- final4$Yi_wiO3[final4$exp==1]-vectora_low[i,1]
  x[[i]]<-final4
}
 


 
#create a matrix with all p-values corresponding to each value a


list_p_low <- matrix(nrow=a, ncol=1)  # nrow= number of a values

# Function for processing a single `i`
compute_p_value <- function(i) {
  id <- x[[i]]$"id2"
  W_obs <- x[[i]]$"exp"
  Yi_CA <- x[[i]]$"Yi_wiCA"
  Yi_O3 <- x[[i]]$"Yi_wiO3"
  Y_obs <- x[[i]]$"Y_obs"
  
  # Compute statistics for all j in one step
  stat <- apply(w0.true, 1, function(w) { #1 represents applying each row in W0.true
    t.test(Yi_O3[w == 1], Yi_CA[w == -1],paired=TRUE)$statistic
  })
  
  # Observed statistics (no absolute value)
  stat_obs <- t.test(Y_obs[W_obs == 1], Y_obs[W_obs == -1],paired=TRUE)$statistic
  
  # Compute p-value
  p_value <- 1*sum(stat>=stat_obs)/n
  
  return(p_value)
}

# Parallel execution
cl <- makeCluster(n_cores)
clusterExport(cl, c("x", "w0.true", "n"))  # Added "n" to the export list
list_p_low[, 1] <- parSapply(cl, 1:a, compute_p_value)
stopCluster(cl)

 




# for upper boundry. This CpG site has a positive effect size -->calculating p value should use <=
 
a=200
vectora = as.matrix(seq(0.002, 0.016888889, length.out = a)) #if changing the number of a values here length.out=# of a values

x<-list()

for(i in 1:a) { # of a value in the loop
  #substitute the NA value in Yi_1 with the correponding value in Yi_0 but "a" more
  final4 <-
    full_join(final1, final2, by = c("id2" = "id2", "exp" = "exp")) %>%
    select(1,2,3,4,6,5)
  
  final4<-
    full_join(final4, final3, by = c("id2" = "id2", "exp" = "exp")) %>%
    select(1,2,3,4,5,7,6)
  
  #substitute the NA value in Yi_1 with the correponding value in Yi_0 but a more
  final4$Yi_wiO3[is.na(final4$Yi_wiO3)] <- final4$Yi_wiCA[final4$exp==-1]+vectora[i,1]
  final4$Yi_wiCA[is.na(final4$Yi_wiCA)] <- final4$Yi_wiO3[final4$exp==1]-vectora[i,1]
  x[[i]]<-final4
}
 


 
#create a matrix with all p-values corresponding to each value a
list_p<- matrix(nrow=a, ncol=1)  # nrow= number of a values

# Function for processing a single `i`
compute_p_value <- function(i) {
  id <- x[[i]]$"id2"
  W_obs <- x[[i]]$"exp"
  Yi_CA <- x[[i]]$"Yi_wiCA"
  Yi_O3 <- x[[i]]$"Yi_wiO3"
  Y_obs <- x[[i]]$"Y_obs"
  
  # Compute statistics for all j in one step
  stat <- apply(w0.true, 1, function(w) { #1 represents applying each row in W0.true
    t.test(Yi_O3[w == 1], Yi_CA[w == -1], paired=TRUE)$statistic
  })
  
  # Observed statistics (no absolute value)
  stat_obs <- t.test(Y_obs[W_obs == 1], Y_obs[W_obs == -1], paired=TRUE)$statistic
  
  # Compute p-value
  p_value <- 1*sum(stat<=stat_obs)/n
  
  return(p_value)
}

# Parallel execution
cl <- makeCluster(n_cores)
clusterExport(cl, c("x", "w0.true", "n"))  # Added "n" to the export list
list_p[, 1] <- parSapply(cl, 1:a, compute_p_value)
stopCluster(cl)
 


# visualization 
# graphing of lower boundry (df1)
 
pdf(file.path(graph_out_dir, "cg12662193_fish_int_low.pdf"), width = 8, height = 6, pointsize = 12, family = "Times")  # Specify the file name

#lower boundry 

# Assuming vectora_low and list_p_low are vectors
df1 <- data.frame(vectora_low, list_p_low = list_p_low[,1])  # Convert to dataframe

ggplot(df1, aes(x = vectora_low, y = list_p_low)) +
  geom_point(color = "blue", size = 2, alpha = 0.6) +  # Scatter plot
  geom_hline(yintercept = 0.025, linetype = "dashed", color = "red", size = 1) +  # Horizontal line
  labs(
    title = "P-values across different hypothetical constant treatment effect for cg12662193",
    x = "Hypothetical constant treatment effect (a)",
    y = "p-value"
  ) +
  theme_classic()   # Use a clean theme

dev.off()

 

#graphing of upper boundry (df2)
 
pdf(file.path(graph_out_dir, "cg12662193_fish_int_up.pdf"), width = 8, height = 6, pointsize = 12, family = "Times")  # Specify the file name
#upper boundry 

# Assuming vectora_low and list_p_low are vectors
df2 <- data.frame(vectora, list_p = list_p[,1])  # Convert to dataframe

ggplot(df2, aes(x = vectora, y = list_p)) +
  geom_point(color = "blue", size = 2, alpha = 0.6) +  # Scatter plot
  geom_hline(yintercept = 0.025, linetype = "dashed", color = "red", size = 1) +  # Horizontal line
  labs(
    title = "P-values across different hypothetical constant treatment effect for cg12662193",
    x = "Hypothetical constant treatment effect (a)",
    y = "p-value"
  ) +
  theme_classic()   # Use a clean theme

dev.off()
 

# confidence CI: [0.004473684, 0.01342105]



# section 2: cg19869469 of C6orf227. effect size= 0.017787267
 

mn <- which(colnames(final_secondary_CpG)=="cg19869469")   #provide the column number for the CpG site 

#all possible randomizations
n <- 2^17

#create a data frame that contains the CpG site of interest
final1<-
  final_secondary_CpG %>%
  select(c(1:4),mn) %>%
  rename(Y_obs = cg19869469 )

#data wrangling
#ozone, 1=O3
final2 <-
  final1%>%
  filter(exp==1) %>%
  rename(Yi_wiO3 = Y_obs)%>%
  select(3,4,5)

#clean air, -1=CA
final3 <-
  final1 %>%
  filter(exp==-1)%>%
  rename(Yi_wiCA = Y_obs)%>%
  select(3,4,5)


#full join table final4 with final2 or final 3, respectively
final4 <-
  full_join(final1, final2, by = c("id2" = "id2", "exp" = "exp")) %>%
  select(1,2,3,4,6,5)

final4<-
  full_join(final4, final3, by = c("id2" = "id2", "exp" = "exp")) %>%
  select(1,2,3,4,5,7,6)
 


#for lower boundry. This CpG site has a positive effect size -->calculating p value should use >=
#should have no absolute value when calculating p value
 
a=200
vectora_low = as.matrix(seq(-0.01, 0.033, length.out = a)) #if changing the number of a values here length.out=# of a values

x<-list()

for(i in 1:a) { # of a value in the loop
  #substitute the NA value in Yi_1 with the correponding value in Yi_0 but "a" more
  final4 <-
    full_join(final1, final2, by = c("id2" = "id2", "exp" = "exp")) %>%
    select(1,2,3,4,6,5)
  
  final4<-
    full_join(final4, final3, by = c("id2" = "id2", "exp" = "exp")) %>%
    select(1,2,3,4,5,7,6)
  
  #substitute the NA value in Yi_1 with the correponding value in Yi_0 but a more
  final4$Yi_wiO3[is.na(final4$Yi_wiO3)] <- final4$Yi_wiCA[final4$exp==-1]+vectora_low[i,1]
  final4$Yi_wiCA[is.na(final4$Yi_wiCA)] <- final4$Yi_wiO3[final4$exp==1]-vectora_low[i,1]
  x[[i]]<-final4
}
 


 
#create a matrix with all p-values corresponding to each value a


list_p_low <- matrix(nrow=a, ncol=1)  # nrow= number of a values

# Function for processing a single `i`
compute_p_value <- function(i) {
  id <- x[[i]]$"id2"
  W_obs <- x[[i]]$"exp"
  Yi_CA <- x[[i]]$"Yi_wiCA"
  Yi_O3 <- x[[i]]$"Yi_wiO3"
  Y_obs <- x[[i]]$"Y_obs"
  
  # Compute statistics for all j in one step
  stat <- apply(w0.true, 1, function(w) { #1 represents applying each row in W0.true
    t.test(Yi_O3[w == 1], Yi_CA[w == -1],paired=TRUE)$statistic
  })
  
  # Observed statistics (no absolute value)
  stat_obs <- t.test(Y_obs[W_obs == 1], Y_obs[W_obs == -1],paired=TRUE)$statistic
  
  # Compute p-value
  p_value <- 1*sum(stat>=stat_obs)/n
  
  return(p_value)
}

# Parallel execution
cl <- makeCluster(n_cores)
clusterExport(cl, c("x", "w0.true", "n"))  # Added "n" to the export list
list_p_low[, 1] <- parSapply(cl, 1:a, compute_p_value)
stopCluster(cl)

 


# for upper boundry. This CpG site has a positive effect size -->calculating p value should use <=
 
a=200
vectora = as.matrix(seq(0.005, 0.03131578, length.out = a)) #if changing the number of a values here length.out=# of a values

x<-list()

for(i in 1:a) { # of a value in the loop
  #substitute the NA value in Yi_1 with the correponding value in Yi_0 but "a" more
  final4 <-
    full_join(final1, final2, by = c("id2" = "id2", "exp" = "exp")) %>%
    select(1,2,3,4,6,5)
  
  final4<-
    full_join(final4, final3, by = c("id2" = "id2", "exp" = "exp")) %>%
    select(1,2,3,4,5,7,6)
  
  #substitute the NA value in Yi_1 with the correponding value in Yi_0 but a more
  final4$Yi_wiO3[is.na(final4$Yi_wiO3)] <- final4$Yi_wiCA[final4$exp==-1]+vectora[i,1]
  final4$Yi_wiCA[is.na(final4$Yi_wiCA)] <- final4$Yi_wiO3[final4$exp==1]-vectora[i,1]
  x[[i]]<-final4
}
 


 
#create a matrix with all p-values corresponding to each value a
list_p<- matrix(nrow=a, ncol=1)  # nrow= number of a values

# Function for processing a single `i`
compute_p_value <- function(i) {
  id <- x[[i]]$"id2"
  W_obs <- x[[i]]$"exp"
  Yi_CA <- x[[i]]$"Yi_wiCA"
  Yi_O3 <- x[[i]]$"Yi_wiO3"
  Y_obs <- x[[i]]$"Y_obs"
  
  # Compute statistics for all j in one step
  stat <- apply(w0.true, 1, function(w) { #1 represents applying each row in W0.true
    t.test(Yi_O3[w == 1], Yi_CA[w == -1], paired=TRUE)$statistic
  })
  
  # Observed statistics (no absolute value)
  stat_obs <- t.test(Y_obs[W_obs == 1], Y_obs[W_obs == -1], paired=TRUE)$statistic
  
  # Compute p-value
  p_value <- 1*sum(stat<=stat_obs)/n
  
  return(p_value)
}

# Parallel execution
cl <- makeCluster(n_cores)
clusterExport(cl, c("x", "w0.true", "n"))  # Added "n" to the export list
list_p[, 1] <- parSapply(cl, 1:a, compute_p_value)
stopCluster(cl)
 


# visualization 
# graphing of lower boundry (df3)
 
pdf(file.path(graph_out_dir, "cg19869469_fish_int_low.pdf"), width = 8, height = 6, pointsize = 12, family = "Times")  # Specify the file name

#lower boundry 

# Assuming vectora_low and list_p_low are vectors
df3 <- data.frame(vectora_low, list_p_low = list_p_low[,1])  # Convert to dataframe

ggplot(df3, aes(x = vectora_low, y = list_p_low)) +
  geom_point(color = "blue", size = 2, alpha = 0.6) +  # Scatter plot
  geom_hline(yintercept = 0.025, linetype = "dashed", color = "red", size = 1) +  # Horizontal line
  labs(
    title = "P-values across different hypothetical constant treatment effect for cg19869469",
    x = "Hypothetical constant treatment effect (a)",
    y = "p-value"
  ) +
  theme_classic()   # Use a clean theme

dev.off()

 

#graphing of upper boundry (df4)
 
pdf(file.path(graph_out_dir, "cg19869469_fish_int_up.pdf"), width = 8, height = 6, pointsize = 12, family = "Times")  # Specify the file name
#upper boundry 

# Assuming vectora_low and list_p_low are vectors
df4 <- data.frame(vectora, list_p = list_p[,1])  # Convert to dataframe

ggplot(df4, aes(x = vectora, y = list_p)) +
  geom_point(color = "blue", size = 2, alpha = 0.6) +  # Scatter plot
  geom_hline(yintercept = 0.025, linetype = "dashed", color = "red", size = 1) +  # Horizontal line
  labs(
    title = "P-values across different hypothetical constant treatment effect for cg19869469",
    x = "Hypothetical constant treatment effect (a)",
    y = "p-value"
  ) +
  theme_classic()   # Use a clean theme

dev.off()
 

# confidence interval: [0.009368421, 0.02615236]



# section 3: cg14795253 of PLSCR1. effect size = 0.006089086
 
mn <- which(colnames(final_secondary_CpG)=="cg14795253")   #provide the column number for the CpG site 

#all possible randomizations
n <- 2^17

#create a data frame that contains the CpG site of interest
final1<-
  final_secondary_CpG %>%
  select(c(1:4),mn) %>%
  rename(Y_obs = cg14795253)

#data wrangling
#ozone, 1=O3
final2 <-
  final1%>%
  filter(exp==1) %>%
  rename(Yi_wiO3 = Y_obs)%>%
  select(3,4,5)

#clean air, -1=CA
final3 <-
  final1 %>%
  filter(exp==-1)%>%
  rename(Yi_wiCA = Y_obs)%>%
  select(3,4,5)


#full join table final4 with final2 or final 3, respectively
final4 <-
  full_join(final1, final2, by = c("id2" = "id2", "exp" = "exp")) %>%
  select(1,2,3,4,6,5)

final4<-
  full_join(final4, final3, by = c("id2" = "id2", "exp" = "exp")) %>%
  select(1,2,3,4,5,7,6)
 


#for lower boundry. This CpG site has a positive effect size -->calculating p value should use >=
#should have no absolute value when calculating p value
 
a=200
vectora_low = as.matrix(seq(-0.004, 0.013, length.out = a)) #if changing the number of a values here length.out=# of a values

x<-list()

for(i in 1:a) { # of a value in the loop
  #substitute the NA value in Yi_1 with the correponding value in Yi_0 but "a" more
  final4 <-
    full_join(final1, final2, by = c("id2" = "id2", "exp" = "exp")) %>%
    select(1,2,3,4,6,5)
  
  final4<-
    full_join(final4, final3, by = c("id2" = "id2", "exp" = "exp")) %>%
    select(1,2,3,4,5,7,6)
  
  #substitute the NA value in Yi_1 with the correponding value in Yi_0 but a more
  final4$Yi_wiO3[is.na(final4$Yi_wiO3)] <- final4$Yi_wiCA[final4$exp==-1]+vectora_low[i,1]
  final4$Yi_wiCA[is.na(final4$Yi_wiCA)] <- final4$Yi_wiO3[final4$exp==1]-vectora_low[i,1]
  x[[i]]<-final4
}
 


 

list_p_low <- matrix(nrow=a, ncol=1)  # nrow= number of a values

# Function for processing a single `i`
compute_p_value <- function(i) {
  id <- x[[i]]$"id2"
  W_obs <- x[[i]]$"exp"
  Yi_CA <- x[[i]]$"Yi_wiCA"
  Yi_O3 <- x[[i]]$"Yi_wiO3"
  Y_obs <- x[[i]]$"Y_obs"
  
  # Compute statistics for all j in one step
  stat <- apply(w0.true, 1, function(w) { #1 represents applying each row in W0.true
    t.test(Yi_O3[w == 1], Yi_CA[w == -1],paired=TRUE)$statistic
  })
  
  # Observed statistics (no absolute value)
  stat_obs <- t.test(Y_obs[W_obs == 1], Y_obs[W_obs == -1],paired=TRUE)$statistic
  
  # Compute p-value
  p_value <- 1*sum(stat>=stat_obs)/n
  
  return(p_value)
}

# Parallel execution
cl <- makeCluster(n_cores)
clusterExport(cl, c("x", "w0.true", "n"))  # Added "n" to the export list
list_p_low[, 1] <- parSapply(cl, 1:a, compute_p_value)
stopCluster(cl)

 

# for upper boundry. This CpG site has a positive effect size -->calculating p value should use <=
 
a=200
vectora = as.matrix(seq(0.001, 0.012, length.out = a)) #if changing the number of a values here length.out=# of a values

x<-list()

for(i in 1:a) { # of a value in the loop
  #substitute the NA value in Yi_1 with the correponding value in Yi_0 but "a" more
  final4 <-
    full_join(final1, final2, by = c("id2" = "id2", "exp" = "exp")) %>%
    select(1,2,3,4,6,5)
  
  final4<-
    full_join(final4, final3, by = c("id2" = "id2", "exp" = "exp")) %>%
    select(1,2,3,4,5,7,6)
  
  #substitute the NA value in Yi_1 with the correponding value in Yi_0 but a more
  final4$Yi_wiO3[is.na(final4$Yi_wiO3)] <- final4$Yi_wiCA[final4$exp==-1]+vectora[i,1]
  final4$Yi_wiCA[is.na(final4$Yi_wiCA)] <- final4$Yi_wiO3[final4$exp==1]-vectora[i,1]
  x[[i]]<-final4
}
 


 
#create a matrix with all p-values corresponding to each value a
list_p<- matrix(nrow=a, ncol=1)  # nrow= number of a values

# Function for processing a single `i`
compute_p_value <- function(i) {
  id <- x[[i]]$"id2"
  W_obs <- x[[i]]$"exp"
  Yi_CA <- x[[i]]$"Yi_wiCA"
  Yi_O3 <- x[[i]]$"Yi_wiO3"
  Y_obs <- x[[i]]$"Y_obs"
  
  # Compute statistics for all j in one step
  stat <- apply(w0.true, 1, function(w) { #1 represents applying each row in W0.true
    t.test(Yi_O3[w == 1], Yi_CA[w == -1], paired=TRUE)$statistic
  })
  
  # Observed statistics (no absolute value)
  stat_obs <- t.test(Y_obs[W_obs == 1], Y_obs[W_obs == -1], paired=TRUE)$statistic
  
  # Compute p-value
  p_value <- 1*sum(stat<=stat_obs)/n
  
  return(p_value)
}

# Parallel execution
cl <- makeCluster(n_cores)
clusterExport(cl, c("x", "w0.true", "n"))  # Added "n" to the export list
list_p[, 1] <- parSapply(cl, 1:a, compute_p_value)
stopCluster(cl)
 


# visualization 
# graphing of lower boundry (df5)
 
pdf(file.path(graph_out_dir, "cg14795253_fish_int_low.pdf"), width = 8, height = 6, pointsize = 12, family = "Times")  # Specify the file name

#lower boundry 

# Assuming vectora_low and list_p_low are vectors
df5 <- data.frame(vectora_low, list_p_low = list_p_low[,1])  # Convert to dataframe

ggplot(df5, aes(x = vectora_low, y = list_p_low)) +
  geom_point(color = "blue", size = 2, alpha = 0.6) +  # Scatter plot
  geom_hline(yintercept = 0.025, linetype = "dashed", color = "red", size = 1) +  # Horizontal line
  labs(
    title = "P-values across different hypothetical constant treatment effect for cg14795253",
    x = "Hypothetical constant treatment effect (a)",
    y = "p-value"
  ) +
  theme_classic()   # Use a clean theme

dev.off()

 

#graphing of upper boundry (df6)
 
pdf(file.path(graph_out_dir, "cg14795253_fish_int_up.pdf"), width = 8, height = 6, pointsize = 12, family = "Times")  # Specify the file name
#upper boundry 

# Assuming vectora_low and list_p_low are vectors
df6 <- data.frame(vectora, list_p = list_p[,1])  # Convert to dataframe

ggplot(df6, aes(x = vectora, y = list_p)) +
  geom_point(color = "blue", size = 2, alpha = 0.6) +  # Scatter plot
  geom_hline(yintercept = 0.025, linetype = "dashed", color = "red", size = 1) +  # Horizontal line
  labs(
    title = "P-values across different hypothetical constant treatment effect for cg14795253",
    x = "Hypothetical constant treatment effect (a)",
    y = "p-value"
  ) +
  theme_classic()   # Use a clean theme

dev.off()
 

# confidence interval: [0.002891579, 0.009277008]


# section 4: cg08301503 of C6orf227. effect size = 0.020336783
 

mn <- which(colnames(final_secondary_CpG)=="cg08301503")   #provide the column number for the CpG site 

#all possible randomizations
n <- 2^17

#create a data frame that contains the CpG site of interest
final1<-
  final_secondary_CpG %>%
  select(c(1:4),mn) %>%
  rename(Y_obs = cg08301503)

#data wrangling
#ozone, 1=O3
final2 <-
  final1%>%
  filter(exp==1) %>%
  rename(Yi_wiO3 = Y_obs)%>%
  select(3,4,5)

#clean air, -1=CA
final3 <-
  final1 %>%
  filter(exp==-1)%>%
  rename(Yi_wiCA = Y_obs)%>%
  select(3,4,5)


#full join table final4 with final2 or final 3, respectively
final4 <-
  full_join(final1, final2, by = c("id2" = "id2", "exp" = "exp")) %>%
  select(1,2,3,4,6,5)

final4<-
  full_join(final4, final3, by = c("id2" = "id2", "exp" = "exp")) %>%
  select(1,2,3,4,5,7,6)
 


#for lower boundry. This CpG site has a positive effect size -->calculating p value should use >=
#should have no absolute value when calculating p value
 
a=200
vectora_low = as.matrix(seq(-0.006, 0.043, length.out = a)) #if changing the number of a values here length.out=# of a values

x<-list()

for(i in 1:a) { # of a value in the loop
  #substitute the NA value in Yi_1 with the correponding value in Yi_0 but "a" more
  final4 <-
    full_join(final1, final2, by = c("id2" = "id2", "exp" = "exp")) %>%
    select(1,2,3,4,6,5)
  
  final4<-
    full_join(final4, final3, by = c("id2" = "id2", "exp" = "exp")) %>%
    select(1,2,3,4,5,7,6)
  
  #substitute the NA value in Yi_1 with the correponding value in Yi_0 but a more
  final4$Yi_wiO3[is.na(final4$Yi_wiO3)] <- final4$Yi_wiCA[final4$exp==-1]+vectora_low[i,1]
  final4$Yi_wiCA[is.na(final4$Yi_wiCA)] <- final4$Yi_wiO3[final4$exp==1]-vectora_low[i,1]
  x[[i]]<-final4
}
 


 

list_p_low <- matrix(nrow=a, ncol=1)  # nrow= number of a values

# Function for processing a single `i`
compute_p_value <- function(i) {
  id <- x[[i]]$"id2"
  W_obs <- x[[i]]$"exp"
  Yi_CA <- x[[i]]$"Yi_wiCA"
  Yi_O3 <- x[[i]]$"Yi_wiO3"
  Y_obs <- x[[i]]$"Y_obs"
  
  # Compute statistics for all j in one step
  stat <- apply(w0.true, 1, function(w) { #1 represents applying each row in W0.true
    t.test(Yi_O3[w == 1], Yi_CA[w == -1],paired=TRUE)$statistic
  })
  
  # Observed statistics (no absolute value)
  stat_obs <- t.test(Y_obs[W_obs == 1], Y_obs[W_obs == -1],paired=TRUE)$statistic
  
  # Compute p-value
  p_value <- 1*sum(stat>=stat_obs)/n
  
  return(p_value)
}

# Parallel execution
cl <- makeCluster(n_cores)
clusterExport(cl, c("x", "w0.true", "n"))  # Added "n" to the export list
list_p_low[, 1] <- parSapply(cl, 1:a, compute_p_value)
stopCluster(cl)

 

# for upper boundry. This CpG site has a positive effect size -->calculating p value should use <=
 
a=200
vectora = as.matrix(seq(-0.004, 0.047368421,length.out = a)) #if changing the number of a values here length.out=# of a values

x<-list()

for(i in 1:a) { # of a value in the loop
  #substitute the NA value in Yi_1 with the correponding value in Yi_0 but "a" more
  final4 <-
    full_join(final1, final2, by = c("id2" = "id2", "exp" = "exp")) %>%
    select(1,2,3,4,6,5)
  
  final4<-
    full_join(final4, final3, by = c("id2" = "id2", "exp" = "exp")) %>%
    select(1,2,3,4,5,7,6)
  
  #substitute the NA value in Yi_1 with the correponding value in Yi_0 but a more
  final4$Yi_wiO3[is.na(final4$Yi_wiO3)] <- final4$Yi_wiCA[final4$exp==-1]+vectora[i,1]
  final4$Yi_wiCA[is.na(final4$Yi_wiCA)] <- final4$Yi_wiO3[final4$exp==1]-vectora[i,1]
  x[[i]]<-final4
}
 


 
#create a matrix with all p-values corresponding to each value a
list_p<- matrix(nrow=a, ncol=1)  # nrow= number of a values

# Function for processing a single `i`
compute_p_value <- function(i) {
  id <- x[[i]]$"id2"
  W_obs <- x[[i]]$"exp"
  Yi_CA <- x[[i]]$"Yi_wiCA"
  Yi_O3 <- x[[i]]$"Yi_wiO3"
  Y_obs <- x[[i]]$"Y_obs"
  
  # Compute statistics for all j in one step
  stat <- apply(w0.true, 1, function(w) { #1 represents applying each row in W0.true
    t.test(Yi_O3[w == 1], Yi_CA[w == -1], paired=TRUE)$statistic
  })
  
  # Observed statistics (no absolute value)
  stat_obs <- t.test(Y_obs[W_obs == 1], Y_obs[W_obs == -1], paired=TRUE)$statistic
  
  # Compute p-value
  p_value <- 1*sum(stat<=stat_obs)/n
  
  return(p_value)
}

# Parallel execution
cl <- makeCluster(n_cores)
clusterExport(cl, c("x", "w0.true", "n"))  # Added "n" to the export list
list_p[, 1] <- parSapply(cl, 1:a, compute_p_value)
stopCluster(cl)
 


# visualization 
# graphing of lower boundry (df7)
 
pdf(file.path(graph_out_dir, "cg08301503_fish_int_low.pdf"), width = 8, height = 6, pointsize = 12, family = "Times")  # Specify the file name

#lower boundry 

# Assuming vectora_low and list_p_low are vectors
df7 <- data.frame(vectora_low, list_p_low = list_p_low[,1])  # Convert to dataframe

ggplot(df7, aes(x = vectora_low, y = list_p_low)) +
  geom_point(color = "blue", size = 2, alpha = 0.6) +  # Scatter plot
  geom_hline(yintercept = 0.025, linetype = "dashed", color = "red", size = 1) +  # Horizontal line
  labs(
    title = "P-values across different hypothetical constant treatment effect for cg08301503",
    x = "Hypothetical constant treatment effect (a)",
    y = "p-value"
  ) +
  theme_classic()   # Use a clean theme

dev.off()

 

#graphing of upper boundry (df8)
 
pdf(file.path(graph_out_dir, "cg08301503_fish_int_up.pdf"), width = 8, height = 6, pointsize = 12, family = "Times")  # Specify the file name
#upper boundry 

# Assuming vectora_low and list_p_low are vectors
df8 <- data.frame(vectora, list_p = list_p[,1])  # Convert to dataframe

ggplot(df8, aes(x = vectora, y = list_p)) +
  geom_point(color = "blue", size = 2, alpha = 0.6) +  # Scatter plot
  geom_hline(yintercept = 0.025, linetype = "dashed", color = "red", size = 1) +  # Horizontal line
  labs(
    title = "P-values across different hypothetical constant treatment effect for cg08301503",
    x = "Hypothetical constant treatment effect (a)",
    y = "p-value"
  ) +
  theme_classic()   # Use a clean theme

dev.off()
 


# confidence interval: [0.007710526, 0.03286842]


# section 5: cg01392313 of C6orf227. effect size = 0.028433911
 

mn <- which(colnames(final_secondary_CpG)=="cg01392313")   #provide the column number for the CpG site 

#all possible randomizations
n <- 2^17

#create a data frame that contains the CpG site of interest
final1<-
  final_secondary_CpG %>%
  select(c(1:4),mn) %>%
  rename(Y_obs = cg01392313)

#data wrangling
#ozone, 1=O3
final2 <-
  final1%>%
  filter(exp==1) %>%
  rename(Yi_wiO3 = Y_obs)%>%
  select(3,4,5)

#clean air, -1=CA
final3 <-
  final1 %>%
  filter(exp==-1)%>%
  rename(Yi_wiCA = Y_obs)%>%
  select(3,4,5)


#full join table final4 with final2 or final 3, respectively
final4 <-
  full_join(final1, final2, by = c("id2" = "id2", "exp" = "exp")) %>%
  select(1,2,3,4,6,5)

final4<-
  full_join(final4, final3, by = c("id2" = "id2", "exp" = "exp")) %>%
  select(1,2,3,4,5,7,6)
 


#for lower boundry. This CpG site has a positive effect size -->calculating p value should use >=
#should have no absolute value when calculating p value
 
a=200
vectora_low = as.matrix(seq(-0.012, 0.06, length.out = a)) #if changing the number of a values here length.out=# of a values

x<-list()

for(i in 1:a) { # of a value in the loop
  #substitute the NA value in Yi_1 with the correponding value in Yi_0 but "a" more
  final4 <-
    full_join(final1, final2, by = c("id2" = "id2", "exp" = "exp")) %>%
    select(1,2,3,4,6,5)
  
  final4<-
    full_join(final4, final3, by = c("id2" = "id2", "exp" = "exp")) %>%
    select(1,2,3,4,5,7,6)
  
  #substitute the NA value in Yi_1 with the correponding value in Yi_0 but a more
  final4$Yi_wiO3[is.na(final4$Yi_wiO3)] <- final4$Yi_wiCA[final4$exp==-1]+vectora_low[i,1]
  final4$Yi_wiCA[is.na(final4$Yi_wiCA)] <- final4$Yi_wiO3[final4$exp==1]-vectora_low[i,1]
  x[[i]]<-final4
}
 


 

list_p_low <- matrix(nrow=a, ncol=1)  # nrow= number of a values

# Function for processing a single `i`
compute_p_value <- function(i) {
  id <- x[[i]]$"id2"
  W_obs <- x[[i]]$"exp"
  Yi_CA <- x[[i]]$"Yi_wiCA"
  Yi_O3 <- x[[i]]$"Yi_wiO3"
  Y_obs <- x[[i]]$"Y_obs"
  
  # Compute statistics for all j in one step
  stat <- apply(w0.true, 1, function(w) { #1 represents applying each row in W0.true
    t.test(Yi_O3[w == 1], Yi_CA[w == -1],paired=TRUE)$statistic
  })
  
  # Observed statistics (no absolute value)
  stat_obs <- t.test(Y_obs[W_obs == 1], Y_obs[W_obs == -1],paired=TRUE)$statistic
  
  # Compute p-value
  p_value <- 1*sum(stat>=stat_obs)/n
  
  return(p_value)
}

# Parallel execution
cl <- makeCluster(n_cores)
clusterExport(cl, c("x", "w0.true", "n"))  # Added "n" to the export list
list_p_low[, 1] <- parSapply(cl, 1:a, compute_p_value)
stopCluster(cl)

 

# for upper boundry. This CpG site has a positive effect size -->calculating p value should use <=
 
a=200
vectora = as.matrix(seq(-0.007, 0.060526316, length.out = a)) #if changing the number of a values here length.out=# of a values

x<-list()

for(i in 1:a) { # of a value in the loop
  #substitute the NA value in Yi_1 with the correponding value in Yi_0 but "a" more
  final4 <-
    full_join(final1, final2, by = c("id2" = "id2", "exp" = "exp")) %>%
    select(1,2,3,4,6,5)
  
  final4<-
    full_join(final4, final3, by = c("id2" = "id2", "exp" = "exp")) %>%
    select(1,2,3,4,5,7,6)
  
  #substitute the NA value in Yi_1 with the correponding value in Yi_0 but a more
  final4$Yi_wiO3[is.na(final4$Yi_wiO3)] <- final4$Yi_wiCA[final4$exp==-1]+vectora[i,1]
  final4$Yi_wiCA[is.na(final4$Yi_wiCA)] <- final4$Yi_wiO3[final4$exp==1]-vectora[i,1]
  x[[i]]<-final4
}
 


 
#create a matrix with all p-values corresponding to each value a
list_p<- matrix(nrow=a, ncol=1)  # nrow= number of a values

# Function for processing a single `i`
compute_p_value <- function(i) {
  id <- x[[i]]$"id2"
  W_obs <- x[[i]]$"exp"
  Yi_CA <- x[[i]]$"Yi_wiCA"
  Yi_O3 <- x[[i]]$"Yi_wiO3"
  Y_obs <- x[[i]]$"Y_obs"
  
  # Compute statistics for all j in one step
  stat <- apply(w0.true, 1, function(w) { #1 represents applying each row in W0.true
    t.test(Yi_O3[w == 1], Yi_CA[w == -1], paired=TRUE)$statistic
  })
  
  # Observed statistics (no absolute value)
  stat_obs <- t.test(Y_obs[W_obs == 1], Y_obs[W_obs == -1], paired=TRUE)$statistic
  
  # Compute p-value
  p_value <- 1*sum(stat<=stat_obs)/n
  
  return(p_value)
}

# Parallel execution
cl <- makeCluster(n_cores)
clusterExport(cl, c("x", "w0.true", "n"))  # Added "n" to the export list
list_p[, 1] <- parSapply(cl, 1:a, compute_p_value)
stopCluster(cl)
 

# confidence interval: [0.01205263, 0.04490304]

# visualization 
# graphing of lower boundry (df9)
 
pdf(file.path(graph_out_dir, "cg01392313_fish_int_low.pdf"), width = 8, height = 6, pointsize = 12, family = "Times")  # Specify the file name

#lower boundry 

# Assuming vectora_low and list_p_low are vectors
df9 <- data.frame(vectora_low, list_p_low = list_p_low[,1])  # Convert to dataframe

ggplot(df9, aes(x = vectora_low, y = list_p_low)) +
  geom_point(color = "blue", size = 2, alpha = 0.6) +  # Scatter plot
  geom_hline(yintercept = 0.025, linetype = "dashed", color = "red", size = 1) +  # Horizontal line
  labs(
    title = "P-values across different hypothetical constant treatment effect for cg01392313",
    x = "Hypothetical constant treatment effect (a)",
    y = "p-value"
  ) +
  theme_classic()   # Use a clean theme

dev.off()

 

#graphing of upper boundry (df10)
 
pdf(file.path(graph_out_dir, "cg01392313_fish_int_up.pdf"), width = 8, height = 6, pointsize = 12, family = "Times")  # Specify the file name
#upper boundry 

# Assuming vectora_low and list_p_low are vectors
df10 <- data.frame(vectora, list_p = list_p[,1])  # Convert to dataframe

ggplot(df10, aes(x = vectora, y = list_p)) +
  geom_point(color = "blue", size = 2, alpha = 0.6) +  # Scatter plot
  geom_hline(yintercept = 0.025, linetype = "dashed", color = "red", size = 1) +  # Horizontal line
  labs(
    title = "P-values across different hypothetical constant treatment effect for cg01392313",
    x = "Hypothetical constant treatment effect (a)",
    y = "p-value"
  ) +
  theme_classic()   # Use a clean theme

dev.off()
 
