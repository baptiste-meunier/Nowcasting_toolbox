# Code ranking potential regressors by informative power with respects to a target variable
# 
# USER INPUTS: 
#     - an input dataset which should:
#           o be organized along the lines of the 'data_xxx.xlsx' files (xxx = country name) in the folder 'dataset'
#           o be named 'data_xxx.xlsx' (xxx = country name) and be located in the folder 'dataset'
#           o have both quarterly and monthly variables (resp. in 'Monthly' and 'Quarterly' tabs)
#           o NB: in this code, the name of input dataset is automatically inferred from the country name
#     - the mnemotechnic for the target variable (generally the Haver or FAME handle)
#     - a choice of pre-selection method. Three are available: 
#           o 1 = SIS: ranking of regressors is based on the correlation of the regressor with the target variable
#           o 2 = LARS: ranking of regressors is based is based on the LARS algorithm
#           o 3 = t-stat based: ranking of regressors is based on the t-stat of its coefficient in a univariate regression on the target variable
#           o NB: The literature generally finds more accurate forecasts if LARS is used. 
#                 Among the three, this is the only multivariate method (i.e. other variables are taken into account when assessing the informative power of a regressor). 
#                 Other methods are univariate: using them might result in selecting very similar variables, therefore with only limited informative power after the 2nd or 3rd variable. 
#     - the period (start and end dates) on which the pre-selection is performed
#     - the numnber of leads or lags to apply to the target series
#           o lead = 0 for a pre-selection on contemporaneous dates
#           o lead > 0 for a pre-selection with the target series taken at n-period ahead (and at n-period back if lead < 0)
#
# OUTPUTS:
#     - A CSV file with the ranking of the regressors (from most to least informative). NB:
#           o The CSV file also includes correlation coefficient (for SIS) and t-stat values (for t-stat based)
#           o The name of the CSV file specifies the country, the pre-selection method used (1 to 3), and the period on which it is conducted
#           o The output file is located in the folder 'eval' in the sub-folder corresponding to the country
#           o The output file also contains information about the frequency, group, and publication lags of all series
#               > NB1: publication delay is relative to the timeliest series (for which delay=0). Other delays are in comparison with this timeliest series.
#               > NB2: publication delay is automatically infered by the code givne the presence of NA. Make sure that all series in the input have been updated at the same time to get meaningful inference.
#
# SOURCE:
# Code taken from Chinn, M. D., Meunier, B., Stumpner, S. (2023).
#                 "Nowcasting World Trade with Machine Learning: a Three-Step Approach".
#                 NBER Working Paper, No 31419, National Bureau of Economic Research
#
# When using techniques, users are kindly requested to cite the original papers:
#   - Least Angle Regression: Efron, B., Hastie, T., Johnstone, I., and Tibshirani, R. (2004). 
#                             “Least angle regression”, 
#                             Annals of Statistics, 32(2), pp. 407–499
#   - Sure Independence Screening: Fan, J., and Lv, J. (2008). 
#                                  “Sure independence screening for ultrahigh dimensional feature space”, 
#                                  Journal of the Royal Statistical Society Series B, 70(5), pp. 849–911
#   - T-stat-based: Bair, E., Hastie, T., Paul, D., and Tibshirani, R. (2006). 
#                   “Prediction by supervised principal components”, 
#                   Journal of the American Statistical Association, 101(473), pp. 119–137
#

rm(list = ls())

#
# 0. SET USER PARAMETERS
#

# User parameters
country <- "Example1"      # country 
                           # NB: it should match the name for the Excel input file (in folder 'dataset')
                           #     it should also  have a sub-folder named after it in folder 'eval'
target <- "Name target"    # name of the target variable (line 4 in input Excel)
start_date <- "2004-03-31" # start date for the pre-selection. 
                           # NB: only quarters starting after this date will be retained
end_date <- "2023-01-01"   # end date for the pre-selection
                           # NB: only quarters ending before this date will be retained
lead <- 0                  # Number of leads for the target variable
                           # 0 = contemporaneous
                           # >0 = for n-period ahead
                           # <0 = for n-period back
select_method <- 3         # choice of pre-selection method
                           # 1 = correlation-based (Sure Independence Screening of Fan and Lv, 2008)
                           # 2 = LARS (Efron et al., 2004)
                           # 3 = t-stat based (Bair et al., 2006)

# If needed for installing packages, uncomment line below (and update the package name)
# install.packages("readr")

# Importing libraries
library(magrittr)
library(dplyr)
library(tidyr)
library(lubridate)
library(tidyverse)
library(readr)
library(readxl)
library(reshape2)
library(zoo)
library(stringr)
library(lars)


#
# 1. IMPORT DATA
#

# Import monthly data in levels
data_mth_init <- read_excel(paste0("./dataset/data_",country,".xlsx"),
                            sheet="Monthly",
                            col_names = FALSE) %>%
  select(-'...1')

transfo_mth <- head(data_mth_init,1)[-1] # NA in first cell because are dates
groups_mth <- data_mth_init[2,][-1] # NA in first cell because are dates

data_mth <- data_mth_init
colnames(data_mth) <- filter(data_mth_init,row_number()==3)
colnames(transfo_mth) <- colnames(data_mth)[-1]
colnames(groups_mth) <- colnames(data_mth)[-1]

if(sum(duplicated(colnames(data_mth)))>0){ # if there are duplicates, the procedure will not work
  xxx <- colnames(data_mth)
  stop(paste0("Monthly input data have duplicates, please delete the variables that appear twice: ", str_c(xxx[duplicated(xxx)],collapse=', ')))
}

data_mth %<>%
  filter(row_number()>=11) %>%
  mutate(across(everything(),
                ~as.numeric(.))) %>%
  rename('date'='.excel_last') %>%
  mutate(date = as_date(as.Date(as.numeric(date), origin="1899-12-30")))

# Get publication delays for monthly data 
# NB: publication delay is computed relative to the series with the later date (set to 0)
#     it means that a publication delay of 1 does NOT mean that the variable is necessarily
#     published with a one-month delay - just that this is published one month after the most
#     timely series of the dataset
temp <- data_mth %>%
  filter(!is.na(date))
last_datem <- temp$date[nrow(temp)]
temp %<>%
  select(-date)
n_obs <- nrow(temp)
nlagsm <- data.frame(matrix(NA,
                            nrow = ncol(temp),
                            ncol = 3))
for (jj in 1:ncol(temp)){
  
  # Select individual variable
  temp_col <- temp[jj]
  temp_name <- colnames(temp_col)
  
  # Check last non-NA position and difference with n_obs
  non_NA_index <- which(!is.na(temp_col))
  last_non_NA <- max(non_NA_index)
  diff_na <- n_obs - last_non_NA 
  
  # Write in nlags
  nlagsm[jj,1] <- temp_name
  nlagsm[jj,2] <- diff_na
  nlagsm[jj,3] <- 'M'
  
}

# Import quarterly data in levels
data_qtr_init <- read_excel(paste0("./dataset/data_",country,".xlsx"),
                            sheet="Quarterly",
                            col_names = FALSE) %>%
  select(-'...1')

transfo_qtr <- head(data_qtr_init,1)[-1] # NA in first cell because are dates
groups_qtr <- data_qtr_init[2,][-1] # NA in first cell because are dates

data_qtr <- data_qtr_init
colnames(data_qtr) <- filter(data_qtr_init,row_number()==3)
colnames(transfo_qtr) <- colnames(data_qtr)[-1]
colnames(groups_qtr) <- colnames(data_qtr)[-1]

if(sum(duplicated(colnames(data_qtr)))>0){ # if there are duplicates, the procedure will not work
  xxx <- colnames(data_qtr)
  stop(paste0("Quarterly input data have duplicates, please delete the variables that appear twice: ", str_c(xxx[duplicated(xxx)],collapse=', ')))
}

data_qtr %<>%
  filter(row_number()>=11) %>%
  mutate(across(everything(),
                ~as.numeric(.))) %>%
  rename('date'='.excel_last') %>%
  mutate(date = as_date(as.Date(as.numeric(date), origin="1899-12-30")))

# Get publication delays for quarterly data (relative to most advanced date = 0)
temp <- data_qtr %>%
  filter(!is.na(date))
last_dateq <- temp$date[nrow(temp)]
temp %<>%
  select(-date, -all_of(target))
n_obs <- nrow(temp)
nlagsq <- data.frame(matrix(NA,
                            nrow = ncol(temp),
                            ncol = 3))
if(ncol(temp)>0){
  for (jj in 1:ncol(temp)){
    
    # Select individual variable
    temp_col <- temp[jj]
    temp_name <- colnames(temp_col)
    
    # Check last non-NA position and difference with n_obs
    non_NA_index <- which(!is.na(temp_col))
    last_non_NA <- max(non_NA_index)
    diff_na <- n_obs - last_non_NA 
    
    # Write in nlags
    nlagsq[jj,1] <- temp_name
    nlagsq[jj,2] <- diff_na
    nlagsq[jj,3] <- 'Q'
    
  }
}

# Merge publication lags for monthly and quarterly variables
diff_qtom <- interval(ymd(last_dateq),ymd(last_datem)) %/% months(1) # Gives the number of months between the two
nlags_all <- rbind(nlagsq,nlagsm)
if(diff_qtom>=0){
  nlags_all %<>%
    mutate(X2=ifelse(X3=='Q',3*X2+diff_qtom,X2))
}else{
  nlags_all %<>%
    mutate(X2=ifelse(X3=='Q',3*X2,X2-diff_qtom))
}
colnames(nlags_all) <- c('variable','publication delay (relative to most timely = 0)','frequency')

# Merge groups names for monthly and quarterly
groups_all <- cbind(groups_mth,groups_qtr) %>%
  t() %>%
  as.data.frame() %>%
  rename(group=V1) %>%
  rownames_to_column("variable")

# Merge with nlags_all
nlags_final <- nlags_all %>%
  merge(groups_all,by="variable",all.x=TRUE)


#
# 2. TRANSFORM DATA
#

# Get transformations (monthly)
# Growth rates are actually dlog - as is the case in the DFM code
transfo_mth %<>% t() %>%
  as.data.frame()
list_level <- rownames(filter(transfo_mth,V1==0)) # 0 = level
list_mom <- rownames(filter(transfo_mth,V1==1)) # 1 = month-on-month growth rate
list_diff <- rownames(filter(transfo_mth,V1==2)) # 2 = first difference
list_qoqAR <- rownames(filter(transfo_mth,V1==3)) # 3 = annualized quarter-on-quarter growth rate
list_yoy <- rownames(filter(transfo_mth,V1==4)) # 4 = year-on-year growth rate

# Run transformations (monthly)
data_mth_transfo <- data_mth %>%
  mutate(across(all_of(list_mom),
                ~(log(.)-log(lag(.))))) %>%
  mutate(across(all_of(list_diff),
                ~(.-lag(.)-1))) %>%
  mutate(across(all_of(list_qoqAR),
                ~(log(.)-log(lag(.,3)))*4)) %>%
  mutate(across(all_of(list_yoy),
                ~(log(.)-log(lag(.,12)))))

# Standardize and match into qoq using the Mariano and Murasawa (2003) approximation
# The same approximation is used in the DFM (see R matrix in DFM code)
data_mth_transfo %<>%
  mutate(across(-date,
                ~scale(.))) %>%
  mutate(across(-date,
                ~1*.+2*lag(.)+3*lag(.,2)+2*lag(.,3)+1*lag(.,4)))

# Get transformations (quarterly)
transfo_qtr %<>% t() %>%
  as.data.frame()
list_level <- rownames(filter(transfo_qtr,V1==0)) # 0 = level
list_qoq <- rownames(filter(transfo_qtr,V1==1)) # 1 = quarter-on-quarter growth rate
list_diff <- rownames(filter(transfo_qtr,V1==2)) # 2 = first difference
list_qoqAR <- rownames(filter(transfo_qtr,V1==3)) # 3 = annualized quarter-on-quarter growth rate
list_yoy <- rownames(filter(transfo_qtr,V1==4)) # 4 = year-on-year growth rate

# Run transformations (quarterly)
data_qtr_transfo <- data_qtr %>%
  mutate(across(all_of(list_qoq),
                ~(log(.)-log(lag(.))))) %>%
  mutate(across(all_of(list_diff),
                ~(.-lag(.)-1))) %>%
  mutate(across(all_of(list_qoqAR),
                ~(log(.)-log(lag(.)))*4)) %>%
  mutate(across(all_of(list_yoy),
                ~(log(.)-log(lag(.,4)))))

# Merge datasets
data_all_transfo <- data_qtr_transfo %>%
  merge(data_mth_transfo,by="date",all.x=TRUE)
  

#
# 3. PREPARE DATA FOR PRE-SELECTION
#

# Prepare target series
if(lead>=0){
  target_data <- data_all_transfo %>%
    select(date, target = all_of(target)) %>%
    mutate(target = lead(target,lead)) %>%
    filter(date > as.Date(start_date)) %>%
    filter(date < as.Date(end_date))
}else{
  target_data <- data_all_transfo %>%
    select(date, target = all_of(target)) %>%
    mutate(target = lag(target,lead)) %>%
    filter(date > as.Date(start_date)) %>%
    filter(date < as.Date(end_date))
}   

# Prepare LHS variables
lhs_data <- data_all_transfo %>%
  select(-all_of(target)) %>%
  filter(date > as.Date(start_date)) %>%
  filter(date < as.Date(end_date)) %>%
  select(-date)

# Remove series with NAs
# NB: It is necessary to run only on LHS (otherwise, if target data has NA it will be discarded)
lhs_data <- lhs_data[ , colSums(is.na(lhs_data)) == 0]
data_final <- target_data %>%
  cbind(lhs_data)
  
# Print list of series excluded
full_series <- colnames(select(data_all_transfo,-all_of(target)))
selected_series <- colnames(data_final)
excluded_series <- setdiff(full_series,selected_series)
print("=================================================================================")
print("                         LIST OF VARIABLES EXCLUDED                              ")
print(excluded_series)
print("=================================================================================")


#
# 4. RUN PRE-SELECTION
#

# Initiate data
smpl <- data_final %>%
  select(-date) %>%
  drop_na()

if(select_method==1){ # 1 = Correlation-based (SIS: Fan and Lv, 2008)
  
  order <- smpl %>%
    mutate(across(-target,~cor(.,
                               target,
                               use = "pairwise.complete.obs",
                               method  = "pearson"))) %>%
    select(-target) %>%
    distinct() %>%
    t() %>%
    as.data.frame() %>%
    rename(corr=V1)

  order %<>%
    mutate(corr = abs(corr)) %>%
    arrange(desc(corr)) %>%
    tibble::rownames_to_column("variable")
  
}else if(select_method==2){ # 2 = LARS (Efron et al., 2004)
  
  # Creating batches for LARS loop
  n_var = ncol(smpl)-1
  batch <- smpl
  count <- 0
  order <- data.frame(V1=c(0))
  
  # LARS loop (have to run in different batchs if number of variables is too high)
  while (count < n_var) {
    
    # Running LARS equation
    x <- as.matrix(select(batch,-target))
    y <- as.matrix(select(batch,target))
    eq <- lars(x = x,
               y = y,
               type="lar")
    
    # Ordering the variables
    out <- as.data.frame(coef(eq)) %>%
      summarise(across(everything(),~sum(.==0))) %>%
      t() %>%
      as.data.frame() %>%
      tibble::rownames_to_column('name') %>%
      arrange(V1) %>%
      group_by(V1) %>%
      filter(n()==1) %>%
      tibble::column_to_rownames('name')
    
    order <- rbind(order,out)
    
    # Deleting variables already ordered from the sample
    var <- row.names(out)
    batch %<>%
      select(!all_of(var))
    
    # Checking if nrow(out) = 0 to avoid infinite loops
    if (nrow(out)==0){
      # Putting all remaining variables at the end
      x <- as.matrix(select(batch,-target))
      n_end <- ncol(x)
      out <- data.frame(V1=rep(1,n_end))
      rownames(out) <- colnames(x)
      order <- rbind(order,out)
      count <- count + n_end
      print("Warning: Using special procedure for LARS")
    }else{
      # Just updating the count      
      count <- count + nrow(out)
    } 
  }
  
  order %<>%
    filter(row_number()>1) %>%
    rename("pecking order"=V1) %>%
    tibble::rownames_to_column("variable")
  

}else if(select_method==3){ # 3 = t-stat based (Bair et al., 2006)
  
  # Create 2 lags of the dependent variable as in Bair et al. (2006)
  init <- smpl %>%
    mutate(L1_target = lag(target,1),
           L2_target = lag(target,2)) %>%
    drop_na()
  
  list_var <- colnames(select(smpl,-ends_with("target")))
  order <- data.frame(variable=NA,
                      tstat=NA)
  
  for (v in seq_along(list_var)){
    
    data_eq <- init %>%
      select(ends_with("target"),
             list_var[v])
    
    eq <- lm(target ~ .,
             data = data_eq,
             na.action = "na.omit")
    
    order[v,1] <- list_var[v]
    order[v,2] <- abs(coef(summary(eq))[4, "t value"])
  }
  
  # Order by t-stat (higher to lower)
  order %<>%
    arrange(desc(tstat))

}


#
# 5. WRITE OUTPUT
#

# Add information on publication delays, frequency, and groups
final_output <- order %>%
  left_join(nlags_final,by="variable")

# Ranking of variables
write.csv(final_output,
          file=paste0("./eval/",country,"/",country,"_preselection_method",select_method,"_lead",lead,"_",start_date,"_to_",end_date,".csv"))
