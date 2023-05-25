library(readxl)
library(tidyverse)
library(gdata)
library(car)

data_raw <- read_excel("data/arabidopsis_data.xlsx", col_types = c("text", rep("numeric", 8), "text")) #making sure every value is numeric

metadata <- tibble("Group" = c("L1", "L2", "S3", "S4", "L5", "L6", "S7", "S8"),
       "Population" = rep(c("Bovra", "Smadalen", "Helin", "Spiterstulen"), each = 2),
       "Size" = rep(c("Large", "Small", "Large", "Small"), each = 2),
       "Elevation" = rep("High", each = 8),
       "Condition" = rep(c("Stable", "Declining", "Declining", "Declining"), each = 2), #difference between "decline" and "decrease"? https://brightspace.ru.nl/d2l/le/content/345973/viewContent/1977818/View
       "Soil" = rep(c("Riverbed", "Riverbed", "Scree", "Rocks"), each = 2),
       "Conc_NaCl" = c(0, 50, 0, 50, 0, 50, 0, 50),
       "Conc_unit" = rep("mM", each = 8))

data <- mutate(data_raw,
               "Length_dif" = data_raw$Length_longest_leaf_a - data_raw$Length_longest_leaf_b,
               "Number_dif" = data_raw$Number_leaves_a - data_raw$Number_leaves_b) #calculating differences between before and after treatment

data_complete <- left_join(data, metadata, by = "Group") #create one complete dataframe

names <- c("Length_longest_leaf_b", "Number_leaves_b", "Number_leaves_a", "Length_longest_leaf_a",  "Wet_weight", "Length_dif", "Number_dif") #collect the names of the columns of interest

sw_results <- vector("list", 7) #prepare an empty list

for (i in 1:7) { #create a list of dataframes where every dataframe contains its Shapiro-Wilk p-values
result <- data_complete %>% group_by(Size) %>% summarize(sw = shapiro.test(!!sym(names[i]))$p.value) %>% rename.vars(to = names[i], from ="sw")
 sw_results[[i]] <- result
}

merged <- sw_results %>% reduce(full_join, by= "Size") #merge all the dataframes from the list

tidy_sw_results <- pivot_longer(data = merged, cols = names,  
                             names_to = "Info",  values_to = "Shapiro-Wilk_p.value") #make the dataframe tidy

tidy_sw_results %>% filter(`Shapiro-Wilk_p.value` > 0.05) #filter for which groups are normally distributed


