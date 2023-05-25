library(readxl)
library(tidyverse)

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
               "Length_dif" = data$Length_longest_leaf_a - data$Length_longest_leaf_b,
               "Number_dif" = data$Number_leaves_a - data$Number_leaves_b) #calculating differences between before and after treatment

#names <- c("Length_longest_leaf_b", "Number_leaves_b", "Number_leaves_a", "Length_longest_leaf_a", "Leaf_shape",  "Wet_weight", "Length_dif", "Number_dif")

#for (i in length(names)) {
# sw_results[[i]] <- data %>% group_by(Group) %>% summarize(p.value_sw = shapiro.test(!!sym(names[i]))$p.value)
#}

#Finding which groups aren't normally distributed
data %>% group_by(Group) %>% summarize(p.value_sw = shapiro.test(Wet_weight)$p.value) %>% filter(p.value_sw < 0.05)
data %>% group_by(Group) %>% summarize(p.value_sw = shapiro.test(Length_dif)$p.value) %>% filter(p.value_sw < 0.05)
data %>% group_by(Group) %>% summarize(p.value_sw = shapiro.test(Number_dif)$p.value) %>% filter(p.value_sw < 0.05)
