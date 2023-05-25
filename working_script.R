library(readxl)
library(tidyverse)
library(gdata)
library(car)
library(cowplot)
library(magick)
library(here)

data_raw <- read_excel(here("data/arabidopsis_data.xlsx"), col_types = c("text", rep("numeric", 8), "text")) #making sure every value is numeric

metadata <- tibble("Group" = c("L1", "L2", "S3", "S4", "L5", "L6", "S7", "S8"),
       "Population" = rep(c("Bovra", "Smadalen", "Helin", "Spiterstulen"), each = 2),
       "Size" = rep(c("Large", "Small", "Large", "Small"), each = 2),
       "Elevation" = rep("High", each = 8),
       "Condition" = rep(c("Stable", "Declining", "Declining", "Declining"), each = 2), #difference between "decline" and "decrease"? https://brightspace.ru.nl/d2l/le/content/345973/viewContent/1977818/View
       "Soil" = rep(c("Riverbed", "Riverbed", "Scree", "Rocks"), each = 2),
       "Treatment" = c("Unsalted", "Salted", "Unsalted", "Salted", "Unsalted", "Salted", "Unsalted", "Salted"),
       "Conc_NaCl" = c(0, 50, 0, 50, 0, 50, 0, 50),
       "Conc_unit" = rep("mM", each = 8))

data <- mutate(data_raw,
               "Length_dif" = data_raw$Length_longest_leaf_a - data_raw$Length_longest_leaf_b,
               "Number_dif" = data_raw$Number_leaves_a - data_raw$Number_leaves_b) #calculating differences between before and after treatment

data_complete <- left_join(data, metadata, by = "Group") #create one complete dataframe

names <- c("Wet_weight", "Length_dif", "Number_dif") #collect the names of the columns of interest

sw_results <- vector("list", 3) #prepare an empty list

for (i in 1:3) { #create a list of dataframes where every dataframe contains its Shapiro-Wilk p-values
result <- data_complete %>% group_by(Size, Treatment) %>% summarize(sw = shapiro.test(!!sym(names[i]))$p.value) %>% rename.vars(to = names[i], from ="sw")
result <- unite(result, "Treatment", Size, Treatment)
 sw_results[[i]] <- result
}

merged <- sw_results %>% reduce(full_join, by= "Treatment") #merge all the dataframes from the list

tidy_sw_results <- pivot_longer(data = merged, cols = names,  
                             names_to = "Info",  values_to = "ShapiroWilk_p.value") #make the dataframe tidy

tidy_sw_results %>% filter(ShapiroWilk_p.value > 0.05) #filter for which groups are normally distributed
tidy_sw_results %>% filter(ShapiroWilk_p.value < 0.05) #filter for which groups are NOT normally distributed

#All groups of interest are normally distributed, except for the difference in number of leaves

leveneTest(data_complete$Wet_weight ~ as.factor(Treatment), data = data_complete) #Pr(>F) = 0.4168, which means equal variance
leveneTest(data_complete$Length_dif ~ as.factor(Treatment), data = data_complete) #Pr(>F) = 0.4793, which means equal variance

t.test(formula = data_complete$Wet_weight ~ data_complete$Treatment, 
       paired = FALSE, var.equal = TRUE)$p.value %>% round(.,3) #t-test, p = 0.011, which means there is a statistical significant difference between the salted and unsalted plants.

summary <- data_complete %>%
  group_by(Treatment) %>%
  summarize(mean_noleaves.increase = mean(Number_dif, na.rm = TRUE),
            mean_lengthleaf.increase = mean(Length_dif, na.rm = TRUE),
            mean_wetweight = mean(Wet_weight, na.rm = TRUE),
            stdevNumber = sd(Number_dif, na.rm = TRUE))

data_salted <- subset(data_complete, Treatment == "Salted")
leveneTest(data_salted$Wet_weight ~ as.factor(Size), data = data_salted) #Pr(>F) = 0.3754 which means equal variance

t.test(formula = data_salted$Wet_weight ~ data_salted$Size, 
       paired = FALSE, var.equal = TRUE)$p.value %>% round(.,3) #t-test, p = 0.848, which means there is no statistical significant difference in wet weight between the large and small salted plants.

data_salted <- subset(data_complete, Treatment == "Salted")
leveneTest(data_salted$Length_dif ~ as.factor(Size), data = data_salted) #Pr(>F) = 0.7135 which means equal variance

t.test(formula = data_salted$Length_dif ~ data_salted$Size, 
       paired = FALSE, var.equal = TRUE)$p.value %>% round(.,3) #t-test, p = 0.563, which means there is no statistical significant difference between the difference in length of the longest leaf of the large and small salted plants.

wilcox.test(Number_dif ~ Size, data = data_salted) #Wilcoxon test, p = 0.4735, which means there is no statistical significant difference between the difference in number of leaves of the large and smal salted plants.

summary_salted <- data_salted %>%
  group_by(Size) %>%
  summarize(mean_noleaves.increase = mean(Number_dif, na.rm = TRUE),
            mean_lengthleaves.increase = mean(Length_dif, na.rm = TRUE),
            mean_wetweight = mean(Wet_weight, na.rm = TRUE),
            stdevlengthleave.increase = sd(Length_dif, na.rm = TRUE),
            stdevwetweight = sd(Wet_weight, na.rm = TRUE),
            stdevNumber = sd(Number_dif, na.rm = TRUE))

plot_length <- summary_salted %>% #plot the increase in longest leaf length grouped by population size
  ggplot(aes(x = Size, y = mean_lengthleaves.increase))+
  geom_col(aes(fill = Size))+
  geom_errorbar(aes(ymin = mean_lengthleaves.increase - stdevlengthleave.increase,
                    ymax = mean_lengthleaves.increase + stdevlengthleave.increase), width=.2)+
  theme_light()+
  labs(title = "Increase in number of leaves,\ngrouped by population size",
       subtitle = "Error bars depict 1 standard deviation",
       x="Population size",
       y="Mean increase in the length\nof the longest leaf in cm")+
  theme(legend.position = "none", text = element_text(size=12), plot.title = element_text(size=15),
        plot.subtitle = element_text(size=12))+
  scale_fill_manual(values=c("#e457b5", "#57e486"))

plot_noleaves <- summary_salted %>% #plot the increase in number of leaves grouped by population size
  ggplot(aes(x = Size, y = mean_noleaves.increase))+
  geom_col(aes(fill = Size))+
  geom_errorbar(aes(ymin = mean_noleaves.increase - stdevNumber, ymax = mean_noleaves.increase + stdevNumber), width=.2)+
  theme_light()+
  labs(title = "Increase in number of leaves,\ngrouped by population size",
       subtitle = "Error bars depict 1 standard deviation",
       x="Population size",
       y="Mean increase in the number of leaves")+
  theme(legend.position = "none", text = element_text(size=12), plot.title = element_text(size=15),
        plot.subtitle = element_text(size=12))+
  scale_fill_manual(values=c("#e457b5", "#57e486"))

plot_wetweight <- summary_salted %>% #plot the increase in number of leaves grouped by treatment
  ggplot(aes(x = Size, y = mean_wetweight))+
  geom_col(aes(fill = Size))+
  geom_errorbar(aes(ymin = mean_wetweight - stdevwetweight, ymax = mean_wetweight + stdevwetweight), width=.2)+
  theme_light()+
  labs(title = "Wet weight of the plants,\ngrouped by population size",
       subtitle = "Error bars depict 1 standard deviation",
       x="Population size",
       y="Mean wet weight in grams")+
  theme(legend.position = "none", text = element_text(size=12), plot.title = element_text(size=15),
        plot.subtitle = element_text(size=12))+
  scale_fill_manual(values=c("#e457b5", "#57e486"))

plot_grid(plot_length, plot_wetweight, plot_noleaves, #combine the 3 plots into 1 figure
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)

data_color <- data_complete %>% #prepare color data for plotting
  mutate(Color_graph = case_when(
    Color == 1 ~ "Green",
    Color == 2 ~ "Yellow",
    Color == 3 ~ "Red",
    Color == 4 ~ "Red and\nyellow"))

plot_color <- data_color %>% #plot number of leaves of each color, grouped by treatment
  ggplot()+ 
  geom_bar(aes(x = Color_graph, fill = Treatment), position = position_dodge2(width = 0.9, preserve = "single"))+
  theme_light()+
  labs(title = "Number of leaves of each color,\ngrouped by treatment",
       x="Color",
       y="Count")+

plot_shape_salt <- data_complete %>% #plot number of leaves of each shape, grouped by treatment
  ggplot()+ 
  geom_bar(aes(x = Leaf_shape, fill = Treatment), position = position_dodge2(width = 0.9, preserve = "single"))+
  theme_light() +
  labs(title = "Number of leaves of each shape,\ngrouped by treatment",
       x="Leaf shape",
       y="Count")

plot_shape_pop <- data_complete %>% #plot number of leaves of each shape, grouped by population size
  ggplot()+ 
  geom_bar(aes(x = Leaf_shape, fill = Size), position = position_dodge2(width = 0.9, preserve = "single"))+
  theme_light() +
  labs(title = "Number of leaves of each shape,\ngrouped by population size",
       x="Leaf shape",
       y="Count")+
  scale_fill_manual(values=c("#e457b5", "#57e486"))+
  guides(fill=guide_legend(title="Population size"))

img_leafshapes <- image_read(here("images/leaf_shapes.png")) %>% image_ggplot() #import the picture showing the different leaf shapes

plot_grid(plot_color, img_leafshapes, plot_shape_salt, plot_shape_pop, #combine the 3 plots and the leaf shape photo into 1 figure
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)
