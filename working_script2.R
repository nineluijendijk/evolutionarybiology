library(readxl)
library(tidyverse)
library(gdata)
library(car)
library(cowplot)
library(magick)
library(here)
library(agricolae)
library(ggsignif)

data_raw <- read_excel(here::here("data/arabidopsis_data.xlsx"), col_types = c("text", rep("numeric", 8), "text")) #making sure every value is numeric

metadata <- tibble("Group" = c("L1", "L2", "S3", "S4", "L5", "L6", "S7", "S8"),
                   "Population" = rep(c("Bovra", "Smadalen", "Helin", "Spiterstulen"), each = 2),
                   "Size" = rep(c("Large", "Small", "Large", "Small"), each = 2),
                   "Elevation" = rep("High", each = 8),
                   "Condition" = rep(c("Stable", "Declining", "Declining", "Declining"), each = 2), #difference between "decline" and "decrease"? https://brightspace.ru.nl/d2l/le/content/345973/viewContent/1977818/View
                   "Soil" = rep(c("Riverbed", "Riverbed", "Scree", "Rocks"), each = 2),
                   "Treatment" = c("Low salinity", "High salinity", "Low salinity", "High salinity", "Low salinity", "High salinity", "Low salinity", "High salinity"),
                   "Conc_NaCl" = c(0, 50, 0, 50, 0, 50, 0, 50),
                   "Conc_unit" = rep("mM", each = 8),
                   "Length_unit" = rep("cm", each = 8),
                   "Weight_unit" = rep("grams"), each = 8)

data <- mutate(data_raw,
               "Length_dif" = data_raw$Length_longest_leaf_a - data_raw$Length_longest_leaf_b,
               "Number_dif" = data_raw$Number_leaves_a - data_raw$Number_leaves_b) #calculating differences between before and after treatment

data_complete <- left_join(data, metadata, by = "Group") #create one complete dataframe

summary <- data_complete %>%
  group_by(Treatment) %>%
  summarize(mean_noleaves.increase = mean(Number_dif, na.rm = TRUE),
            mean_lengthleaf.increase = mean(Length_dif, na.rm = TRUE),
            mean_wetweight = mean(Wet_weight, na.rm = TRUE),
            stdevLength = sd(Length_dif, na.rm = TRUE),
            stdevWetweight = sd(Wet_weight, na.rm = TRUE),
            stdevNumber = sd(Number_dif, na.rm = TRUE))

data_salted <- subset(data_complete, Treatment == "High salinity")

summary_salted <- data_salted %>%
  group_by(Size) %>%
  summarize(mean_noleaves.increase = mean(Number_dif, na.rm = TRUE),
            mean_lengthleaves.increase = mean(Length_dif, na.rm = TRUE),
            mean_wetweight = mean(Wet_weight, na.rm = TRUE),
            stdevlengthleave.increase = sd(Length_dif, na.rm = TRUE),
            stdevwetweight = sd(Wet_weight, na.rm = TRUE),
            stdevNumber = sd(Number_dif, na.rm = TRUE))

summary_complete_pop <- data_salted %>%
  group_by(Population) %>%
  summarize(mean_noleaves.increase = mean(Number_dif, na.rm = TRUE),
            mean_lengthleaves.increase = mean(Length_dif, na.rm = TRUE),
            mean_wetweight = mean(Wet_weight, na.rm = TRUE),
            stdevlengthleave.increase = sd(Length_dif, na.rm = TRUE),
            stdevwetweight = sd(Wet_weight, na.rm = TRUE),
            stdevNumber = sd(Number_dif, na.rm = TRUE))





# ANALYSIS OF LENGTH DIF

leveneTest(data_complete$Length_dif ~ as.factor(Treatment), data = data_complete) #Pr(>F) = 0.4793, which means equal variance

anno <- t.test(formula = data_complete$Length_dif ~ data_complete$Treatment, 
       paired = FALSE, var.equal = TRUE)$p.value %>% round(.,3) #t-test, p = 0.005, which means there is a statistically significant difference in increase in longest leaf length between the High salinity and Low salinity plants.

plot_length1 <- summary %>% #plot the increase in longest leaf length grouped by population size
  ggplot(aes(x = Treatment, y = mean_lengthleaf.increase))+
  geom_col(aes(fill = Treatment))+
  geom_errorbar(aes(ymin = mean_lengthleaf.increase - stdevLength,
                    ymax = mean_lengthleaf.increase + stdevLength), width=.2)+
  theme_light()+
  labs(x="Treatment",
       y="Mean increase in the length\nof the longest leaf in cm")+
  theme(legend.position = "none", text = element_text(size=12), plot.title = element_text(size=14),
        plot.subtitle = element_text(size=12))+
  geom_signif(comparisons = list(c("High salinity", "Low salinity")),
              annotation = formatC(anno, digits = 2))

leveneTest(data_salted$Length_dif ~ as.factor(Size), data = data_salted) #Pr(>F) = 0.7135 which means equal variance

anno <- t.test(formula = data_salted$Length_dif ~ data_salted$Size, 
       paired = FALSE, var.equal = TRUE)$p.value %>% round(.,3) #t-test, p = 0.563, which means there is no statistically significant difference between the difference in length of the longest leaf of the large and small High salinity plants.

plot_length2 <- summary_salted %>% #plot the increase in longest leaf length grouped by population size
  ggplot(aes(x = Size, y = mean_lengthleaves.increase))+
  geom_col(aes(fill = Size))+
  geom_errorbar(aes(ymin = mean_lengthleaves.increase - stdevlengthleave.increase,
                    ymax = mean_lengthleaves.increase + stdevlengthleave.increase), width=.2)+
  theme_light()+
  labs(x="Population size")+
  theme(legend.position = "none", text = element_text(size=12), plot.title = element_text(size=14),
        plot.subtitle = element_text(size=12))+
  ylab(NULL)+
  scale_fill_manual(values=c("#e457b5", "#57e486"))+
  geom_signif(comparisons = list(c("Large", "Small")),
              annotation = formatC(anno, digits = 2))

leveneTest(data_salted$Length_dif ~ as.factor(Population), data = data_salted) #Pr(>F) = 0.8361, which means equal variance

##we perform ANOVA and compare means using the Tukeys HSD
aov <- aov(Length_dif ~ Population, data=data_salted) %>%
  HSD.test("Population", console=TRUE, group=TRUE) %>%
  .$groups %>% as_tibble(rownames="Population")

plot_length3 <- summary_complete_pop %>% #plot the increase in longest leaf length grouped by population
  ggplot(aes(x = Population, y = mean_lengthleaves.increase))+
  geom_col(aes(fill = Population), color = c("#e457b5", "#e457b5", "#57e486", "#57e486"), linewidth = 2)+
  geom_errorbar(aes(ymin = mean_lengthleaves.increase - stdevlengthleave.increase,
                    ymax = mean_lengthleaves.increase + stdevlengthleave.increase), width=.2)+
  theme_light()+
  labs(x="Population")+
  theme(legend.position = "none", text = element_text(size=12), plot.title = element_text(size=14),
        plot.subtitle = element_text(size=12))+
  ylab(NULL)+
  scale_fill_manual(values=c("#8111ee", "#7eee11", "#f0860f", "#0f79f0"))+
  geom_text(data=aov,aes(x=Population,y=Length_dif,label=groups),vjust=-1, hjust=3)


plotgridlength <- plot_grid(plot_length1, plot_length2, plot_length3,
                      labels = c("A", "B", "C"),
                      ncol = 3, nrow = 1)


# ANALYSIS OF WET WEIGHT DIF

leveneTest(data_complete$Wet_weight ~ as.factor(Treatment), data = data_complete) #Pr(>F) = 0.4168, which means equal variance

anno <- t.test(formula = data_complete$Wet_weight ~ data_complete$Treatment, 
       paired = FALSE, var.equal = TRUE)$p.value %>% round(.,3) #t-test, p = 0.011, which means there is a statistically significant difference in wet weight between the High salinity and Low salinity plants.

plot_weight1 <- summary %>% #plot the wet weight grouped by population size
  ggplot(aes(x = Treatment, y = mean_wetweight))+
  geom_col(aes(fill = Treatment))+
  geom_errorbar(aes(ymin = mean_wetweight - stdevWetweight,
                    ymax = mean_wetweight + stdevWetweight), width=.2)+
  theme_light()+
  labs(x="Treatment",
       y="Mean net weight in grams")+
  theme(legend.position = "none", text = element_text(size=12), plot.title = element_text(size=14),
        plot.subtitle = element_text(size=12))+
  geom_signif(comparisons = list(c("High salinity", "Low salinity")),
              annotation = formatC(anno, digits = 2))

leveneTest(data_salted$Wet_weight ~ as.factor(Size), data = data_salted) #Pr(>F) = 0.3754 which means equal variance

anno <- t.test(formula = data_salted$Wet_weight ~ data_salted$Size, 
       paired = FALSE, var.equal = TRUE)$p.value %>% round(.,3) #t-test, p = 0.848, which means there is no statistically significant difference in wet weight between the large and small High salinity plants.


plot_weight2 <-summary_salted %>% #plot the increase in number of leaves grouped by population
  ggplot(aes(x = Size, y = mean_wetweight))+
  geom_col(aes(fill = Size))+
  geom_errorbar(aes(ymin = mean_wetweight - stdevwetweight,
                    ymax = mean_wetweight + stdevwetweight), width=.2)+
  theme_light()+
  labs(x="Population size")+
  theme(legend.position = "none", text = element_text(size=12), plot.title = element_text(size=15),
        plot.subtitle = element_text(size=12))+
  ylab(NULL)+
  scale_fill_manual(values=c("#e457b5", "#57e486"))+
  geom_signif(comparisons = list(c("Large", "Small")),
              annotation = formatC(anno, digits = 2))

leveneTest(data_salted$Wet_weight ~ as.factor(Population), data = data_salted) #Pr(>F) = 0.7846, which means equal variance

##we perform ANOVA and compare means using the Tukeys HSD
aov <- aov(Wet_weight ~ Population, data=data_salted) %>%
  HSD.test("Population", console=TRUE, group=TRUE) %>%
  .$groups %>% as_tibble(rownames="Population")

plot_weight3 <- summary_complete_pop %>% #plot the increase in number of leaves grouped by population
  ggplot(aes(x = Population, y = mean_wetweight))+
  geom_col(aes(fill = Population), color = c("#e457b5", "#e457b5", "#57e486", "#57e486"), linewidth = 2)+
  geom_errorbar(aes(ymin = mean_wetweight - stdevwetweight,
                    ymax = mean_wetweight + stdevwetweight), width=.2)+
  theme_light()+
  labs(x="Population")+
  theme(legend.position = "none", text = element_text(size=12), plot.title = element_text(size=15),
        plot.subtitle = element_text(size=12))+
  ylab(NULL)+
  scale_fill_manual(values=c("#8111ee", "#7eee11", "#f0860f", "#0f79f0"))+
  geom_text(data=aov,aes(x=Population,y=Wet_weight,label=groups),vjust=-1, hjust=3)

plotgridweight <- plot_grid(plot_weight1, plot_weight2, plot_weight3,
                            labels = c("A", "B", "C"),
                            ncol = 3, nrow = 1)



# ANALYSIS OF NUMBER DIF

leveneTest(data_complete$Number_dif ~ as.factor(Treatment), data = data_complete) #Pr(>F) = 0.6024, which means equal variance

anno <- t.test(formula = data_complete$Number_dif ~ data_complete$Treatment, 
       paired = FALSE, var.equal = TRUE)$p.value %>% round(.,3) #t-test, p = 0.868, which means there is a statistically significant difference in wet weight between the High salinity and Low salinity plants.

plot_number1 <- summary %>% #plot the increase in number of leaves grouped by population size
  ggplot(aes(x = Treatment, y = mean_noleaves.increase))+
  geom_col(aes(fill = Treatment))+
  geom_errorbar(aes(ymin = mean_noleaves.increase - stdevNumber,
                    ymax = mean_noleaves.increase + stdevNumber), width=.2)+
  theme_light()+
  labs(x="Treatment",
       y="Mean increase in the number of leaves")+
  theme(legend.position = "none", text = element_text(size=12), plot.title = element_text(size=14),
        plot.subtitle = element_text(size=12))+
  geom_signif(comparisons = list(c("High salinity", "Low salinity")),
              annotation = formatC(anno, digits = 2))

leveneTest(data_salted$Number_dif ~ as.factor(Size), data = data_salted) #Pr(>F) = 0.7076 which means equal variance

anno <- t.test(formula = data_salted$Number_dif ~ data_salted$Size, 
               paired = FALSE, var.equal = TRUE)$p.value %>% round(.,3)

plot_number2 <- summary_salted %>% #plot the increase in number of leaves grouped by population
  ggplot(aes(x = Size, y = mean_noleaves.increase))+
  geom_col(aes(fill = Size))+
  geom_errorbar(aes(ymin = mean_noleaves.increase - stdevNumber,
                    ymax = mean_noleaves.increase + stdevNumber), width=.2)+
  theme_light()+
  labs(x="Population size")+
  theme(legend.position = "none", text = element_text(size=12), plot.title = element_text(size=15),
        plot.subtitle = element_text(size=12))+
  ylab(NULL)+
  scale_fill_manual(values=c("#e457b5", "#57e486"))+
  geom_signif(comparisons = list(c("Large", "Small")),
              annotation = formatC(anno, digits = 2))

leveneTest(data_salted$Number_dif ~ as.factor(Population), data = data_salted) #Pr(>F) = 0.1503, which means equal variance

##we perform ANOVA and compare means using the Tukeys HSD
aov <- aov(Number_dif ~ Population, data=data_salted) %>%
  HSD.test("Population", console=TRUE, group=TRUE) %>%
  .$groups %>% as_tibble(rownames="Population")

pairwise.t.test(data_salted$Number_dif, data_salted$Population) #p-value is 0.017

plot_number3 <- summary_complete_pop %>% #plot the increase in number of leaves grouped by population
  ggplot(aes(x = Population, y = mean_noleaves.increase))+
  geom_col(aes(fill = Population), color = c("#e457b5", "#e457b5", "#57e486", "#57e486"), linewidth = 2)+
  geom_errorbar(aes(ymin = mean_noleaves.increase - stdevNumber,
                    ymax = mean_noleaves.increase + stdevNumber), width=.2)+
  theme_light()+
  labs(x="Population")+
  theme(legend.position = "none", text = element_text(size=12), plot.title = element_text(size=15),
        plot.subtitle = element_text(size=12))+
  ylab(NULL)+
  scale_fill_manual(values=c("#8111ee", "#7eee11", "#f0860f", "#0f79f0"))+
  geom_text(data=aov,aes(x=Population,y=Number_dif,label=groups),vjust=-1, hjust=1.5)+
  geom_signif(comparisons = list(c("Bovra", "Helin")),
              annotation = formatC(0.017, digits = 2))

plotgridnumber <- plot_grid(plot_number1, plot_number2, plot_number3,
                            labels = c("A", "B", "C"),
                            ncol = 3, nrow = 1)








# ANALYZING COLOR DIFFERENCE

chi1 <- chisq.test(data_complete$Treatment, data_complete$Color, correct=FALSE)$p.value %>% round(digits = 3) # p = 0.8044
chi2 <- chisq.test(data_salted$Size, data_salted$Color, correct=FALSE)$p.value  %>% round(digits = 3) # p = 0.3334
chi3 <- chisq.test(data_salted$Population, data_salted$Color, correct=FALSE)$p.value  %>% round(digits = 3) # p = 0.01822

data_color <- data_complete %>% #prepare color data for plotting
  mutate(Color_graph = case_when(
    Color == 1 ~ "Green",
    Color == 2 ~ "Yellow",
    Color == 3 ~ "Red",
    Color == 4 ~ "Red and\nyellow"))

plot_color1 <- data_color %>% #plot number of leaves of each color, grouped by treatment
  ggplot()+ 
  geom_bar(aes(x = Color_graph, fill = Treatment), position = position_dodge2(width = 0.9, preserve = "single"))+
  theme_light()+
  theme(legend.position = "none")+
  labs(y="Count")+
  xlab(NULL)+
  annotate("text", x = 3, y = 22.5, 
           label = paste("Chi-squared p-value =", chi1))

plot_color2 <- data_color %>% #plot number of leaves of each color, grouped by population
  ggplot()+ 
  geom_bar(aes(x = Color_graph, fill = Size), linewidth = 1,
           position = position_dodge2(width = 0.9, preserve = "single"))+
  theme_light()+
  theme(legend.position = "none")+
  xlab(NULL)+
  ylab(NULL)+
  scale_fill_manual(values = c("#e457b5", "#57e486"))+
  guides(color=guide_legend(title="Population size"))+
  annotate("text", x = 3, y = 22.5, 
           label = paste("Chi-squared p-value =", chi2))

plot_color3 <- data_color %>% #plot number of leaves of each color, grouped by population
  ggplot()+ 
  geom_bar(aes(x = Color_graph, fill = Population, color = Size), linewidth = 1,
           position = position_dodge2(width = 0.9, preserve = "single"))+
  theme_light()+
  theme(legend.position = "none")+
  xlab(NULL)+
  ylab(NULL)+
  scale_fill_manual(values=c("#8111ee", "#7eee11", "#f0860f", "#0f79f0"))+
  scale_color_manual(values = c("#e457b5", "#57e486"))+
  guides(color=guide_legend(title="Population size"))+
  annotate("text", x = 3, y = 13.1, 
           label = paste("Chi-squared p-value =", chi3))



plot_color01 <- data_color %>% #plot number of leaves of each color, grouped by treatment
  ggplot()+ 
  geom_bar(aes(x = Color_graph, fill = Treatment), position = position_dodge2(width = 0.9, preserve = "single"))+
  theme_light()+
  labs(x="Color",
       y="Count")

plot_color02 <- data_color %>% #plot number of leaves of each color, grouped by population
  ggplot()+ 
  geom_bar(aes(x = Color_graph, fill = Size), linewidth = 1,
           position = position_dodge2(width = 0.9, preserve = "single"))+
  labs(fill = "Population size")+
  theme_light()+
  xlab(NULL)+
  ylab(NULL)+
  scale_fill_manual(values = c("#e457b5", "#57e486"))+
  annotate("text", x = 3, y = 22.5, 
           label = paste("Chi-squared p-value =", chi2))

plot_color03 <- data_color %>% #plot number of leaves of each color, grouped by population
  ggplot()+ 
  geom_bar(aes(x = Color_graph, fill = Population), linewidth = 1,
           position = position_dodge2(width = 0.9, preserve = "single"))+
  theme_light()+
  labs(x="Color",
       y="Count")+
  ylab(NULL)+
  scale_fill_manual(values=c("#8111ee", "#7eee11", "#f0860f", "#0f79f0"))+
  scale_color_manual(values = c("#e457b5", "#57e486"))+
  guides(color=guide_legend(title="Population size"))

legend1 <- get_legend(plot_color01)
legend2 <- get_legend(plot_color02)
legend3 <- get_legend(plot_color03)


plotgridcolor0 <- plot_grid(plot_color1, plot_color2, plot_color3,
          plot_grid(legend1, legend2, legend3, ncol = 1, nrow = 3),
                      labels = c("A", "B", "C"),
                      ncol = 4, nrow = 1, rel_widths = c(1, 1, 1, 0.4))

plotgridcolor <- ggdraw(add_sub(plotgridcolor0, "Color", vpadding=grid::unit(0,"lines"),y=6, x=0.5, vjust=6))








plot_shape_salt <- data_complete %>% #plot number of leaves of each shape, grouped by treatment
  ggplot()+ 
  geom_bar(aes(x = Leaf_shape, fill = Treatment), position = position_dodge2(width = 0.9, preserve = "single"))+
  theme_light() +
  labs(title = "Number of leaves of each shape,\ngrouped by treatment",
       x="Leaf shape",
       y="Count")

plot_shape_popsize <- data_salted %>% #plot number of leaves of each shape, grouped by population size
  ggplot()+ 
  geom_bar(aes(x = Leaf_shape, fill = Size), position = position_dodge2(width = 0.9, preserve = "single"))+
  theme_light() +
  labs(title = "Number of leaves of each shape,\ngrouped by population size",
       subtitle = "Plot only showing high salinity groups",
       x="Leaf shape",
       y="Count")+
  scale_fill_manual(values=c("#e457b5", "#57e486"))+
  guides(fill=guide_legend(title="Population size"))

plot_shape_pop <- data_salted %>% #plot number of leaves of each shape, grouped by population
  ggplot()+ 
  geom_bar(aes(x = Leaf_shape, fill = Population, color = Size), linewidth = 2,
           position = position_dodge2(width = 0.9, preserve = "single"))+
  theme_light() +
  labs(title = "Number of leaves of each shape,\ngrouped by population",
       subtitle = "Plot only showing high salinity groups",
       x="Leaf shape",
       y="Count")+
  scale_fill_manual(values=c("#8111ee", "#7eee11", "#f0860f", "#0f79f0"))+
  scale_color_manual(values = c("#e457b5", "#57e486"))+
  guides(color=guide_legend(title="Population size"))




plotgrid <- plot_grid(plot_color, plot_color2, plot_shape_salt, img_leafshapes, plot_shape_popsize, plot_shape_pop, #combine the 3 plots and the leaf shape photo into 1 figure
                       labels = c("A", "B", "C", "D", "E", "F"),
                       ncol = 2, nrow = 3)



aov(Length_dif ~ Population*Treatment, data=data_complete) %>% summary() # p = 0.5056
aov(Wet_weight ~ Population*Treatment, data=data_complete) %>% summary() # p = 0.6773
aov(Number_dif ~ Population*Treatment, data=data_complete) %>% summary() # p = 0.845


chisq.test(data_complete$Treatment, data_complete$Leaf_shape, correct=FALSE) # p = 0.003457
chisq.test(data_salted$Size, data_salted$Leaf_shape, correct=FALSE) # p = 0.01702
chisq.test(data_salted$Population, data_salted$Leaf_shape, correct=FALSE) # p = 0.0002744


plot_color01 <- data_color %>% #plot number of leaves of each color, grouped by treatment
  ggplot()+ 
  geom_bar(aes(x = Color_graph, fill = Treatment), position = position_dodge2(width = 0.9, preserve = "single"))+
  theme_light()+
  labs(x="Color",
       y="Count")

plot_color03 <- data_color %>% #plot number of leaves of each color, grouped by population
  ggplot()+ 
  geom_bar(aes(x = Color_graph, fill = Population, color = Size), linewidth = 1,
           position = position_dodge2(width = 0.9, preserve = "single"))+
  theme_light()+
  labs(x="Color",
       y="Count")+
  ylab(NULL)+
  scale_fill_manual(values=c("#8111ee", "#7eee11", "#f0860f", "#0f79f0"))+
  scale_color_manual(values = c("#e457b5", "#57e486"))+
  guides(color=guide_legend(title="Population size"))


