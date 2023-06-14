library(readxl)
library(tidyverse)
library(gdata)
library(car)
library(cowplot)
library(magick)
library(here)
library(FSA)
library(agricolae)

##Import your long format data
data <- structure(list(model = structure(c(5L, 5L, 5L, 10L, 10L, 10L, 15L, 15L, 15L, 20L, 20L, 20L, 25L, 25L, 25L, 30L, 30L, 30L), .Label = c("g121", "g122", "g123", "g124", "g125", "g171", "g172", "g173", "g174", "g175", "n121", "n122", "n123", "n124", "n125", "n171", "n172", "n173", "n174", "n175", "r121", "r122", "r123", "r124", "r125", "r171", "r172", "r173", "r174", "r175"), class = "factor"), metric = structure(c(5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L), .Label = c("a", "k", "r", "ss", "sp"), class = "factor"), value = c(0.185, 0.248, 0.253, 0.284, 0.327, 0.286, 0.16, 0.222, 0.217, 0.29, 0.247, 0.255, 0.247, 0.248, 0.222, 0.21, 0.23, 0.278)), class = "data.frame", row.names = c(NA, -18L))
##duplicate the dataframe
data1<-as.data.frame(data)
##we perform ANOVA and compare means using the Tukeys HSD
aov <- aov(value~model, data=data1) %>%
  HSD.test("model", console=TRUE, group=TRUE) %>%
  .$groups %>% as_tibble(rownames="model")
##use the function below to reorder the ANOVA output such that the smallest mean is assigned the letter a and so forth for the next means
reorder<-function(inV){
  collapsed <- paste(inV,sep="",collapse = "")
  u <- unique(strsplit(collapsed,"")[[1]])
  if(length(u)<2){
    return(inV)
  }
  u <- u[order(u)]
  m <- matrix(nrow=NROW(inV),ncol=length(u))
  m[]<-F
  for(i in 1:length(inV)){
    s <- strsplit(inV[i],"")[[1]]
    index <- match(s,u)
    m[i,index] <- T
  }
  for(i in 1:(length(u)-1)){
    firstColT <- match(T,m[,i])[1] #first row with true in current column
    firstT <- match(T,rowSums(m[,i:length(u)] > 0))[1] #first row with true in rest
    if(firstT < firstColT){
      colT <- match(T,m[firstT,i:length(u)])[1]
      colT <- colT + i - 1 #correct index for leftout columns in match
      tmp <- m[,colT]
      m[,colT] <- m[,i]
      m[,i] <- tmp
    }
  }
  res <- vector(mode = "character", length=length("trt"))
  for(i in 1:length(inV)){
    l <- u[m[i,]]
    res[i] <- paste(l,sep="",collapse = "")
  }
  return(res)
}
aov2 <- aov[rev(rownames(aov)),] #order the result the way you want
aov2$groups <- reorder(as.character(aov2$groups))
aov2
#PLOT BOXPLOT WITH LETTERS
ggplot(data1,aes(x=model,y=value)) + 
  geom_boxplot()+
  stat_summary(fun=mean, geom="point", size=2, color="red")+labs(x="\nModels", y="Specificity\n") +
  theme_bw(base_family = "Arial")+ 
  labs_pubr(base_size=12)+
  theme(legend.position="none")+
  geom_text(data=aov2,aes(x=model,y=value,label=groups),vjust=-1)





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

data_salted <- subset(data_complete, Treatment == "High salinity")

##duplicate the dataframe
data1<-as.data.frame(data_complete)
##we perform ANOVA and compare means using the Tukeys HSD
aov <- aov(Number_dif ~ Population, data=data_complete) %>%
  HSD.test("Population", console=TRUE, group=TRUE) %>%
  .$groups %>% as_tibble(rownames="Population")

aov(Number_dif ~ Population*Treatment, data=data_complete) %>% summary()


summary_complete_pop <- data_salted %>%
  group_by(Population) %>%
  summarize(mean_noleaves.increase = mean(Number_dif, na.rm = TRUE),
            mean_lengthleaves.increase = mean(Length_dif, na.rm = TRUE),
            mean_wetweight = mean(Wet_weight, na.rm = TRUE),
            stdevlengthleave.increase = sd(Length_dif, na.rm = TRUE),
            stdevwetweight = sd(Wet_weight, na.rm = TRUE),
            stdevNumber = sd(Number_dif, na.rm = TRUE))

summary_complete_pop %>% #plot the increase in number of leaves grouped by population
  ggplot(aes(x = Population, y = mean_noleaves.increase))+
  geom_col(aes(fill = Population), color = c("#e457b5", "#e457b5", "#57e486", "#57e486"), linewidth = 2)+
  geom_errorbar(aes(ymin = mean_noleaves.increase - stdevNumber,
                    ymax = mean_noleaves.increase + stdevNumber), width=.2)+
  theme_light()+
  labs(title = "Increase in number of leaves,\ngrouped by population",
       subtitle = "Error bars depict 1 standard deviation,\nplot only showing high salinity groups",
       x="Population",
       y="Mean increase in the number of leaves")+
  theme(legend.position = "none", text = element_text(size=12), plot.title = element_text(size=15),
        plot.subtitle = element_text(size=12))+
  scale_fill_manual(values=c("#8111ee", "#7eee11", "#f0860f", "#0f79f0"))+
  geom_text(data=aov,aes(x=Population,y=Number_dif,label=groups),vjust=-1, hjust=3)
