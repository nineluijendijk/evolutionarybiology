library(readxl)
library(tidyverse)

data <- read_excel("data/arabidopsis_data.xlsx", col_types = c("text", rep("numeric", 8), "text"))

           