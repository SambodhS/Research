library(tidyverse)
library(openxlsx2)
library(mixOmics)

setwd("~/Downloads")
set.seed(123)

data <- read_xlsx("~/Downloads/Research/Microbial/Data/Metabolomics Layer of PAMPer.xlsx")
responder <- read_xlsx("~/Downloads/Research/Microbial/Data/Responder.xlsx")
metabolites <- read_xlsx("~/Downloads/Research/Microbial/Data/Metabolites.xlsx")
origins <- read_xlsx("~/Downloads/Research/Microbial/Data/Origins.xlsx")

#### Cleaning ####

names <- data$BIOCHEMICAL
kegg <- dplyr::select(data, c(2, 12))
pathways <- dplyr::select(data, c(2:4))
data <- dplyr::select(data, -c(1:13)) %>%
  t() %>%
  as.data.frame()
data[,1] <- replace(data[,1], is.na(data[,1]), 0)
colnames(data) <- names
colnames(data)[1] <- "TIME"

CClusters <- responder %>%
  mutate(RESPONDER = case_when(
    !is.na(tdiff) & tdiff <= 3 ~ "EarlyNon-Survivors",
    (!is.na(tdiff) & tdiff > 3) | icu_los >= 7 ~ "NonResolvers",
    is.na(tdiff) ~ "Resolvers"
  ))
data_row_names <- substr(rownames(data), 1, 8)
matching_rows <- match(data_row_names, CClusters$`Study ID`)
data$RESPONDER <- NA
matching_indices <- !is.na(matching_rows)
data$RESPONDER[matching_indices] <- CClusters$RESPONDER[matching_rows[matching_indices]]
Responder <- data$RESPONDER

data <- data[,-c(900:901)] %>%
  mutate(STATUS = case_when(
    startsWith(rownames(data), "M") ~ "Control",
    startsWith(rownames(data), "P") ~ "Trauma")) %>%
  mutate(GROUP = ifelse(STATUS == "Control", "HC", "Trauma")) %>%
  mutate(GROUP = ifelse(STATUS == "Trauma" & TIME == 0, "T0", GROUP)) %>%
  mutate(GROUP = ifelse(STATUS == "Trauma" & TIME == 24, "T24", GROUP)) %>%
  mutate(GROUP = ifelse(STATUS == "Trauma" & TIME == 72, "T72", GROUP)) %>%
  dplyr::select(-c(1, 900))
groups <- data$GROUP
data <- data[, -c(899, 900)] %>%
  t() %>%
  log2() %>%
  scale() %>%
  as.data.frame()

microbial <- data %>%
  rownames_to_column() %>%
  right_join(metabolites, by = c("rowname" = "metabolites")) %>%
  column_to_rownames("rowname") %>%
  t() %>%
  as.data.frame()
microbial$GROUP <- groups
microbial$RESPONDER <- Responder

whole <- data %>%
  t() %>%
  as.data.frame()
whole$GROUP <- groups
whole$RESPONDER <- Responder

#### sPLS-DA #### 

# Assigning MetOrigin-derived origins to top contributing metabolites of sPLS-DA model

test <- whole %>%
  filter(GROUP == "T0" & !is.na(RESPONDER))
test_groups <- test$RESPONDER
test <- test[,-c(899:900)]

model <- mixOmics::splsda(test, test_groups)
plotIndiv(model, style = "ggplot2", ind.names = F, ellipse = T, legend = T, title = "PAMPer Metabolome sPLS-DA")
auroc(model)
loadings <- plotLoadings(model, 
                         contrib = 'max', 
                         method = 'mean', 
                         title = " ")
contributions <- as.data.frame(loadings$X) %>%
  rownames_to_column() %>%
  left_join(origins, by = c("rowname" = "BIOCHEMICAL")) %>%
  left_join(pathways, by = c("rowname" = "BIOCHEMICAL")) %>%
  dplyr::select(c(1, 9, 11, 14, 18)) %>%
  filter(GroupContrib == "EarlyNon-Survivors") %>%
  arrange(desc(importance))

contributions <- contributions[,(c(1, 3, 4, 2, 5))]
write_xlsx(contributions, "contributions.xlsx")
