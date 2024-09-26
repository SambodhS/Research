library(tidyverse)
library(openxlsx2)
library(mixOmics)
library(caret)
library(ROCR)

setwd("~/Downloads")

data <- read_xlsx("~/Downloads/Research/Microbial/Data/Metabolomics Layer of PAMPer.xlsx")
outcomes <- read_xlsx("~/Downloads/Research/Microbial/Data/Outcomes.xlsx")
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

CClusters <- outcomes %>%
  mutate(RESPONDER = case_when(
    !is.na(tdiff) & tdiff <= 3 ~ "EarlyNon-Survivors",
    (!is.na(tdiff) & tdiff > 3) | icu_los >= 7 ~ "NonResolvers",
    is.na(tdiff) ~ "Resolvers" # & !is.na(icu_los) & icu_los < 7 ????
  )) %>%
  mutate(TDM = case_when(
    !is.na(alive_at_30) & alive_at_30 == 1 ~ "Alive",
    !is.na(alive_at_30) & alive_at_30 == 2 ~ "Dead",
    alive_at_30 == 3 | is.na(alive_at_30) ~ NA
  ))
data_row_names <- substr(rownames(data), 1, 8)
matching_rows <- match(data_row_names, CClusters$`Study ID`)
data$RESPONDER <- NA
matching_indices <- !is.na(matching_rows)
data$RESPONDER[matching_indices] <- CClusters$RESPONDER[matching_rows[matching_indices]]
Responder <- data$RESPONDER
data$TDM <- NA
data$TDM[matching_indices] <- CClusters$TDM[matching_rows[matching_indices]]
TDM <- data$TDM

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
microbial$TDM <- TDM

whole <- data %>%
  t() %>%
  as.data.frame()
whole$GROUP <- groups
whole$RESPONDER <- Responder
whole$TDM <- TDM

#### sPLS-DA #### 

# Micromets classification

set <- microbial %>%
  filter(GROUP == "T72" & !is.na(RESPONDER))
set_groups <- set$RESPONDER
set <- set[,-c(151:152)]

model <- mixOmics::splsda(set, set_groups)
plotIndiv(model, style = "ggplot2", ind.names = F, ellipse = T, legend = T, title = "PAMPer Metabolome sPLS-DA")
auroc(model)
loadings <- plotLoadings(model, 
                         contrib = 'max', 
                         method = 'mean', 
                         title = " ")
contributions <- as.data.frame(loadings$X) %>%
  arrange(desc(importance))

# Assigning MetOrigin-derived origins to top contributing metabolites for outcomes

test <- whole %>%
  filter(GROUP == "T72" & !is.na(RESPONDER))
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

# Train/test split for outcome prediction performance

set <- microbial %>%
  filter(GROUP == "T72" & !is.na(RESPONDER))
set_groups <- set$RESPONDER
train_idx <- createDataPartition(set$RESPONDER, times=1, p=0.8, list=FALSE)
train <- set[train_idx,]
train_groups <- train$RESPONDER
train <- train[,-c(151:152)]
test <- set[-train_idx,]
test_groups <- test$RESPONDER
test <- test[,-c(151:152)]

model <- mixOmics::splsda(train, train_groups)
loadings <- plotLoadings(model, 
                         contrib = 'max', 
                         method = 'mean', 
                         title = " ")
contributions <- as.data.frame(loadings$X) %>%
  arrange(desc(importance)) %>%
  dplyr::select(8)
auroc(model)
auroc(model, test, test_groups)

#### Logistic Regression ####

# Using microbial-only metabolites as features in outcome prediction

pure <- data %>%
  rownames_to_column() %>%
  right_join(metabolites, by = c("rowname" = "metabolites")) %>%
  left_join(origins, by = c("rowname" = "BIOCHEMICAL")) %>%
  filter(rowname == "1-hydroxy-2-naphthalenecarboxylate" | rowname == "N-methylproline" |
           rowname == "6-oxopiperidine-2-carboxylate") %>%
  #filter(grepl("Microbiota", Origin) & !grepl("Host", Origin) & !grepl("Food related", Origin) &
  #         !grepl("Drug related", Origin) & !grepl("Environment", Origin)) %>%
  column_to_rownames() %>%
  dplyr::select(-c(504:508)) %>%
  t() %>%
  as.data.frame()
pure$GROUP <- groups
pure$RESPONDER <- Responder
pure <- filter(pure, GROUP == "T72") %>%
  mutate(RESPONDER = ifelse(RESPONDER == "Resolvers", 1, 0)) %>% # resolvers = 1 
  #dplyr::select(-9)
  dplyr::select(-4)
pure <- pure[-136,] # removing patient with NA responder value 
pure$RESPONDER <- as.factor(pure$RESPONDER)
train_idx <- createDataPartition(pure$RESPONDER, times=1, p=0.8, list=FALSE)
train <- pure[train_idx,]
test <- pure[-train_idx,]

model <- glm(RESPONDER ~ ., family=binomial, data=train)
summary(model) 
probs <- data.frame(probs = predict(model, test, type="response")) %>%
  mutate(pred = ifelse(probs > 0.5, "1", "0"))
prediction <- prediction(predict(model, test, type="response"), test$RESPONDER)
performance(prediction, measure="auc")@y.values[[1]]
print(sum(data.frame(probs$pred == test$RESPONDER), na.rm=TRUE))
#perf <- performance(prediction, "tpr", "fpr")
#plot(perf, col="red")

(0.7236842+0.7894737+0.7631579+0.8421053+0.9013158+0.75+0.7039474+0.6907895+0.6381579+0.7960526)/10
((18+20+20+21+21+18+19+20+21+18)/10)/27

# Using sPLS-DA-selected top 20 contributing micromets in T72 outcome classification

select <- data %>%
  rownames_to_column() %>%
  filter(rowname == "pentose acid*")%>%# | rowname == "dihydroferulate" | rowname == "3-(3-hydroxyphenyl)propionate" |
           # rowname == "dimethylglycine" | rowname == "3-hydroxyhippurate" | rowname == "isovalerate (C5)" |
           # rowname == "N-methylproline" | rowname == "stachydrine" | rowname == "S-methylcysteine" |
           # rowname == "glycerol 3-phosphate" | rowname == "alpha-hydroxycaproate" |
           # rowname == "ergothioneine" | rowname == "2-isopropylmalate" | rowname == "salicylate" |
           # rowname == "4-hydroxyglutamate" | rowname == "beta-cryptoxanthin" |
           # rowname == "quinate" | rowname == "ascorbate (Vitamin C)" | rowname == "1,6-anhydroglucose" |
           # rowname == "deoxycholate") %>%
  column_to_rownames() %>%
  t() %>%
  as.data.frame()
select$GROUP <- groups
select$RESPONDER <- Responder
select <- filter(select, GROUP == "T72" & !is.na(RESPONDER)) %>%
  mutate(RESPONDER = ifelse(RESPONDER == "Resolvers", 1, 0)) %>% # resolvers = 1
  #dplyr::select(-21)
  dplyr::select(-2)
train_idx <- createDataPartition(select$RESPONDER, times=1, p=0.8, list=FALSE)
train <- select[train_idx,]
test <- select[-train_idx,]

model <- glm(RESPONDER ~ ., family=binomial, data=train)
summary(model)
probs <- data.frame(probs = predict(model, test, type="response")) %>%
  mutate(pred = ifelse(probs > 0.5, "1", "0"))
prediction <- prediction(predict(model, test, type="response"), test$RESPONDER)
performance(prediction, measure="auc")@y.values[[1]]
print(sum(data.frame(probs$pred == test$RESPONDER), na.rm=TRUE))
perf <- performance(prediction, "tpr", "fpr")
plot(perf, col="red")

(0.9319728+0.7575758+0.9387755+0.8181818+0.7282051+0.6782609+0.8947368+0.8571429+0.8777778+0.58125)/10
((24+24+26+24+18+21+23+23+23+21)/10)/27

