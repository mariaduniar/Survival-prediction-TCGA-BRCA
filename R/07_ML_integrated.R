###############################################################################
#########                   ANÁLISIS INTEGRATIVO               ################
###############################################################################

rm(list = ls())

setwd("C:/Users/Usuario/Desktop/Bioinformatica/TFM")

library(caret)
library(pROC)


# - CARGA DE DATOS ----

# Datos clínicos
clin_data<-readRDS("clinical_survival.rds")

# Matriz de expresión normalizada
exp_data<-readRDS("counts_survival.rds")

# Matriz binaria de vías alteradas 
path_data<-read.csv("mut_data_binary_pathways.csv", row.names = 1)

# Lista genes diferencialmente expresados 
DEG<-read.csv("DESeq2_train/list_DEG.csv")

# Lista vías seleccionadas por LASSO
paths_lasso<-read.csv("pathways_LASSO.csv")


# - Preparación datos ----

# - Filtrado genes DEG
exp_data_DEG<-exp_data[,DEG$x]

# - Filtrado rutas seleccionadas por LASSO 
path_data_subset<-path_data[,paths_lasso$x]

# - Unión 
integration<-merge(exp_data_DEG, path_data_subset, by=0)
rownames(integration)<-integration$Row.names
integration<-integration[,-1]

# - Añadimos Risk group y convertimos a factor

integration$risk_group <- clin_data$risk_group[
  match(rownames(integration), clin_data$submitter_id)
]

table(integration$risk_group)


# - División Train/Test ----

train_ids<-readRDS("train_IDs.rds")
test_ids<-readRDS("test_IDs.rds")


train_data<-integration[train_ids, ]
test_data<-integration[test_ids, ]

# - Comprobamos proporciones 
table(integration$risk_group)
table(train_data$risk_group)
table(test_data$risk_group)


 # Renombramos niveles
train_data$risk_group <- factor(
  train_data$risk_group,
  levels = c("High risk", "Low risk"),
  labels = c("High_risk", "Low_risk")
)

test_data$risk_group <- factor(
  test_data$risk_group,
  levels = c("High risk", "Low risk"),
  labels = c("High_risk", "Low_risk")
)

# - SVM ----

# Creamos objeto para el parámetro trControl
ctrl <- trainControl(
  method = "repeatedcv", #validación cruzada
  number = 10, #10 folds
  repeats = 5,
  classProbs = TRUE, #obtener probabilidades de clase
  savePredictions = "final", #guardar predicciones optimas
  summaryFunction = twoClassSummary #calcular metricas de rendimiento de problemas de 2 clases, es decir ROC, sensibilidad y especificidad
)

## - SVM LINEAL ----
set.seed(123)
svmModelLineal <- train(risk_group ~.,
                        data = train_data,
                        method = "svmLinear",
                        trControl = ctrl,
                        preProcess = NULL,
                        tuneLength = expand.grid(C = c(0.01, 0.1, 1, 10)),
                        prob.model = TRUE,
                        metric="ROC") 

svmModelLineal

# - Predicciones sobre el test
pred_smvl<-predict(svmModelLineal, newdata=test_data)
table(pred_smvl)

confusionMatrix(pred_smvl, test_data$risk_group, mode="everything")

# - Probabilidades de clase
svmModelLineal$pred$High_risk #train

test_prob_svml<-predict(svmModelLineal, newdata = test_data, type = "prob") #test

# - ROC

roc_svml<-roc(
  response = test_data$risk_group, #Clases reales
  predictor = test_prob_svml$High_risk, #Probabilidad de ser high risk
  levels = c("Low_risk", "High_risk")
)

#AUC
auc(roc_svml)
ci.auc(roc_svml)



## - SVM RADIAL ----
set.seed(123)
svmModelRadial <- train(risk_group ~.,
                        data = train_data,
                        method = "svmRadial",
                        trControl = ctrl,
                        preProcess = NULL,
                        tuneLength = 4, 
                        prob.model = TRUE,
                        metric="ROC") 

svmModelRadial
plot(svmModelRadial)

# - Predicciones sobre el test
pred_svmr<-predict(svmModelRadial, newdata = test_data)
table(pred_svmr)

confusionMatrix(pred_svmr, test_data$risk_group, mode = "everything")

# Importancia de las variables 
varImp(svmModelRadial)

# - Probabilidades de clase 
svmModelRadial$pred$High_risk #train

test_prob_svmr<-predict(svmModelRadial, newdata = test_data, type = "prob") #test

# - ROC

roc_svmr<-roc(
    response = test_data$risk_group, #Clases reales
    predictor = test_prob_svmr$High_risk, #Probabilidad de ser high risk
    levels = c("Low_risk", "High_risk")
)

#AUC
auc(roc_svmr)
ci.auc(roc_svmr)


## - RANDOM FOREST ----

set.seed(123)
rfModel <- train(risk_group ~ .,
                 data = train_data,
                 method = "rf",
                 trControl = ctrl,
                 metric = "ROC",
                 tuneLenght = 3 #vemos 3 valores para que no tarde demasiado
)

rfModel

# - Predicciones sobre test
pred_rf<-predict(rfModel, newdata = test_data)
table(pred_rf)

confusionMatrix(pred_rf, test_data$risk_group, mode = "everything")

# - Probabilidades de clase 
rfModel$pred$High_risk #train

test_prob_rf<-predict(rfModel, newdata = test_data, type = "prob") #test

# - ROC

roc_rf<-roc(
  response = test_data$risk_group, #Clases reales
  predictor = test_prob_rf$High_risk, #Probabilidad de ser high risk
  levels = c("Low_risk", "High_risk")
)

#AUC
auc(roc_rf)
ci.auc(roc_rf)



# - CURVAS ROC PARA LOS 3 MODELOS

png(filename = "plots/ROC_multimodal.png", #para guardar
    width = 3000,     #en píxeles 
    height = 2400,    
    res = 300)

plot(roc_svml,
     col = "steelblue",
     lwd = 2,
     cex.lab=1.8)

plot(roc_svmr,
     col = "firebrick",
     lwd = 2,
     add = TRUE)

plot(roc_rf,
     col = "seagreen",
     lwd = 2,
     add = TRUE)

legend("bottomright",
       legend = c(
         paste0("SVM lineal (AUC=", round(auc(roc_svml),3),")"),
         paste0("SVM radial (AUC=", round(auc(roc_svmr),3),")"),
         paste0("Random Forest (AUC=", round(auc(roc_rf),3),")")
       ),
       col = c("steelblue","firebrick","seagreen"),
       lwd = 2,
       cex = 1.5)
dev.off()


