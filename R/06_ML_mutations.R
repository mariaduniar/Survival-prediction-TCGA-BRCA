################################################################################
############       ANÁLISIS COMPARATIVOS MAF Y TRANSFORMACIÓN      #############
################################################################################

rm(list = ls())

setwd("C:/Users/Usuario/Desktop/Bioinformatica/TFM")

library(maftools)
library(dplyr)
library(tidyr)
library(tibble)
library(caret)
library(msigdbr) #para construir vías a partir de genes
library(Matrix)
library(glmnet)
library(pROC)

# 1. CARGA DE DATOS ----

mut_data<-readRDS("mutations_survival.rds") #MAF

clin_data<-readRDS("clinical_survival.rds") #datos clínicos


# 2. Transformación a matriz binaria para ML ----  

#paciente x gen (1/0 mutado o no)

# - Eliminamos mutaciones silenciosas
mut_binary<-mut_data %>%
  filter(Variant_Classification!="Silent")

# - Eliminamos mutaciones clasificadas como benignas (PolyPhen)
mut_binary<- mut_binary %>% 
  filter(!grepl("^benign", PolyPhen))

# - Eliminamos mutaciones repetidas en el mismo gen por paciente para construir la matriz binaria
mut_binary<-mut_binary %>%
  distinct(patient_id, Hugo_Symbol)

# - Añadimos nueva columna con valor=1
mut_binary<-mut_binary %>%
  mutate(mutated=1)

# - Convertimos a formato ancho 
mut_matrix<-mut_binary %>%
  pivot_wider(
    names_from = Hugo_Symbol,
    values_from = mutated,
    values_fill = 0
  )
mut_matrix<-column_to_rownames(mut_matrix, var="patient_id")
  #hay 8961 genes 

write.csv(mut_matrix, "mut_data_binary.csv")

# - Vemos la frecuencia de los genes 
gene_freq<-colSums(mut_matrix)

summary(gene_freq)
hist(gene_freq, breaks = 150)
 #hay muchos que solo aparecen mutados en 1 muestra (freq=1)


# - 3. Preparamos para ML ----

#Añadimos risk group y convertimos a factor

mut_matrix$risk_group <- clin_data$risk_group[
  match(rownames(mut_matrix), clin_data$submitter_id)]

table(mut_matrix$risk_group)
mut_matrix$risk_group<-as.factor(mut_matrix$risk_group)

str(mut_matrix$risk_group)


# - 4. División Train/Test ----

set.seed(1999)

#Índices para las observaciones
trainIndex<-createDataPartition(mut_matrix$risk_group, p=0.8, list=FALSE)

train_data<-mut_matrix[trainIndex, ]
test_data<-mut_matrix[-trainIndex, ]

# - Comprobamos proporciones 
table(mut_matrix$risk_group)
table(train_data$risk_group)
table(test_data$risk_group)


## - Preprocesamiento SOLO TRAIN ----

# - Volvemos a ver frecuencia en train 
gene_freq_train<-colSums(train_data[,-8962])

summary(gene_freq_train)

sum(gene_freq_train==0) #hay 1745 genes que en el conjunto de entrenamiento no aparecen mutados en ninguna muestra

### - Filtramos mutaciones que aparecen en al menos 4 muestras (~ 1% de la cohorte) ----
genes_mut<-names(gene_freq_train[gene_freq_train>=4]) #lista de genes que retenemos (404)

train_data <- train_data[, c(genes_mut, "risk_group"), drop= FALSE]

### - Seleccionamos genes con LASSO ----

# - Entrenamos modelo 
set.seed(123)
lasso_model <- cv.glmnet(x=as.matrix(train_data[,-405]),
                         y=train_data$risk_group,
                         family = "binomial", alpha = 1)
#binomial: indica problema de clasificación binaria
#alpha: indica LASSO puro
#por defecto nfolds=10

#  Obtenemos coeficientes
coef_lasso <- coef(lasso_model, s = "lambda.min")

#  Seleccionamos paths lambda distinto de 0
genes_lasso <- rownames(coef_lasso)[coef_lasso[, 1] != 0]
genes_lasso

 ## No ha retenido ningún gen

# - 7. CAMBIAR GENES -> VÍAS ALTERADAS ----


# Lista de genes de nuestra matriz binaria (TODOS, no solo los que filtramos por frecuencia en train)

genes<-colnames(mut_matrix[,-8962])

#Vemos las colecciones de MSigBD
head(msigdbr_collections(),20)

#Elegimos KEGG (186 VÍAS)
pathways <- msigdbr(db_species = "HS",
                    species = "human",
                    collection = "C2",
                    subcollection = "CP:KEGG_LEGACY" 
)

# Vemos cuántos de mis genes están anotados en pathways
sum(genes %in% pathways$gene_symbol) #2552

## - CONSTRUIMOS MATRIZ [PACIENTES X VÍAS]

#Filtramos pathways con los genes que tenemos mutados 
pathway_genes <- pathways[
  pathways$gene_symbol %in% colnames(mut_matrix),
  c("gs_name", "gene_symbol")
]

length(unique(pathway_genes$gs_name)) #nuestros genes aparecen en las 186 vías
pathway_names <- unique(pathway_genes$gs_name) #nombres de las vías obtenidas

   #Vamos a explorar cuantos genes hay por vía 
genes_por_via <- table(pathway_genes$gs_name)
View(genes_por_via)


#construimos matriz gen x vía
gene_pathway <- sparseMatrix(
  i = match(pathway_genes$gene_symbol, genes), #coordenada genes
  j = match(pathway_genes$gs_name, pathway_names), #coordenadas pathways
  x = 1,
  dims = c(length(genes), length(pathway_names)),
  dimnames = list(genes, pathway_names)
)

# - Multiplicamos mut_matriz (paciente x gen) por matriz gene_pathway (gen x vías)

mut_matrix_genes<-as.matrix(mut_matrix[,-8962]) #matriz de mutaciones sin risk group

dim(mut_matrix_genes)
dim(gene_pathway) #el número de cols en la primera matriz(genes) debe ser igual al numero de observaciones en la otra 

# OBTENEMOS MATRIZ FINALMENTE:
mut_pathway <- mut_matrix_genes %*% gene_pathway
dim(mut_pathway) #Obtenemos 444 observaciones (pacientes) x 186 vías

mut_pathway<-as.matrix(mut_pathway)

# La convertimos a binaria (0/1 = ausencia/presencia de un gen mutado en esa vía)
mut_pathway[mut_pathway > 0] <- 1

# FILTRAMOS vías que tengan +10 genes para eliminar ruido 
# (lo comento porque no lo he hecho)
#mut_pathway_filt<-mut_pathway[, names(genes_por_via[genes_por_via >= 10])]
#dim(mut_pathway_filt)


## - 7.1 Preparamos nueva matrix para ML ----

mut_pathway<-as.data.frame(mut_pathway)

# Añadimos risk group 
mut_pathway$risk_group <- clin_data$risk_group[
  match(rownames(mut_pathway), clin_data$submitter_id)
]

table(mut_pathway$risk_group)
mut_pathway$risk_group<-as.factor(mut_pathway$risk_group)

str(mut_pathway$risk_group)
write.csv(mut_pathway, "mut_data_binary_pathways.csv")

# Dividimos Train/Test usando los IDs que tenemos guardados
train_ids<-readRDS("train_IDs.rds")
test_ids<-readRDS("test_IDs.rds")

train_data_path<-mut_pathway[train_ids, ]
test_data_path<-mut_pathway[test_ids, ]

#Proporciones
table(mut_pathway$risk_group)
table(train_data_path$risk_group)
table(test_data_path$risk_group)

## - 7.2 Preprocesamiento TRAIN DATA ----

# - Vemos frecuencia de los pathways en la cohorte train

path_freq<-colSums(train_data_path[, -187])
summary(path_freq)

hist(path_freq, breaks = 30)

# Renombramos niveles
train_data_path$risk_group <- factor(
  train_data_path$risk_group,
  levels = c("High risk", "Low risk"),
  labels = c("High_risk", "Low_risk")
)

test_data_path$risk_group <- factor(
  test_data_path$risk_group,
  levels = c("High risk", "Low risk"),
  labels = c("High_risk", "Low_risk")
)


##- 7.3 LASSO ----
# para seleccionar vías 

# - Entrenamos modelo LASSO
set.seed(1)
lasso_model <- cv.glmnet(x=as.matrix(train_data_path[,-187]),
                         y=train_data_path$risk_group,
                         family = "binomial", alpha = 1)
#binomial: indica problema de clasificación binaria
#alpha: indica LASSO puro

#  Obtenemos coeficientes
coef_lasso <- coef(lasso_model, s = "lambda.min")

#  Seleccionamos paths lambda distinto de 0
pathways_lasso <- rownames(coef_lasso)[coef_lasso[, 1] != 0]
pathways_lasso

#  Eliminamos intercepto 
pathways_lasso <- pathways_lasso[pathways_lasso != "(Intercept)"]
pathways_lasso #se han seleccionado 10 vías 

write.csv(pathways_lasso, "pathways_LASSO.csv")

# Filtramos train y test de nuevo 
train_data_path<-train_data_path[,c(pathways_lasso, "risk_group")]
test_data_path<-test_data_path[,c(pathways_lasso, "risk_group")]


## - 7.4. SVM LINEAL ----

# Parámetro trcontrol
ctrl <- trainControl(
  method = "repeatedcv", #validación cruzada 
  number = 10, #10 folds
  repeats = 5, #5 repeticiones
  classProbs = TRUE, #obtener probabilidades de clase
  savePredictions = "final", #guardar predicciones óptimas
  summaryFunction = twoClassSummary #calcular metricas de rendimiento de problemas de 2 clases, es decir ROC, sensibilidad y especificidad
)

set.seed(123)
svmModelLineal_path <- train(risk_group ~.,
                        data = train_data_path,
                        method = "svmLinear",
                        trControl = ctrl,
                        preProcess = NULL,
                        tuneGrid = expand.grid(C = c(0.01, 0.1, 1, 10)), #prueba 4 valores de C
                        metric="ROC", 
                        prob.model = TRUE) 

svmModelLineal_path

#Predicciones de clase sobre el test
pred_SVML_path <- predict(svmModelLineal_path, newdata = test_data_path )
table(pred_SVML_path)

confusionMatrix(pred_SVML_path, test_data_path$risk_group, positive = "High_risk", mode="everything")
#no se ha predcho ninguna observacion como high risk (sensibilidad 0)

# Probabilidades de clase sobre test
test_prob_svml_path<-predict(svmModelLineal_path, newdata = test_data_path, type = "prob")
View(test_prob_svml_path) #las probabilidades son prácticamente el mismo valor para todas las observaciones (pasa lo mismo que con geneS)
sd(test_prob_svml_path$High_risk) #la desviacion es muy baja
(svmModelLineal_path$pred) #igual para el train durante la cv


## 7.5. SVM KERNEL GAUSSIANO ----
set.seed(123)
svmModelRadial_path <- train(risk_group ~.,
                        data = train_data_path,
                        method = "svmRadial",
                        trControl = ctrl,
                        preProcess = NULL,
                        tuneLength = 4, #se deja que el modelo escoja los hiperparametros
                        prob.model = TRUE,
                        metric="ROC") 

svmModelRadial_path
plot(svmModelRadial_path)

# - Predicciones sobre el test 
pred_SVMR_path <- predict(svmModelRadial_path, newdata = test_data_path )
table(pred_SVMR_path)
confusionMatrix(pred_SVMR_path, test_data_path$risk_group, positive = "High_risk", mode="everything")


# - Probabilidades de clase 
test_prob_svmr_path<-predict(svmModelRadial_path, newdata = test_data_path, type = "prob")
View(test_prob_svmr_path)
summary(test_prob_svmr$High_risk) #las probabilidades tienen algo mas de variabilidad que en el modelo lineal

varImp(svmModelRadial_path)

## - 7.6 RF ----
set.seed(123)
rf_model_path <- train(risk_group ~ .,
                  data = train_data_path,
                  method = "rf",
                  trControl = ctrl,
                  preProcess = NULL,  
                  tuneLength = 3, #explora 3 valores distintos
                  importance = TRUE, #calcula la importancia de las variables
                  metric="ROC" )

## Resumen del modelo entrenado
rf_model_path
#mrty=2
plot(rf_model_path)

# - Predicciones sobre el test 
pred_rf_path<-predict(rf_model_path, newdata = test_data_path)
table(pred_rf_path)
confusionMatrix(pred_rf_path, test_data_path$risk_group, positive = "High_risk")
 #sensibilidad 0

test_prob_rf_path<-predict(rf_model_path, newdata = test_data_path, type="prob")

## 7.8. CURVAS ROC ----

roc_test_svml<-roc(
  response = test_data_path$risk_group, #Clases reales
  predictor = test_prob_svml_path$High_risk, #Probabilidad de ser high risk
  levels = c("Low_risk", "High_risk")
)

roc_test_svmr<-roc(
  response = test_data_path$risk_group, #Clases reales
  predictor = test_prob_svmr_path$High_risk, #Probabilidad de ser high risk
  levels = c("Low_risk", "High_risk")
)

roc_test_rf<-roc(
  response = test_data_path$risk_group, #Clases reales
  predictor = test_prob_rf_path$High_risk, #Probabilidad de ser high risk
  levels = c("Low_risk", "High_risk")
)

# - AUC con su IC 95%

auc(roc_test_svml) # SVM lineal
ci.auc(roc_test_svml)

auc(roc_test_svmr) # SVM radial
ci.auc(roc_test_svmr)

auc(roc_test_rf) # RF
ci.auc(roc_test_rf)

# - Representamos gráficamente

png(filename = "plots/ROC_mut_pathways.png", #para guardar
    width = 3000,     #en píxeles 
    height = 2400,    
    res = 300)

plot(roc_test_svml,
     col = "steelblue",
     lwd = 2,
     cex.lab=1.8)

plot(roc_test_svmr,
     col = "firebrick",
     lwd = 2,
     add = TRUE)

plot(roc_test_rf,
     col = "seagreen",
     lwd = 2,
     add = TRUE)

legend("bottomright",
       legend = c(
         paste0("SVM lineal (AUC=", round(auc(roc_test_svml),3),")"),
         paste0("SVM radial (AUC=", round(auc(roc_test_svmr),3),")"),
         paste0("Random Forest (AUC=", round(auc(roc_test_rf),3),")")
       ),
       col = c("steelblue","firebrick","seagreen"),
       lwd = 2,
       cex = 1.5)
dev.off()
