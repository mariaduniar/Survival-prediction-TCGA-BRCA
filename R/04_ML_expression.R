################################################################################
############       ANÁLISIS EXPRESIÓN DIFERENCIAL Y MODELOS ML     #############
################################################################################

rm(list = ls())

setwd("C:/Users/Usuario/Desktop/Bioinformatica/TFM")

library(caret)
library(DESeq2)
library(dplyr)
library(EnhancedVolcano)
library(pheatmap)
library(pROC)

# - 1. CARGA DE DATOS ----

# - Datos de expresion, con los conteos CRUDOS de lecturas 
exp_data_raw<-readRDS("count_matrix_raw_survival.rds") #(ya filtrada por pacientes incluidos en la clasificación, n=444)
  
# - Datos clínicos 
clin_data<-readRDS("clinical_survival.rds")


# - 2. DIVISIÓN TRAIN/TEST ----

# - Añadimos columna risk group a los datos de expresión
exp_data_raw$risk_group <- clin_data$risk_group[
  match(rownames(exp_data_raw), clin_data$submitter_id)]

# - Índices para las observaciones
set.seed(1999)
trainIndex<-createDataPartition(exp_data_raw$risk_group, p=0.8, list=FALSE)

train_data<-exp_data_raw[trainIndex, ]
test_data<-exp_data_raw[-trainIndex, ]

# - Comprobamos proporciones 
table(exp_data_raw$risk_group)
table(train_data$risk_group)
table(test_data$risk_group)

# - GUARDAMOS LOS IDs de cada grupo, así me aseguro que train y test siempre son los mismos pacientes
train_ids <- rownames(train_data)
test_ids  <- rownames(test_data)

saveRDS(train_ids, "train_IDs.rds")
saveRDS(test_ids, "test_IDs.rds")

# - 3. ANÁLISIS EXPRESIÓN DIFERENCIAL ----
# - Sólo sobre el conjunto de Entrenamiento (train_data)

# - Cargamos metadatos (necesario para construir objeto dds)
metadata<-readRDS("metadata_survival.rds")

# -Filtramos metadatos por train_data y añadimos risk_group
metadata<-metadata %>%
  filter(rownames(metadata) %in% rownames(train_data))

all(rownames(train_data)==rownames(metadata))

metadata$risk_group <- train_data$risk_group[
  match(rownames(metadata), rownames(train_data))]

table(metadata$risk_group) #comprobamos
table(train_data$risk_group)

#Preparamos matriz de expresión para objeto dds
exp_matrix<-train_data[,-59428]
exp_matrix<-t(exp_matrix)

all(colnames(exp_matrix)==rownames(metadata))

## - 3.1. Creamos objeto dds ---- 

dds<-DESeqDataSetFromMatrix(
  countData = exp_matrix,
  colData = metadata, 
  design= ~ risk_group 
)

# Especificamos el grupo de referencia
dds$risk_group <- relevel(dds$risk_group, ref = "Low risk")


## - 3.2  Filtrado genes con baja expresión ----

# seguimos el mismo criterio que escogimos en el preprocesamiento 
smallestGroupSize <- 10 #no es obligatorio, hace referencia al mínimo de muestras por el que filtrar las counts
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

## - 3.3  Análisis expresión diferencial ----

dds<-DESeq(dds)

res<-results(dds, contrast = c("risk_group", "High risk", "Low risk"), alpha = 0.01) 
#alpha= umbral FDR

summary(res)

sum(res$padj < 0.01 & res$log2FoldChange >= 1, na.rm=TRUE)
#hay 108 genes sobreexpresados en high risk sobre low risk
sum(res$padj < 0.01 & res$log2FoldChange <= -1, na.rm=TRUE)
#en cambio hay 391 subexpresados 


## - 3.4 Representaciones gráficas ---- 

plotMA(res, ylim=c(-1,1))

### - Volcano plot 
res_df<-as.data.frame(res) #convertimos objeto res a data frame

EnhancedVolcano(res_df, 
                lab=rownames(res_df), 
                x='log2FoldChange', 
                y='pvalue', 
                pCutoff = 0.01, #umbral de significancia p value
                FCcutoff = 1, #umbral de log2Fold change
                labSize = 3, axisLabSize = 10)
ggsave("plots/DESeq2_train/volcano_logfc_1.png", dpi=300, width = 7, height = 7)


### - Heatmap 

# Normalizamos con VST
vsd<-vst(dds, blind=FALSE)

train_vst<-assay(vsd) #matriz de expresión normalizada 

# Obtenemos la matriz sólo con los genes diferencialmente expresados 
res_dif_genes<- res_df[(res$padj < 0.01 & res$log2FoldChange >= 1) | ((res$padj < 0.01 & res$log2FoldChange <= -1)), ]

train_vst_difgenes<-train_vst[rownames(res_dif_genes), ] #(recuerda los genes estan en las filas)

# Representamos heatmap 

annot<-data.frame("Grupo de riesgo (supervivencia 3 años)" = metadata$risk_group,
                  check.names = FALSE,
                  row.names = rownames(metadata))

table(colnames(train_vst_difgenes)==rownames(annot)) #comprobamos alineacion (nunca esta de mas)

annot_colors<- list(
  "Grupo de riesgo (supervivencia 3 años)" = c(`High risk`="orchid1", `Low risk`="cadetblue1")
)

heatmap<-pheatmap(train_vst_difgenes,
                  clustering_method = "ward.D2",
                  color=colorRampPalette(c("deepskyblue4", "white", "firebrick2"))(100),
                  show_colnames = FALSE,
                  show_rownames = FALSE,
                  main="Genes diferencialmente expresados (datos de entrenamiento)",
                  scale="row", 
                  #cutree_rows = 3 #separar filas (genes) según clusters formados
                  #cutree_cols = 5, #separar pacientes según clusters formados
                  annotation_col =  annot,
                  annotation_names_col = FALSE,
                  annotation_colors =  annot_colors,
                  filename = "plots/DESeq2_train/heatmap_DEG.png", #guardar
                  width = 8.5,
                  height = 5.5,
                  res=300
)

#Top 50 genes 

dif_genes_ordered<-res_dif_genes[order(res_dif_genes$padj, decreasing = FALSE), ]

train_vst_top50<-train_vst_difgenes[head(rownames(dif_genes_ordered), 50), ]

table(colnames(train_vst_top50)==rownames(annot))

pheatmap(train_vst_top50,
         clustering_method = "ward.D2",
         color=colorRampPalette(c("deepskyblue4", "white", "firebrick2"))(100),
         show_colnames = FALSE,
         show_rownames = TRUE,
         main="Top 50 genes diferencialmente expresados",
         scale="row", 
         cutree_rows = 4, #separar filas (genes) según clusters formados
         annotation_col =  annot,
         annotation_names_col = FALSE,
         annotation_colors = annot_colors,
         filename = "plots/DESeq2_train/heatmap_top50_DEG.png", #guardar
         width = 10,
         height = 7,
         res=300
)

#Sin clusterizar  pacientes 

annot <- annot[order(annot$`Grupo de riesgo (supervivencia 3 años)`), , drop = FALSE] #ordenamos los grupos
train_vst_top50 <- train_vst_top50[ ,rownames(annot)]

pheatmap(train_vst_top50,
         clustering_method = "ward.D2",
         color=colorRampPalette(c("deepskyblue4", "white", "firebrick2"))(100),
         cluster_cols = FALSE, #no clusterizamos columnas (pacientes)
         show_colnames = FALSE,
         show_rownames = TRUE,
         main="Top 50 genes diferencialmente expresados",
         scale="row", 
         cutree_rows = 4, #separar filas (genes) según clusters formados
         annotation_col =  annot,
         annotation_names_col = FALSE,
         annotation_colors = annot_colors,
         filename = "plots/DESeq2_train/heatmap_top50_DEG_noclust.png", #guardar
         width = 10,
         height = 6,
         res=300
)

## - 3.5 Guardamos datos ----

# - Resultados DESeq2
write.csv(res_df, file = "DESeq2_train/results_deseq.csv")

# - Lista genes DEG
write.csv(rownames(res_dif_genes), file="DESeq2_train/list_DEG.csv")


# - 4. MACHINE LEARNING ----

# Cargamos todo el conjunto de datos de expresión ya normalizado 
exp_data<-readRDS("counts_survival.rds")

# Seleccionamos genes diferencialmente expresados para filtrar la matriz
list_DEG<-rownames(res_dif_genes)

exp_data_DEG<-exp_data[,list_DEG]

#Añadimos risk group
exp_data_DEG$risk_group <- clin_data$risk_group[
  match(rownames(exp_data_DEG), clin_data$submitter_id)]

exp_data_DEG$risk_group<-as.factor(exp_data_DEG$risk_group)
str(exp_data_DEG$risk_group)

# Volvemos a dividir train/test data con los MISMOS ÍNDICES
train_data<-exp_data_DEG[trainIndex,]
test_data<-exp_data_DEG[-trainIndex,]

# Renombramos niveles porque algún algoritmo se ha quejado de los espacios 
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


## - 4.1. SVM ----

### - SVM lineal ----

# Creamos objeto para el parámetro trControl
ctrl <- trainControl(
  method = "repeatedcv", #validación cruzada 
  number = 10, #10 folds
  repeats = 5, #5 repeticiones
  classProbs = TRUE, #obtener probabilidades de clase
  savePredictions = "final", #guardar predicciones óptimas
  summaryFunction = twoClassSummary #calcular metricas de rendimiento de problemas de 2 clases, es decir ROC, sensibilidad y especificidad
)

# Entrenamos modelo
set.seed(123)
svmModelLineal <- train(risk_group ~.,
                        data = train_data,
                        method = "svmLinear",
                        trControl = ctrl,
                        preProcess = c("center", "scale"), #escalado de datos 
                        tuneGrid = expand.grid(C = c(0.01, 0.1, 1, 10)), #prueba 4 valores de C
                        metric="ROC" )

# - Evaluamos resultados
svmModelLineal
   #se ha escogido C=0,01 

# - Predicciones sobre el test con umbral normal (P=0,5)
pred_SVML <- predict(svmModelLineal, newdata = test_data )
table(pred_SVML)

confusionMatrix(pred_SVML, test_data$risk_group, positive = "High_risk", mode = "everything")
  #16% sensibilidad y 89% especificidad

head(svmModelLineal$pred) #predicciones de clase (probabilidades) para cada paciente obtenidas durante la VALIDACIÓN CRUZADA

# - creamos objeto ROC sobre las predicciones OOF de TRAIN 
roc_train_svml <- roc(
  response = svmModelLineal$pred$obs, #Clases reales
  predictor = svmModelLineal$pred$High_risk, #Probabilidad de ser high risk
  levels = c("Low_risk", "High_risk")
)

# - Vemos mejor umbral de decisión según Youden
best_svml <- coords(
  roc_train_svml,
  x = "best",
  best.method = "youden",
  ret = c("threshold", "sensitivity", "specificity"))

best_svml

plot(roc_train_svml,
     col = "steelblue",
     lwd = 2,
     print.auc = TRUE,
     print.thres=TRUE) #representar umbral decisión óptimo


# PREDICCIONES UTILIZANDO UMBRAL DE YOUDEN
# Probabilidades de clase sobre el test
test_prob_svm<- predict(svmModelLineal, newdata = test_data, type = "prob")

threshold <- best_svml$threshold #umbral de decision youden 

#Cambiamos las predicciones con este umbral y evaluamos modelo 
pred_svml_youden <- ifelse(
  test_prob_svm$High_risk >= threshold,
  "High_risk",
  "Low_risk"
)

pred_svml_youden <- factor(pred_svml_youden, levels = c("Low_risk","High_risk"))

confusionMatrix(pred_svml_youden, test_data$risk_group, positive = "High_risk", mode="everything")
   #41% sensibilidad y 69% especificidad


### - SVM con kernel gaussiano ----

# - Utilizamos el mismo objeto para trControl

# - Entrenamos modelo 
set.seed(123)
svmModelRadial <- train(risk_group ~ .,
                        data = train_data,
                        method = "svmRadial",
                        trControl = ctrl,
                        preProcess = c("center", "scale"),
                        metric = "ROC",
                        tuneLength = 4
)

svmModelRadial

# - Predicciones sobre el test con umbral normal
pred_SVMR <- predict(svmModelRadial, newdata = test_data )
table(pred_SVMR)

confusionMatrix(pred_SVMR, test_data$risk_group, positive = "High_risk", mode = "everything")
  #sensibilidad 25% y especificidad 90%

# - creamos objeto ROC sobre TRAIN 
roc_train_svmr <- roc(
  response = svmModelRadial$pred$obs, #Clases reales
  predictor = svmModelRadial$pred$High_risk, #Probabilidad de ser high risk
  levels = c("Low_risk", "High_risk")
)

# - Vemos mejor umbral de decisión según Youden
best_svmr <- coords(
  roc_train_svmr,
  x = "best",
  best.method = "youden",
  ret = c("threshold", "sensitivity", "specificity"))

best_svmr

plot(roc_train_svmr,
     col = "firebrick",
     lwd = 2,
     print.auc = TRUE,
     print.thres=TRUE) #representar umbral decisión óptimo


# PREDICCIONES UTILIZANDO UMBRAL DE YOUDEN

# Probabilidades de clase sobre el test
test_prob_svmr<- predict(svmModelRadial, newdata = test_data, type = "prob")

# Cambiamos las predicciones con el umbral youden y evaluamos modelo 
pred_svmr_youden <- ifelse(
  test_prob_svmr$High_risk >= best_svmr$threshold,
  "High_risk",
  "Low_risk"
)

pred_svmr_youden <- factor(pred_svmr_youden, levels = c("Low_risk","High_risk"))

confusionMatrix(pred_svmr_youden, test_data$risk_group, positive = "High_risk", mode = "everything")
    #41% sensibilidad y 72% especificidad


## - 4.2. RANDOM FOREST ----

modelLookup("rf")
set.seed(123)
rfModel <- train(risk_group ~ .,
                 data = train_data,
                 method = "rf",
                 trControl = ctrl,
                 metric = "ROC",
                 tuneGrid= expand.grid(mtry=c(5,25,50)) #vemos 3 valores para que no tarde demasiado
)
#No hacemos preprocesamiento pues RF no lo necesita

rfModel 

# - Predicciones sobre el test (con umbral normal)
pred_RF<-predict(rfModel, newdata = test_data)
table(pred_RF)

confusionMatrix(pred_RF, test_data$risk_group, positive = "High_risk", mode = "everything")
 #8% sensibilidad y 97% especificidad 


# - creamos objeto ROC sobre TRAIN 
roc_train_rf <- roc(
  response = rfModel$pred$obs, #Clases reales
  predictor = rfModel$pred$High_risk, #Probabilidad de ser high risk
  levels = c("Low_risk", "High_risk")
)

# - Vemos mejor umbral de decisión según Youden
best_rf <- coords(
  roc_train_rf,
  x = "best",
  best.method = "youden",
  ret = c("threshold", "sensitivity", "specificity"))

best_rf 

plot(roc_train_rf,
     col = "seagreen",
     lwd = 2,
     print.auc = TRUE,
     print.thres=TRUE) #representar umbral decisión óptimo


# PREDICCIONES UTILIZANDO UMBRAL DE YOUDEN

# Probabilidades de clase sobre el test
test_prob_rf<- predict(rfModel, newdata = test_data, type = "prob")

# Cambiamos las predicciones con el umbral youden y evaluamos modelo 
pred_rf_youden <- ifelse(
  test_prob_rf$High_risk >= best_rf$threshold,
  "High_risk",
  "Low_risk"
)

pred_rf_youden <- factor(pred_rf_youden, levels = c("Low_risk","High_risk"))

confusionMatrix(pred_rf_youden, test_data$risk_group, positive = "High_risk", mode = "everything")
     #75% sensibilidad y 61% especificidad



## - 4.3. CURVAS ROC (test) ----

#ya tenemos calculadas las probabilidades de clase sobre el test, así que solo
# hay que construir el objeto ROC para el test
# y representar las curvas ROC juntas 

#-ROC

roc_test_svml<-roc(
  response = test_data$risk_group, #Clases reales
  predictor = test_prob_svm$High_risk, #Probabilidad de ser high risk
  levels = c("Low_risk", "High_risk")
)

roc_test_svmr<-roc(
  response = test_data$risk_group, #Clases reales
  predictor = test_prob_svmr$High_risk, #Probabilidad de ser high risk
  levels = c("Low_risk", "High_risk")
)

roc_test_rf<-roc(
  response = test_data$risk_group, #Clases reales
  predictor = test_prob_rf$High_risk, #Probabilidad de ser high risk
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

png(filename = "plots/DESeq2_train/ROC_expression.png", #para guardar
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


