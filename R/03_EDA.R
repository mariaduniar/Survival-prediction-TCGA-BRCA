###############################################################################
#######                   ANÁLISIS EXPLORATORIOS                       ########
##############################################################################

rm(list = ls())

setwd("C:/Users/Usuario/Desktop/Bioinformatica/TFM")

library(tidyverse)
library(ggplot2)
library(stats) #para PCA
library(uwot) #para UMAP
library(factoextra) #representaciones
library(pheatmap) #heatmap
library(maftools) #visualizacion mutaciones
library(RColorBrewer)


## - CARGA DE DATOS ----

# - Matriz de expresión normalizada (VST) n=965 pacientes comunes en los 3 datasets
counts<-read.csv("count_matrix_norm_common.csv", row.names = 1) 
#nombres de fila como ID del paciente (no ID de la muestra)
rownames(counts)<- substr(rownames(counts), 1, 12)

# - Metadatos de expresión
metadata<-readRDS("expression_metadata.rds")
rownames(metadata) <- substr(rownames(metadata), 1, 12)

#- filtramos pacientes comunes
metadata<-metadata[rownames(metadata) %in% rownames(counts), ]

all(rownames(metadata)==rownames(counts)) #comprobamos alineación

# - Datos clínicos 
clin_data<-readRDS("clin_data_common.rds") 

# - Datos de mutaciones 
mut_data<-readRDS("mut_data_common.rds")


#############       1. EXPLORACIÓN DATOS CLÍNICOS    ######################


### - Estado vital ----


table(clin_data$vital_status) #832 pacientes vivos y 133 fallecidos

table(clin_data$vital_status, metadata$paper_BRCA_Subtype_PAM50)

# Diagrama de barras con estado vital
ggplot(clin_data, aes(x = vital_status, fill = vital_status)) +
  geom_bar(width = 0.6) +
  geom_text(
    stat = "count",
    aes(label = after_stat(count)),
    vjust = -0.3) +
  labs(x = "Estado vital",
    y = "Número de pacientes") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title = element_text(face="bold", color="gray20"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())+
  scale_fill_manual(values = c("aquamarine3", "indianred"))

ggsave("plots/barplot_vital_status.png", dpi = 300)

ggplot(clin_data, aes(x="", fill = vital_status)) +
  geom_bar() +
  coord_polar(theta = "y") +
  geom_text(
    stat="count",
    aes(label = after_stat(count)),
    position = position_stack(vjust = 0.5))+
  theme_void()+
  scale_fill_manual(values = c("aquamarine3", "indianred"))+
  guides(fill = guide_legend(title = "Estado vital"))

ggsave("plots/sector_vital_status.png", dpi = 300)

### - Subtipo molecular ----

metadata$paper_BRCA_Subtype_PAM50[is.na(metadata$paper_BRCA_Subtype_PAM50)]<-"Unknown"


ggplot(metadata, aes(x = paper_BRCA_Subtype_PAM50, fill = paper_BRCA_Subtype_PAM50)) +
  geom_bar(width = 0.6) +
  geom_text(
    stat = "count",
    aes(label = after_stat(count)),
    vjust = -0.3) +
  labs(x = "Subtipo molecular",
    y = "Nº pacientes") +
  scale_fill_manual(values=c("indianred2", "darkgoldenrod1", "aquamarine3", "palegreen",  "orchid1", "azure3"))+
  theme_minimal() +
  theme(legend.position = "none",
        axis.title = element_text(face="bold", color="gray20"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

ggsave("plots/barplot_subtype.png", dpi = 300)


### - Days to last follow up (seguimiento) ----

summary(clin_data$days_to_last_follow_up)
 #hay un valor negativo

#      -  Distribución general

ggplot(clin_data, aes(x=days_to_last_follow_up))+
  geom_histogram(fill="aquamarine4")+
  labs(x= "Días hasta último seguimiento",
       y= "Nº pacientes")+
  theme_minimal() +
  theme(axis.title = element_text(face="bold", color = "gray20"))


#     - Pacientes vivos 
clin_data_alive<-filter(clin_data, vital_status=="Alive")

summary(clin_data_alive$days_to_last_follow_up)

ggplot(clin_data_alive, aes(x=days_to_last_follow_up))+
  geom_histogram(fill="aquamarine4")+
  labs(x= "Días hasta último seguimiento",
       y= "Nº pacientes")+
  theme_minimal() +
  theme(axis.title = element_text(face="bold", color = "gray20"))

#        - En años
ggplot(clin_data_alive, aes(x=days_to_last_follow_up/365))+
  geom_histogram(fill="aquamarine4")+
  labs(x= "Tiempo seguimiento (años)",
       y= "Nº pacientes",
       title= "Pacientes vivos")+
  geom_vline(      # añadimos umbral 3 años
    xintercept = 3,
    linetype = "dashed",
    linewidth = 1,
    color = "gold"
  ) +
  annotate(
    "text",
    x = 3,
    y = 175,
    label = "3 años",
    color="darkgoldenrod3",
    fontface="bold",
    hjust=-0.2)+
  theme_minimal() +
  theme(axis.title = element_text(face="bold", color = "gray20"),
        plot.title = element_text(face="bold", hjust=0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

ggsave(file="plots/hist_followup_alive.png", dpi=300)

#En el histograma se observa que hay valores cercanos a 0 
table(
  clin_data$vital_status == "Alive",
  clin_data$days_to_last_follow_up == 0
)
  #hay 16 pacientes vivos sin seguimiento (días hasta último seguimiento=0)
  #hay 2 pacientes fallecidos con días=0

clin_data[clin_data$days_to_last_follow_up == 0,
          c("submitter_id", "vital_status", "days_to_last_follow_up", "days_to_death")]


### - Days to death (pacientes fallecidos) ----

summary(clin_data$days_to_death) #hay 833 NA (debería ser 832)

table(
  clin_data$vital_status == "Dead",
  is.na(clin_data$days_to_death)
)

#Nuevo df solo con los pacientes fallecidos
clin_data_dead<-filter(clin_data, vital_status=="Dead")

## Para los pacientes fallecidos, el día de ultimo seguimiento coincide con los 
# días hasta fallecimiento. Pero he visto que hay un NA en una columna de days_to death

clin_data_dead$days_to_last_follow_up==clin_data_dead$days_to_death 
#todo TRUE menos posición 66

#Sin embargo
clin_data_dead[66, "days_to_death"] # NA
clin_data_dead[66, "days_to_last_follow_up"] #da 26

#He decidido imputar ese dato con el valor de days to last follow up 
idx <- clin_data$vital_status == "Dead" &
  is.na(clin_data$days_to_death)

clin_data$days_to_death[idx] <-
  clin_data$days_to_last_follow_up[idx]

#Comprobamos
sum(clin_data$vital_status == "Dead" &
      is.na(clin_data$days_to_death))

#   - Distribución de days to death (sólo pacientes fallecidos)
clin_data_dead<-filter(clin_data, vital_status=="Dead") #df solo con dead

summary(clin_data_dead$days_to_death)

ggplot(clin_data_dead, aes(x=days_to_death))+
  geom_histogram(fill="indianred")+
  labs(x= "Días hasta fallecimiento",
       y= "Nº pacientes")+
  theme_minimal()+
  theme(axis.title = element_text(face="bold", color = "gray20"))

#      En años 
ggplot(clin_data_dead, aes(x=days_to_death/365))+
  geom_histogram(fill="indianred3")+
  labs(x= "Años hasta fallecimiento",
       y= "Nº pacientes",
       title= "Pacientes fallecidos")+
  geom_vline(      # añadimos umbral 3 años
    xintercept = 3,
    linetype = "dashed",
    linewidth = 1,
    color = "gold"
  ) +
  annotate(
    "text",
    x = 3,
    y = 20,
    label = "3 años",
    color="darkgoldenrod3",
    fontface="bold",
    hjust=-0.2)+
  theme_minimal()+
  theme(axis.title = element_text(face="bold", color = "gray20"),
        plot.title = element_text(face="bold", hjust=0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

ggsave(file="plots/hist_followup_dead.png", dpi=300)


##### - Seguimiento según estado vital ----

ggplot(clin_data, aes(x=vital_status , y=days_to_last_follow_up,))+
  geom_boxplot(alpha=0.7, aes(fill= vital_status))+
  geom_jitter(size=0.9, aes(color=vital_status))+
  labs(x="Estado vital",
       y="Días hasta último seguimiento / fallecimiento") +
  theme_minimal()+
  theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
        axis.title = element_text(face = "bold", color= "gray20"))+
  scale_fill_manual(values=c("aquamarine3","indianred2"))+
  scale_color_manual(values=c("aquamarine4","indianred3"))+
  theme(legend.position = "none")

#        En años 
ggplot(clin_data, aes(x=vital_status , y=days_to_last_follow_up/365))+
  geom_boxplot(alpha=0.7, aes(fill= vital_status))+
  geom_jitter(size=0.9, aes(color=vital_status))+
  labs(x="Estado vital",
       y="Tiempo seguimiento (años)") +
  theme_minimal()+
  theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
        axis.title = element_text(face = "bold", color= "gray20"))+
  scale_fill_manual(values=c("aquamarine3","indianred2"))+
  scale_color_manual(values=c("aquamarine4","indianred3"))+
  theme(legend.position = "none")

ggsave("plots/boxplot_follow_up.png", dpi=300)


summary(clin_data$days_to_last_follow_up/365)
summary(clin_data_alive$days_to_last_follow_up/365)

summary(clin_data$days_to_death/365)


### - Edad ----

summary(clin_data$age_at_index)

ggplot(clin_data, aes(x=age_at_index))+
  geom_histogram(fill="palegreen3", bins = 20)+
  labs(x= "Edad al diagnóstico",
       y= "Nº pacientes")+
  theme_minimal()+
  theme(axis.title = element_text(face="bold", color = "gray20"))

ggsave("plots/hist_age.png", dpi=300)

#   Edad según estado vital
ggplot(clin_data, aes(x=vital_status , y=age_at_index)) +
  geom_boxplot(alpha=0.7, aes(fill= vital_status)) +
  geom_jitter(size=0.9, aes(color=vital_status)) +
  labs(x="Estado vital",
       y="Edad en diagnóstico") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
        axis.title = element_text(face = "bold", color= "gray20")) +
  scale_fill_manual(values=c("aquamarine3","indianred2")) +
  scale_color_manual(values=c("aquamarine4","indianred3")) +
  theme(legend.position = "none")

ggsave("plots/boxplot_age.png", dpi=300)


### Estadio ----

table(clin_data$ajcc_pathologic_stage)
table(metadata$paper_pathologic_stage)

# Convertir los textos "NA" a NA reales
metadata$paper_pathologic_stage[metadata$paper_pathologic_stage == "NA"] <- NA
sum(is.na(metadata$paper_pathologic_stage))


metadata$paper_pathologic_stage[is.na(metadata$paper_pathologic_stage)]<-"Unknown"
table(metadata$paper_pathologic_stage)

ggplot(metadata, aes(x = paper_pathologic_stage, fill = paper_pathologic_stage)) +
  geom_bar(width = 0.6) +
  geom_text(
    stat = "count",
    aes(label = after_stat(count)),
    vjust = -0.3) +
  labs(x = "Estadio",
       y = "Nº pacientes") +
  scale_fill_manual(values=c("lavenderblush2", "rosybrown2", "indianred", "indianred4", "azure3"))+
  theme_minimal() +
  theme(legend.position = "none",
        axis.title = element_text(face="bold", color="gray20"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

ggsave("plots/barplot_stage.png", dpi = 300)



### - Supervivencia a 3 años (1095 días) ----

sum(clin_data_dead$days_to_death<1095) #64/133 pacientes fallecidos mueren antes de 3 años

sum(clin_data_alive$days_to_last_follow_up>=1095) #311/832 pacientes vivos tienen seguimiento pasados 3 años


#### - Nueva variable binaria High Risk/Low Risk basada en la supervivencia a 3 años -----

clin_data<- clin_data %>%
  mutate(risk_group = case_when(
    vital_status=="Dead" & days_to_death < 1095 ~ "High risk",
    days_to_last_follow_up >= 1095 ~ "Low risk",
    TRUE ~ "Excluded"
  ))

table(clin_data$risk_group)

# - Nuevo DF clínico sin los pacientes excluidos

clin_data_surv<-clin_data %>% filter(risk_group != "Excluded")

table(clin_data_surv$vital_status)
table(clin_data_surv$risk_group)

# - Filtramos matriz de conteo y metadatos 

surv_ids<-clin_data_surv$submitter_id

counts_surv<-counts[surv_ids, ] #filtramos
metadata_surv<-metadata[surv_ids, ]

write.csv(counts_surv, file = "counts_survival.csv") #guardamos 
saveRDS(counts_surv, file="counts_survival.rds")
saveRDS(metadata_surv, "metadata_survival.rds")
saveRDS(clin_data_surv, file= "clinical_survival.rds")

# - Añadimos nombre de fila a los datos clínicos y comprobamos que todos los DF están alineados 
rownames(clin_data_surv)<-clin_data_surv$submitter_id

table(rownames(clin_data_surv)==rownames(counts_surv))
table(rownames(counts_surv)==rownames(metadata_surv))
table(rownames(metadata_surv)==rownames(clin_data_surv))


#### Representaciones ----

metadata_surv$risk_group<-clin_data_surv$risk_group

# - Supervivencia y subtipo molecular 
ggplot(clin_data_surv, aes(x =risk_group, fill = metadata_surv$paper_BRCA_Subtype_PAM50)) +
  geom_bar(width = 0.9, position = "dodge") +
  geom_text(
    stat = "count",
    aes(label = after_stat(count)),
    position = position_dodge(width = 0.9), vjust=-0.5) +
  labs(x = "Grupo de riesgo",
       y = "Número de pacientes",
       fill = "Subtipo molecular (PAM50)") +
  theme_minimal() +
  theme(legend.position = "right",
        axis.title = element_text(face="bold", color="gray20"),
        legend.title = element_text(face="bold", color="gray20"))+
  scale_fill_manual(values=c("indianred2", "darkgoldenrod1", "aquamarine3", "palegreen",  "orchid1", "azure3"))

ggsave("plots/barplot_subtype_riskgroup_v2.png", width = 7, height = 4.5, dpi=300)

metadata_surv %>%
  filter(paper_BRCA_Subtype_PAM50 != "Unknown") %>%
  ggplot(aes(x =paper_BRCA_Subtype_PAM50, fill = risk_group)) +
  geom_bar(width = 0.9, position = "dodge") +
  geom_text(
    stat = "count",
    aes(label = after_stat(count)),
    position =position_dodge(width = 0.9), vjust=-0.5) +
  labs(x = "Subtipo molecular (PAM50)",
       y = "Nº pacientes",
       fill = "Grupo de riesgo") +
  theme_minimal() +
  theme(legend.position = "right",
        axis.title = element_text(face="bold", color="gray20"),
        legend.title = element_text(face="bold", color="gray20"))+
 scale_fill_manual(values=c("indianred3", "aquamarine3"))

ggsave("plots/barplot_subtype_riskgroup_v4.png", width = 7, height = 4.5, dpi=300)


# - Supervivencia y estado vital 
ggplot(clin_data_surv, aes(x = risk_group, fill = vital_status)) +
  geom_bar(width = 0.6) +
  geom_text(
    stat = "count",
    aes(label = after_stat(count)),
    position = position_stack(vjust = 0.5)) +
  labs(x = "Grupo de riesgo (supervivencia a 3 años)",
       y = "Nº pacientes",
       fill = "Estado vital") +
  theme_minimal() +
  theme(legend.position = "right",
        axis.title = element_text(face="bold", color="gray20"))+
  scale_fill_manual(values = c("aquamarine3", "indianred3"))

# - Supervivencia y edad 

summary(clin_data_surv$age_at_index)

#   Distribución edad según grupo riesgo 

  #summary de la edad por grupo
aggregate(age_at_index ~ risk_group, data = clin_data_surv, FUN = summary)

ggplot(clin_data_surv, aes(x=risk_group , y=age_at_index)) +
  geom_boxplot(alpha=0.7, aes(fill= risk_group)) +
  geom_jitter(size=0.9, aes(color=risk_group)) +
  labs(x="Grupo de riesgo",
       y="Edad al diagnóstico") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
        axis.title = element_text(face = "bold", color= "gray20")) +
  scale_fill_manual(values=c("indianred2","aquamarine3")) +
  scale_color_manual(values=c("indianred3", "aquamarine4")) +
  theme(legend.position = "none")

ggsave("plots/boxplot_age_riskgroup.png", dpi=300)

# - Supervivencia y estadio 

metadata_surv %>%
  filter(paper_pathologic_stage != "Unknown") %>%
  ggplot(aes(x =paper_pathologic_stage, fill = risk_group)) +
  geom_bar(width = 0.9, position = "dodge") +
  geom_text(
    stat = "count",
    aes(label = after_stat(count)),
    position =position_dodge(width = 0.9), vjust=-0.5) +
  labs(x = "Estadio",
       y = "Nº pacientes",
       fill = "Grupo de riesgo") +
  theme_minimal() +
  theme(legend.position = "right",
        axis.title = element_text(face="bold", color="gray20"),
        legend.title = element_text(face="bold", color="gray20"))+
  scale_fill_manual(values=c("indianred3", "aquamarine3"))

ggsave("plots/barplot_stage_riskgroup.png", height=4, dpi=300)


###############    2. EXPLORACIÓN DATOS EXPRESIÓN    ######################

### -  Selección genes más variables ----

# - Cálculo de la varianza de cada gen
gene_var<-apply(counts_surv, MARGIN=2, FUN=var)

# - Ver frecuencia de varianzas
hist(gene_var,breaks = 100) 

# - Nos quedamos con los más variables
top_n<-1000
top_genes<-names(sort(gene_var, decreasing = TRUE))[1:top_n]
counts_top_1000<-counts_surv[, top_genes] #matriz de conteo reducida
dim(counts_top_1000)


### - PCA ----

pca_result<-prcomp(counts_top_1000, center=TRUE, scale=FALSE) 
#los datos ya están normalizados, pero centramos a 0  

# - Varianza explicada por cada componente 

var_explicada<-pca_result$sdev^2/sum(pca_result$sdev^2)
cumsum(var_explicada)

  # Scree plot 

fviz_eig(pca_result, addlabels = TRUE)

# - Contribución de las variables

fviz_pca_var(
  pca_result,    
  select.var = list(contrib=20), #representamos las 20 variables mas contribuyentes a los PC
  col.var = "cos2",             # coloreamos variables según cos² (calidad representación)
  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
  repel = TRUE
)
ggsave("plots/PCA_contrib_top20var.png", dpi=300)


fviz_contrib(
  pca_result,
  choice = "var",
  axes = 1, # PC1
  top = 20  # número de variables mostradas
)
ggsave("plots/PCA_contrib_to_DIM1.png", dpi=300)

fviz_contrib(
  pca_result,
  choice = "var",
  axes = 2, # PC1
  top = 20  # número de variables mostradas
)

# - Representación gráfica

# Obtenemos las coordenadas de nuestras observaciones en el nuevo espacio de componentes
pca_x<-pca_result$x

# Etiquetas de los ejes del gráfico
x_label <- paste0(paste('PC1', round(var_explicada[1] * 100, 2)), '%')
y_label <- paste0(paste('PC2', round(var_explicada[2] * 100, 2)), '%')

#la función paste0 une dos elementos de texto sin separador (el %)
#la función paste une dos elementos de texto con separador, en este caso PC1 y la varianza explicada 

# Representación gráfica de las primeras dos componentes principales respecto a los datos

ggplot(pca_x, aes(x=PC1, y=PC2, color=metadata_surv$paper_BRCA_Subtype_PAM50)) + #color= colorea según variable categorica, en este caso por subtipo molecular
  geom_point(size=2) +
  labs(title='PCA', x=x_label, y=y_label, color="Subtipo Molecular (PAM50)") + #titulo del grafico y nombre de los ejes. 
  theme_classic() + 
  #a partir de aquí, elementos gráficos específicos dentro de theme()
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        plot.title = element_text(hjust = 0.5, face="bold", size=15),
        axis.title = element_text(face="bold", color="gray20"),
        legend.title = element_text(face="bold", color="gray20"))+
  scale_color_manual(values=c("indianred2", "darkgoldenrod1", "aquamarine3", "palegreen3",  "orchid1", "azure3"))

ggsave("plots/PCA.png", dpi=300)

#Según grupo de riesgo 
ggplot(pca_x, aes(x=PC1, y=PC2, color=clin_data_surv$risk_group)) + #color= colorea según variable categorica, en este caso por subtipo molecular
  geom_point(size=2) +
  labs(title='PCA', x=x_label, y=y_label, color="Grupo de riesgo") + #titulo del grafico y nombre de los ejes. 
  theme_classic() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        plot.title = element_text(hjust = 0.5, face="bold", size=15),
        axis.title = element_text(face="bold", color="gray20"),
        legend.title = element_text(face="bold", color="gray20"))+
  scale_color_manual(values=c("indianred2", "aquamarine3"))

ggsave("plots/PCA_riskgroup.png", dpi=300)

#   Con el paquete factoextra
fviz_pca_ind(
  pca_result,
  geom.ind = "point",
  col.ind = metadata_surv$paper_BRCA_Subtype_PAM50,
  addEllipses = TRUE
)

fviz_pca_ind(
  pca_result,
  geom.ind = "point",
  col.ind = clin_data_surv$risk_group,
  addEllipses = TRUE
) #no hay separación clara

#La fuente de variabilidad global viene más determinada por los subtipos moleculares


### - UMAP ----

#UMAP tiene muchas opciones para ajustar hiperparámetros. 
# La forma en la que se ha conseguido separar mejor los grupos ha sido con un nº vecinos bajo
# y min dist también bajo para que los grupos se compacten
umap_result <- umap(counts_top_1000, n_neighbors=5, 
                     n_components = 2, min_dist =0.001,  local_connectivity=1,
                     ret_model = TRUE, verbose = TRUE)

umap_df<-data.frame(umap_result$embedding) #data frame con las coordenadas

#Representtamos según subtipo molecular
ggplot(umap_df, aes(x=X1, y=X2, color=metadata_surv$paper_BRCA_Subtype_PAM50))+
  geom_point(size=3)+
  #scale_color_manual(values=c("indianred2", "aquamarine3", "palegreen", "darkgoldenrod1", "orchid1"))+
  labs(title='UMAP', x="Dim 1", y="Dim 2", color="Subtipo molecular (PAM50)") +
  theme_classic() +
  theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold", size = 15),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face="bold", color="gray20"))+
  scale_color_manual(values=c("indianred2", "darkgoldenrod1", "aquamarine3", "palegreen3",  "orchid1", "azure3"))
  
ggsave("plots/UMAP.png", dpi=300)

#Representamos según grupo riesgo 
ggplot(umap_df, aes(x=X1, y=X2, color=clin_data_surv$risk_group))+
  geom_point(size=3)+
  #scale_color_manual(values=c("indianred2", "aquamarine3", "palegreen", "darkgoldenrod1", "orchid1"))+
  labs(title='UMAP', x="Dim 1", y="Dim 2", color="Grupo de riesgo") +
  theme_classic() +
  theme(panel.grid.major=element_line (color="gray90"),panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="gray98"),
        plot.title = element_text(hjust = 0.5, face="bold", size = 15),
        axis.title = element_text(face = "bold") 
  )


### - CLUSTERING ----

#### - Clustering no jerárquico (k-means)

#escalamos los datos (z-score), impotante porque el k-means calcula distancias
counts_top_scaled<-scale(counts_top_1000)

# Determinar nº óptimo clusteres
fviz_nbclust(counts_top_scaled, kmeans, method = "wss") +
  ggtitle("Optimal number of clusters", subtitle = "") +
  theme_classic()

set.seed(1999)


# k means 
km<-kmeans(counts_top_scaled, centers=3, iter.max = 100, nstart=25)

# Visualización 
fviz_cluster(km, counts_top_scaled, xlab = '', ylab = '', geom="point") + #se escribe el resultado del clustering y los datos utilizados para el mismo
  ggtitle("Clustering no jerárquico con kmeans", subtitle = "") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, margin = margin(b = -10)))

table(km$cluster)

### - HEATMAP ----
## Con clustering jerárquico 

counts_top_var_t<-t(counts_top_1000) #trasponer para poner filas=genes y col=pacientes

annot_col<-data.frame("Subtipo molecular (PAM50)" = metadata_surv$paper_BRCA_Subtype_PAM50,
                      check.names = FALSE)
row.names(annot_col)<-rownames(metadata_surv)

all(colnames(counts_top_var_t)==rownames(annot_col)) #comprobamos alineacion

annot_colors<- list(
  "Subtipo molecular (PAM50)" = c(Basal = "indianred2",
                                  Her2 = "darkgoldenrod1",
                                  LumA = "aquamarine3",
                                  LumB = "palegreen3", 
                                  Normal = "orchid1",
                                  Unknown ="azure3")
)

heatmap<- pheatmap(counts_top_var_t,
                   clustering_method = "ward.D2",
                   color=colorRampPalette(c("deepskyblue4", "white", "firebrick2"))(100),
                   show_colnames = FALSE,
                   show_rownames = FALSE,
                   main="Heatmap de la expresión génica",
                   scale="row", 
                   #cutree_rows = 3 #separar filas (genes) según clusters formados
                   cutree_cols = 5, #separar pacientes según clusters formados
                   annotation_col =  annot_col,
                   annotation_names_col = FALSE,
                   annotation_colors = annot_colors, #color para las anotaciones
                   filename = "plots/heatmap_top_1000genes.png", #guardar
                   width = 7,
                   height = 5,
                   res=300
)

# - Segun risk group (NO resulta muy INFORMATIVO)

annot_col_riskgroup<-data.frame("Grupo de riesgo (supervivencia 3 años)" = clin_data_surv$risk_group,
                                           check.names = FALSE)
row.names(annot_col_riskgroup)<-rownames(clin_data_surv)

heatmap2<- pheatmap(counts_top_var_t,
                   clustering_method = "ward.D2",
                   color=colorRampPalette(c("deepskyblue4", "white", "firebrick2"))(100),
                   show_colnames = FALSE,
                   show_rownames = FALSE,
                   main="Heatmap de la expresión génica",
                   scale="row", 
                   #cutree_rows = 3 #separar filas (genes) según clusters formados
                   #cutree_cols = 2, #separar pacientes según clusters formados
                   annotation_col =  annot_col_riskgroup,
                   annotation_names_col = FALSE,
                   #annotation_colors =  #color para las anotaciones
                   filename = "heatmap_top_1000genes_surv.png", #guardar
                   width = 7,
                   height = 5,
                   res=300
)


#############       3. EXPLORACIÓN DATOS MUTACIONES    ######################

### - Filtrado ----

# - Filtramos pacientes incluidos en la clasificación riesgo (n=444) 

mut_surv<- mut_data %>%
  filter(sample_id %in% metadata_surv$sample) 

# - Comprobamos 
length(unique(mut_surv$Tumor_Sample_Barcode)) #444

# - Guardamos 
saveRDS(mut_surv, file = "mutations_survival.rds")

write.table(
  mut_surv,
  file = "mutations_surv.maf",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

### - Exploramos clasificación 

table(mut_surv$Variant_Type) #tipo de mutacion

table(mut_surv$Variant_Classification) #tipo de efecto en la proteina

### - VISUALIZACIONES  ----

#Convertimos a clase MAF 
 #para poder utilizar las funciones de maftools (https://bioconductor.org/packages/release/bioc/vignettes/maftools/inst/doc/maftools.html#7_Visualization)

maf <- read.maf(maf = mut_surv)
   #-Silent variants: 9185 
   #--Possible FLAGS among top ten genes:
   #  TTN
   # MUC16


#Graficamos
png(filename = "plots/summary_maf.png", #para guardar
    width = 3000,     #en píxeles 
    height = 2400,    
    res = 300)
plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()

png(filename = "plots/oncoplot.png", 
    width = 3000,     #en píxeles 
    height = 2400,    
    res = 300)
oncoplot(maf = maf, top = 10)
dev.off()


#curiodsidad: qué paciente es el que presenta una carga mutacional tan alta
mut_surv %>%
  count(Tumor_Sample_Barcode) %>%
  filter(n > 3000)



