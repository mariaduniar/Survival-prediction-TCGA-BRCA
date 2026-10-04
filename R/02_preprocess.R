##############################################################################
################### PREPROCESAMIENTO DE DATOS ###############################
##############################################################################

rm(list = ls())

setwd("C:/Users/Usuario/Desktop/Bioinformatica/TFM")

library(SummarizedExperiment)
library(tibble) #para pasar rownames a columnas
library(dplyr)
library(DESeq2) #filtrado y normalización VST


## 1. - Cargamos datos guardados en objeto R ----

exp_data<-readRDS("expression.rds") #datos de expresión
clin_data<-readRDS("clinical.rds") #datos clínicos
mut_data<-readRDS("mutations.rds") #datos de mutaciones somáticas


## 2. - Preparamos matrices ----

### 2.1. - MATRIZ DE EXPRESIÓN ----

counts<-as.data.frame(assay(exp_data)) #matriz de conteo
annot<-as.data.frame(rowRanges(exp_data)) #matriz de anotaciones de los genes
metadata<-as.data.frame(colData(exp_data)) #metadatos


#### 2.1.1. - Asignamos los IDs a nombres de genes ----

# - Comprobamos que los nombres de fila (ID Ensembl) están en el mismo orden 
all(rownames(counts) == rownames(annot))

# - Cambiamos por el nombre del gen y agrupamos 
anyNA(annot$gene_name) #comprobamos que no hay NA
any(annot$gene_name=="") #comprobamos que no hay vacíos

counts_genes<-rowsum(counts, group = annot$gene_name)

  # *otra forma de hacerlo: counts_genes2<-aggregate(counts, by=list(annot$gene_name), FUN="sum")


#### 2.1.2. - Filtramos muestras según datos clínicos ----

count_matrix<-as.data.frame(t(counts_genes)) #trasponemos 
 #Ahora tenemos las observaciones como pacientes (1111) y las variables como genes (60660)

#Pasamos los nombres de fila a una nueva columna para poder manejarlos mejor
count_matrix<-rownames_to_column(count_matrix, var="barcode")

#En esta matriz algunos pacientes están duplicados (de una misma muestra se extraen varios viales A,B o C).
#Por ejemplo tenemos los barcodes: TCGA-A7-A0DC-01A-11R-A00Z-07 y TCGA-A7-A0DC-01B-04R-A22O-07 cuando 
#en los datos clínicos solo tenemos TCGA-A7-A0DC. (https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/#:~:text=The%20following%20figure%20illustrates%20how,Plate)

#Por tanto, vamos a quedarnos sólo con 1 muestra por paciente: 

# - Creamos una nueva columna donde extraemos los 12 primeros caracteres (ID de pacientes) con la función substring
count_matrix$patient_id<-substr(count_matrix$barcode, 1, 12)

# -La movemos a la 2nda posicion
count_matrix <- count_matrix %>% relocate(patient_id, .after = 1)

# - Eliminamos los duplicados
count_matrix <- count_matrix[!duplicated(count_matrix$patient_id), ]
  #Observamos que ahora tenemos 1095 pacientes

# - Comprobamos que no haya duplicados
anyDuplicated(count_matrix$barcode)

# - Guardamos la matriz de expresión sin duplicados 
write.csv(count_matrix, file = "count_matrix_raw.csv")

# - Finalmente, filtramos los pacientes que tengan datos clínicos 

clinical_ids<-clin_data$submitter_id #vector con los ids que aparecen en la matriz clinical

count_matrix_filt<-count_matrix[count_matrix$patient_id %in% clinical_ids, ]

#No se ha eliminado ninguna fila, por lo que todos los pacientes que tienen datos de expresión tienen datos clínicos
# (lo cual era un poco lógico ya que la matriz clínica tiene 1098 pacientes)
# Se puede comprobar con:
all(count_matrix$patient_id %in% clin_data$submitter_id)

#Eliminamos entonces la matriz filtrada pues es la misma
rm(count_matrix_filt)


#### 2.1.3.- Filtrado y normalización con DESeq2 ----

# - Necesitamos la matriz de metadatos para obtener un objeto DESeqDataSet. 
#Al igual que la matriz de expresión, contiene todas la muestras analizadas (1111), no de los pacientes, 
#por lo que habría que filtrarla como hemos hecho con la matriz de conteo

metadata<-metadata[metadata$barcode %in% count_matrix$barcode, ]

saveRDS(metadata, "expression_metadata.rds") #lo guardamos como objeto R ya que tiene columnas complejas (listas)


# - Finalmente, construimos el objeto DESeqDataSet

   # Hay que volver a trasponer la matriz de conteo ya que DESeq2 sólo acepta los barcodes como nombres de columna 
count_matrix_dds<-count_matrix[, -(1:2)] #quitamos columnas 1y 2
rownames(count_matrix_dds)<-count_matrix$barcode #añadimos barcodes como nombres de fila
count_data_dds <- t(count_matrix_dds) #trasponemos
dim(count_data_dds)
view(count_data_dds)

dds<-DESeqDataSetFromMatrix(
  countData = count_data_dds,
  colData = metadata, #matriz de metadatos
  design= ~ 1 #no ponemos ninguna condición todavía porque no se va a hacer análisis diferencial (https://rdrr.io/bioc/DESeq2/man/DESeqDataSet.html)
  ) 

# - Filtrado de genes con baja expresión 

smallestGroupSize <- 10 #no es obligatorio, hace referencia al mínimo de muestras por el que filtrar las counts.
# Un valor de 10 indica en la siguiente línea que retenemos los genes con valores counts mayor o igual a 10 en al menos 10 muestras 
#Normalmente se elige según el número de muestras que tengas en tus grupos experimentales (no es el caso ahora) 
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

dim(dds) #han pasado el filtro 34.440 genes


# - Normalización 

# Vamos a normalizar los datos con el método VST (variance Stabilizing Tranformation), más
# recomendado para tamaños muestrales grandes por la eficiencia computacional.  
# Calcula una transformación que estabiliza la varianza a lo largo de todo el rango de medias.

vsd <- vst(dds, blind=TRUE)

#El argumento blind se utiliza si se quiere que la normalización sea ajena al diseño experimental. 
#Nosotros con ~1 no especificamos diseño pero aún así vamos a poner TRUE.

# Extraemos la matriz filtrada y normalizada: 

vst_counts<-as.data.frame(t(assay(vsd))) #volvemos a trasponer para tener filas=pacientes 

write.csv(vst_counts, file = "count_matrix_norm.csv")
   

### 2.2 - MATRIZ DE MUTACIONES (MAF) ----

View(mut_data)

# - Vemos cuantas muestras hay 
length(unique(mut_data$Tumor_Sample_Barcode)) #(n=991)

# - Filtramos para que sólo haya mutaciones del tumor primario (codigo 01) 
mut_data <- mut_data %>%
  filter(substr(Tumor_Sample_Barcode, 14, 15) == "01")

# Extraemos los ID de pacientes
mut_data$patient_id <- substr(mut_data$Tumor_Sample_Barcode, 1, 12)

length(unique(mut_data$Tumor_Sample_Barcode)) #985 muestras
length(unique(mut_data$patient_id)) #968 pacientes 
#También hay pacientes duplicados

# - Filtramos mismas muestras que conservamos en la matriz de expresión

mut_data$sample_id <- substr(mut_data$Tumor_Sample_Barcode, 1, 16)

mut_data_filt<- mut_data %>%
  filter(sample_id %in% metadata$sample) 

length(unique(mut_data_filt$Tumor_Sample_Barcode)) #n=971 sigue habiendo duplicados porque el barcode no es el mismo

# - Vemos los pacientes duplicados
duplicados <- mut_data_filt %>%
  distinct(patient_id, Tumor_Sample_Barcode) %>%
  count(patient_id) %>%
  filter(n > 1)

duplicados

# Seleccionar un único barcode por paciente
barcode_keep <- mut_data_filt %>%
  distinct(patient_id, Tumor_Sample_Barcode) %>%
  group_by(patient_id) %>%
  slice(1) %>%
  ungroup()

# Conservar todas las mutaciones de esos barcodes
mut_data_filt <- mut_data_filt %>%
  semi_join(barcode_keep,
            by = c("patient_id", "Tumor_Sample_Barcode"))

# Comprobamos
length(unique(mut_data_filt$patient_id))
length(unique(mut_data_filt$Tumor_Sample_Barcode))



# - FILTRAR IDS COMUNES ----

#Vector con los ID de la matriz de mutaciones
mutation_ids <- unique(mut_data$patient_id) 

  #Vemos IDs comunes
common_ids<- intersect(count_matrix$patient_id, mutation_ids) #hay 965

  #Vemos diferencias 
setdiff(mutation_ids, count_matrix$patient_id) #pacientes con datos de mutaciones pero sin datos de expresión (3)
setdiff(count_matrix$patient_id, mutation_ids) #pacientes con datos de expresión sin datos de mutaciones (130)
setdiff(mutation_ids, clin_data$submitter_id) #comprobamos que no hay pacientes con datos de mutaciones sin datos clínicos

 #Filtramos  y guardamos
vst_counts_common <- vst_counts[count_matrix$patient_id %in% common_ids, ]
write.csv(vst_counts_common, file="count_matrix_norm_common.csv")

mut_common<-mut_data_filt[mut_data_filt$patient_id %in% common_ids, ]
saveRDS(mut_common, "mut_data_common.rds")

clin_common<-clin_data[clin_data$submitter_id %in% common_ids, ]
saveRDS(clin_common, "clin_data_common.rds")
