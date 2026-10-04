#############################################################################################
################ DESCARGA DE DATOS DE LA BASE DE DATOS GENOMIC DATA COMMONS ###################
#############################################################################################

# info obtenida de https://www.bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/query.html#Harmonized_data_options
setwd("C:/Users/Usuario/Desktop/Bioinformatica/TFM/")

## 1.- Carga librería TCGAbiolinks ----

library(TCGAbiolinks)

#Funciones clave:
 #- GDCquery(): definir qué quieres
 #- GDCdownload(): descargarlo
 #- GDCprepare(): convertirlo en objeto usable


## 2.- Listamos proyectos disponibles en la BBDD Genomic Data Commons ----

projects <- getGDCprojects()
  #data frame con los 91 proyectos

head(projects$project_id)


## 3.- Búsqueda ----

#Primero vamos a ver que tipo de datos tenemos dentro del proyecto TCGA-BRCA 

getProjectSummary("TCGA-BRCA")

#Finalmente hacemos la búsqueda de los datasets de interés con la función GDCquery 

query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling", #buscamos la expresión génica
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts", #conteo de las lecturas tras el alineamiento con STAR
  sample.type = "Primary Tumor" #tipo de muestra
)

## 4.- Descarga ----

GDCdownload(query, method = "api", files.per.chunk = 10)


## 5.- Preparar los datos ----

exp_data<-GDCprepare(query)
#Nos devuelve un objeto tipo SummarizedExperiment, que contiene 3 matrices principales:
 # - información de las muestras (metadatos): se accede con la función colData()
 # - matriz del ensayo: se accede con la función assay()
 # - matriz de información de las variables (genes): se accede con la función rowRanges()

#Para ejecutar estas funciones necesitamos cargar la librería SumarizedExperiment:
library(SummarizedExperiment)

### 5.1.- Sacar matriz de expresión ----

count_matrix <- assay(exp_data)  
View(count_matrix)
#se observa que la matriz tiene los genes como observaciones (filas) y los pacientes como columnas (n=1111).

### 5.2.- Obtener metadatos ----

metadata<-as.data.frame(colData(exp_data))
 #contiene información sobre las pacientes y la toma de muestras (tejido y preservación, edad al diagnóstico, estadio, etc)

### 5.2.- Obtener información genes ----

genes<-as.data.frame(rowRanges(exp_data))
 #de aquí podemos obtener el nombre del gen a partir del ID del transcrito 


## 6-. REPETIMOS PROCESO para los datos de MUTACIONES SOMÁTICAS e INFORMACIÓN CLÍNICA ----

### 6.1.- Búsqueda ----

getProjectSummary("TCGA-BRCA")

query_mut<-GDCquery(
  project = "TCGA-BRCA",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation", #buscamos las mutaciones somáticas
)

query_clinical<-GDCquery(
  project = "TCGA-BRCA",
  data.category = "Clinical", #buscamos la información clínica
  data.type = "Clinical Supplement",
  data.format = "BCR XML"
)

#OTRA FORMA PARA OBTENER LOS DATOS CLÍNICOS: con la función GDCquery_clinic 
#Te devuelve directamente un dataframe listo para usar
clinical <- GDCquery_clinic("TCGA-BRCA", type = "clinical")

### 6.2.- Descargar ----

GDCdownload(query_mut, method = "api", files.per.chunk = 10)
 
GDCdownload(query_clinical, method = "api", files.per.chunk = 10)

### 6.3. - Preparar datos ----

mut_data<-GDCprepare(query_mut)
#devuelve un dataframe con todas las mutaciones encontradas (+89000).en cada paciente hay mas de 1 gen mutado

clin_data<-GDCprepare_clinic(query_clinical, clinical.info = "patient") 
#devuelve un dataframe con información clínica para cada paciente

####!!! Hay información faltante en este dataframe en comparación al obtenido en la variable clinical


## 7-. GUARDAMOS datos en objeto RDS ----
# - Para poder recuperarlos fácilmente luego sin tener que volver a hacer todo lo anterior

saveRDS(exp_data, "expression.rds")
saveRDS(clinical, "clinical.rds")
saveRDS(mut_data, "mutations.rds")




