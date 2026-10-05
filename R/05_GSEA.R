## =========================================================
## GSEA preranked de resultados DESeq2
## =========================================================

# Tutoriales: https://www.zubairkhalid.com/knowledge/bioinformatics/how-to-run-gene-set-enrichment-analysis-gsea-with-fgsea-in-r-a-complete-workflow-for-rna-seq
#             https://biostatsquid.com/gene-set-enrichment-analysis/ 
#            https://biostatsquid.com/fgsea-tutorial-gsea/


rm(list = ls())

setwd("C:/Users/Usuario/Desktop/Bioinformatica/TFM")

library(dplyr)
library(msigdbr)
library(fgsea)
library(ggplot2)


## 1. Cargar resultados de DESeq2 ----

res_deseq <- read.csv("DESeq2_train/results_deseq.csv")


## 2. Crear el ranking de genes ----

#a partir de "stat"

sum(is.na(res_deseq))

ranking <- res_deseq %>%
  filter(
    !is.na(stat), #filtramos valores vacíos o na
    !is.na(X),
    X != ""
  ) %>%
  group_by(X) %>% #agrupamos por gen
  slice_max(
    order_by = abs(stat), #ordenamos por el valor de stat
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()

geneList <- ranking$stat #lista de genes con su stat
names(geneList) <- ranking$X

geneList <- sort(geneList, decreasing = TRUE)


## Comprobar el ranking
head(geneList)
tail(geneList)


## 3. Obtener conjuntos génicos Hallmark de MSigDB ---- 

hallmark <- msigdbr(
  db_species = "HS",
  species = "Homo sapiens",
  collection = "H"
)


## 4. Convertir MSigDB al formato que necesita fgsea ----

pathways <- split(  #agrupamos los genes por vía (gs_name)
  hallmark$gene_symbol,
  hallmark$gs_name
)


## - Comprobar número de vías
length(pathways)


## - Comprobar cuántos genes de nuestro ranking están presentes en MSigDB -----

sum(names(geneList) %in% hallmark$gene_symbol)

length(intersect(
  names(geneList),
  hallmark$gene_symbol
))


## 5. Ejecutar GSEA preranked ----

set.seed(123)

gsea_res <- fgsea(
  pathways = pathways, #lista de pathways con los genes
  stats = geneList, #ranking de nuestros genes
  minSize = 15, #minimo de genes por set
  maxSize = 500 #maximo de genes por set
)


## 6. Ordenar resultados por FDR ----

gsea_res <- gsea_res %>%
  arrange(padj)


## Ver resultados
gsea_res[, c(
  "pathway",
  "NES",
  "pval",
  "padj",
  "size"
)]


## 7. Seleccionar vías significativas ----

gsea_sig <- gsea_res %>%
  filter(padj < 0.05)

gsea_sig #29 vías significativas


## Guardar resultados gsea
write.csv(
  gsea_res[,-"leadingEdge"], #no guardamos la ultima columna que es una lista
  "GSEA_results.csv",
  row.names = FALSE
)


## 8. Representar las vías significativas ----

top_gsea <- gsea_sig %>%
  arrange(desc(abs(NES))) %>% #ordenar pot valor absoluto de NES
  slice_head(n = 15)

#Limpiamos la columna pathway para que se vea mejor en la gráfica
top_gsea<-top_gsea %>%
  mutate(
    pathway_clean = pathway %>%
      gsub("^HALLMARK_", "", .) %>%
      gsub("_", " ", .)
  )


ggplot(top_gsea, 
       aes(x = NES, 
           y = reorder(pathway_clean, NES), #ordenamos paths según NES
           size = -log10(padj), color=NES)) +
  geom_point() +
  scale_color_gradient2(
    low = "steelblue4",
    mid = "white",
    high = "firebrick3",
    midpoint = 0 ) +
  theme_minimal() +
  labs( x = "Normalized Enrichment Score (NES)",
        y = NULL,
        size = "-log10(FDR)") 

ggsave("plots/plot_GSEA.png", dpi=300)


# Vamos a explorar los genes de algún pathway

#LeadingEdge hace referencia a la columna con los genes que mas contribuyen al
#enriquecimiento de la vía

gsea_res[1]$leadingEdge


# 9. Típico gráfico de barras de enriquecimiento ----

#Para cada pathway individual

plotEnrichment(
     pathways[["HALLMARK_E2F_TARGETS"]],
     geneList
 )
 #el"codigo de barras" representa la posición en el ranking de cada gen. como esta vía
# está sobrerrepresentada se observan más barras en la derecha

plotEnrichment(
  pathways[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]],
  geneList
)
