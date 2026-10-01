# Simulazione: Dataset GSE213192 (Matrici 2x2)
# Obiettivo: Testare mixMVPLN su matrici a tre vie ridotte (Tempi x Pazienti)
# per verificare la separazione tra tipi cellulari (CD4 vs CD8).

library(mixMVPLN)
library(mclust) # Per Adjusted Rand Index

# 1. CARICAMENTO DATI ----------------------------------------------------------
# Sostituire con il percorso corretto del file scaricato da GEO
# GSE213192_raw_counts <- read.csv("percorso/a/GSE213192_raw_counts.csv", row.names=1)

head(GSE213192_raw_counts)

# Selezione delle colonne di interesse (da adattare se gli indici cambiano)
subset_data_bambino1 <- GSE213192_raw_counts[, 4:7]
subset_data_bambino2 <- GSE213192_raw_counts[, 10:13]

data <- cbind(subset_data_bambino1, subset_data_bambino2)

# Filtro esplorativo: manteniamo solo i geni molto espressi (somma >= 100)
data <- data[rowSums(data >= 100) > 0, ]

# 2. CREAZIONE DELLE MATRICI ---------------------------------------------------
# Estraiamo i primi 100 geni per limitare i tempi di calcolo
num_geni <- 100
nomi_geni <- rownames(data)[1:num_geni]

matrices_list1 <- list()
matrices_list2 <- list()

for (i in 1:num_geni) {
  # Matrice per CD4 (Tempi: AC, FU x Pazienti: bambino1, bambino2)
  matrix_cd4 <- matrix(
    c(data$X4580.AC.CD4[i], data$X9630.AC.CD4[i], 
      data$X4580.FU.CD4[i], data$X9630.FU.CD4[i]),
    nrow = 2, byrow = TRUE,
    dimnames = list(c("AC", "FU"), c("bambino1", "bambino2"))
  )
  matrices_list1[[i]] <- matrix_cd4
  
  # Matrice per CD8
  matrix_cd8 <- matrix(
    c(data$X4580.AC.CD8[i], data$X9630.AC.CD8[i], 
      data$X4580.FU.CD8[i], data$X9630.FU.CD8[i]),
    nrow = 2, byrow = TRUE,
    dimnames = list(c("AC", "FU"), c("bambino1", "bambino2"))
  )
  matrices_list2[[i]] <- matrix_cd8
}

# Aggiungiamo i nomi alle liste per mantenere traccia del gruppo reale (Ground Truth)
names(matrices_list1) <- paste0(nomi_geni, "_CD4")
names(matrices_list2) <- paste0(nomi_geni, "_CD8")

# Uniamo in un unico dataset di 200 matrici
dataset_completo <- c(matrices_list1, matrices_list2)

# Estraiamo i gruppi reali per il calcolo dell'accuratezza
nomi_completi <- names(dataset_completo)
gruppi_veri <- substr(nomi_completi, nchar(nomi_completi) - 2, nchar(nomi_completi))

# 3. ESECUZIONE DEL MODELLO (VGA) ----------------------------------------------
# Stimiamo 1 o 2 cluster per vedere se il modello riconosce CD4 e CD8
modello_VGA <- mvplnVGAclus(dataset = dataset_completo,
                            gmin = 1,
                            gmax = 2,
                            normalize = "Yes")

# 4. VALUTAZIONE DEI RISULTATI -------------------------------------------------
# Estraiamo le etichette predette per l'esecuzione con 2 cluster
cluster_predetti <- modello_VGA$allResults[[2]]$clusterlabels

# Calcoliamo l'Adjusted Rand Index
ari_score <- adjustedRandIndex(cluster_predetti, gruppi_veri)

cat("\n=== RISULTATI CLUSTERING ===\n")
cat("Adjusted Rand Index (CD4 vs CD8):", round(ari_score, 4), "\n")
table(Veri = gruppi_veri, Predetti = cluster_predetti)