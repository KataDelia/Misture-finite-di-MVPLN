# Progetto: Misture finite di MVPLN per matrici a tre vie nell'RNA-seq
# Autori: Delia Anamaria Bogdan, Maria Laura Fazio, Samuele Confalone

# 1. INSTALLAZIONE E CARICAMENTO PACCHETTI -------------------------------------
# Decommentare le righe seguenti per installare i pacchetti necessari
# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(version = "3.20")
# install.packages("devtools")
# BiocManager::install("edgeR")
# devtools::install_github("anjalisilva/mixMVPLN", build_vignettes = TRUE)
# install.packages("mclust") # Per il calcolo dell'Adjusted Rand Index

library(devtools)
library(edgeR)
library(mixMVPLN)
library(MASS)
library(clusterGeneration)
library(mclust) # aggiunto per valutare i cluster

# 2. CARICAMENTO DATI ----------------------------------------------------------
# INSERIRE QUI IL CODICE PER CARICARE I DATI NTS, RVLM, CVLM
# Esempio: 
# load("percorso/ai/tuoi/dati.RData")
# NTS <- read.csv("dati_NTS.csv", row.names=1) ...

# Verifica dimensioni preliminari
dim(NTS); dim(RVLM); dim(CVLM)

nomi_geni <- rownames(NTS) # 32800 geni
tempi <- c("Settimana 8", "Settimana 10", "Settimana 12", "Settimana 16", "Settimana 24")
nuclei <- c("CVLM", "RVLM", "NTS")

# 3. PREPARAZIONE DATASET (SHR e WKY) ------------------------------------------
# Seleziono una sola replica per ciascun tempo/nucleo
SHR <- data.frame(nomi_geni,
                  SHR8_N = NTS[, 1], SHR10_N = NTS[, 6], SHR12_N = NTS[, 11], SHR16_N = NTS[, 17], SHR24_N = NTS[, 23],
                  SHR8_R = RVLM[, 1], SHR10_R = RVLM[, 6], SHR12_R = RVLM[, 11], SHR16_R = RVLM[, 17], SHR24_R = RVLM[, 23],
                  SHR8_C = CVLM[, 1], SHR10_C = CVLM[, 6], SHR12_C = CVLM[, 11], SHR16_C = CVLM[, 17], SHR24_C = CVLM[, 23])

WKY <- data.frame(nomi_geni,
                  WKY8_N = NTS[, 29], WKY10_N = NTS[, 34], WKY12_N = NTS[, 38], WKY16_N = NTS[, 41], WKY24_N = NTS[, 44],
                  WKY8_R = RVLM[, 29], WKY10_R = RVLM[, 34], WKY12_R = RVLM[, 38], WKY16_R = RVLM[, 41], WKY24_R = RVLM[, 44],
                  WKY8_C = CVLM[, 29], WKY10_C = CVLM[, 34], WKY12_C = CVLM[, 38], WKY16_C = CVLM[, 41], WKY24_C = CVLM[, 44])

# 4. CREAZIONE MATRICI A TRE VIE -----------------------------------------------
# Funzione helper per creare le matrici 5x3 per ogni gene
crea_matrice <- function(data_row, prefix) {
  mat <- matrix(
    as.numeric(data_row[2:16]), # Evita di prendere il nome del gene
    nrow = 5, ncol = 3, byrow = FALSE,
    dimnames = list(tempi, c("NTS", "RVLM", "CVLM"))
  )
  return(mat)
}

lista_matrici_SHR <- lapply(1:nrow(SHR), function(i) crea_matrice(SHR[i, ]))
names(lista_matrici_SHR) <- paste0(SHR$nomi_geni, "_SHR")

lista_matrici_WKY <- lapply(1:nrow(WKY), function(i) crea_matrice(WKY[i, ]))
names(lista_matrici_WKY) <- paste0(WKY$nomi_geni, "_WKY")

# Unione delle due liste
lista_matrici <- c(lista_matrici_WKY, lista_matrici_SHR)

# 5. PULIZIA E SOTTOCAMPIONAMENTO ----------------------------------------------
# Manteniamo solo interi positivi 
# NOTA: mixMVPLN richiede valori > 0. Questo step filtra i geni con conteggi pari a zero.
lista_matrici_interi_positivi <- lista_matrici[sapply(lista_matrici, function(mat) {
  all(mat == floor(mat)) && all(mat > 0) 
})]

num_matrici <- length(lista_matrici_interi_positivi)
cat("Geni rimasti dopo il filtro > 0:", num_matrici, "\n")

# Estrazione casuale di un sottoinsieme per limiti computazionali (es. 200 geni)
set.seed(123) 
matrici_casuali <- lista_matrici_interi_positivi[sample(1:num_matrici, 200)]

# Proporzioni di campionamento
nomi_casuali <- names(matrici_casuali)
gruppi_veri <- substr(nomi_casuali, nchar(nomi_casuali) - 2, nchar(nomi_casuali))

cat("Proporzione WKY estratti:", mean(gruppi_veri == "WKY"), "\n")
cat("Proporzione SHR estratti:", mean(gruppi_veri == "SHR"), "\n")

# 6. ESECUZIONE MODELLI CLUSTERING ---------------------------------------------

# --- 1. Metodo VGA (Variational Gaussian Approximation) ---
VGA <- mvplnVGAclus(dataset = matrici_casuali,
                    gmin = 1, gmax = 2,
                    normalize = "Yes")

cluster_VGA <- VGA$allResults[[2]]$clusterlabels

# Valutazione tramite Adjusted Rand Index per evitare il problema del label switching
ari_VGA <- adjustedRandIndex(cluster_VGA, gruppi_veri)
cat("VGA - Adjusted Rand Index:", round(ari_VGA, 4), "\n")


# --- 2. Metodo MCMC ---
# ATTENZIONE: nIterations impostato a 40 solo a scopo esplorativo/didattico.
# Per inferenza reale aumentare significativamente (es. 2000+).
MCMC <- mvplnMCMCclus(dataset = matrici_casuali,
                      gmin = 1, gmax = 2,
                      nIterations = 40,
                      normalize = "Yes")


# --- 3. Metodo Ibrido ---
HYB <- mvplnHybriDclus(dataset = matrici_casuali,
                       gmin = 1, gmax = 2,
                       normalize = "Yes")