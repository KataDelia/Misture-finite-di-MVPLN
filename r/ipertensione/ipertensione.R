rm(list=ls())

library("devtools")
library("edgeR")
library("mixMVPLN")
library("MASS")
library("clusterGeneration")

# verifico i tre dataset per i tre nuclei
NTS <- read.csv("NTS.txt", sep = "\t")
CVLM <- read.csv("CVLM.txt", sep = "\t")
RVLM <- read.csv("RVLM.txt", sep = "\t")

head(NTS)
dim(NTS)
names(NTS)

nomi_geni <- rownames(NTS) # 32800 geni
tempi <- c("Settimana 8", "Settimana 10", "Settimana 12", "Settimana 16", "Settimana 24")
nuclei <- c("CVLM", "RVLM", "NTS")

# due dataset per ciascuna tipologia di topo.
# Seleziono una sola replica per ciascun tempo/nucleo
# possiamo usare il valore mediano per ogni replicato
SHR <- data.frame(nomi_geni,
                  SHR8_N = NTS[, 1],
                  SHR10_N = NTS[, 6],
                  SHR12_N = NTS[, 11],
                  SHR16_N = NTS[, 17],
                  SHR24_N = NTS[, 23],
                  SHR8_R = RVLM[, 1],
                  SHR10_R = RVLM[, 6],
                  SHR12_R = RVLM[, 11],
                  SHR16_R = RVLM[, 17],
                  SHR24_R = RVLM[, 23],
                  SHR8_C = CVLM[, 1],
                  SHR10_C = CVLM[, 6],
                  SHR12_C = CVLM[, 11],
                  SHR16_C = CVLM[, 17],
                  SHR24_C = CVLM[, 23])

WKY <- data.frame(nomi_geni,
                  WKY8_N = NTS[, 29],
                  WKY10_N = NTS[, 34],
                  WKY12_N = NTS[, 38],
                  WKY16_N = NTS[, 41],
                  WKY24_N = NTS[, 44],
                  WKY8_R = RVLM[, 29],
                  WKY10_R = RVLM[, 34],
                  WKY12_R = RVLM[, 38],
                  WKY16_R = RVLM[, 41],
                  WKY24_R = RVLM[, 44],
                  WKY8_C = CVLM[, 29],
                  WKY10_C = CVLM[, 34],
                  WKY12_C = CVLM[, 38],
                  WKY16_C = CVLM[, 41],
                  WKY24_C = CVLM[, 44])

# funzioni che mi creano una matrice in riferimento a ciascun gene
# Otteniamo due liste di matrici 5 x 2
# aggiungo un'etichetta per i geni che riporta alla tipologia di topo

lista_matrici_N <- lapply(1:nrow(SHR), function(i) {
  matrice <- matrix(
    c(
      SHR[i, "SHR8_N"], SHR[i, "SHR10_N"], SHR[i, "SHR12_N"], SHR[i, "SHR16_N"], SHR[i, "SHR24_N"],
      WKY[i, "WKY8_N"], WKY[i, "WKY10_N"], WKY[i, "WKY12_N"], WKY[i, "WKY16_N"], WKY[i, "WKY24_N"]
    ),
    nrow = 5, ncol = 2, byrow = FALSE,
    dimnames = list(tempi, c("RHS", "WKY"))
  )
  matrice
})

names(lista_matrici_N) <- lista_matrici_N$nomi_geni

# pulisco le matrici estraendo quelle intere positive da tutta la lista
lista_matrici_interi_positivi <- lista_matrici_N[sapply(lista_matrici_N, function(matrice) {
  all(matrice == floor(matrice)) && all(matrice > 100) })]

length(lista_matrici_interi_positivi)
num_matrici <- length(lista_matrici_interi_positivi)

# dei geni rimasti ne estraggo casualmente un tot. (in questo caso 200) 
set.seed(1234) 
matrici_casuali <- lista_matrici_interi_positivi[sample(1:num_matrici, 200)]
head(matrici_casuali)

# provo con 200 geni e due cluster, l'algoritmo gira
# provo con 1000 geni e due cluster, l'algoritmo non gira 

VGA = mvplnVGAclus(dataset = matrici_casuali,
                   gmin = 1,
                   gmax = 2,
                   normalize = "Yes")
table(VGA$AICAll$AICmodelselectedLabels)

# provo con 200 geni e tre cluster, l'algoritmo non gira
VGA = mvplnVGAclus(dataset = matrici_casuali,
                   gmin = 3,
                   gmax = 3,
                   normalize = "Yes")


# provo con 1000 geni e due cluster, l'algoritmo non gira
matrici_casuali <- lista_matrici_interi_positivi[sample(1:num_matrici, 1000)]
head(matrici_casuali)

VGA = mvplnVGAclus(dataset = matrici_casuali,
                   gmin = 1,
                   gmax = 2,
                   normalize = "Yes")
