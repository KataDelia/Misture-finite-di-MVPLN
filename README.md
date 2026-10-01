# Misture finite di MVPLN per matrici a tre vie nell'RNA-seq

![R](https://img.shields.io/badge/r-%23276DC3.svg?style=for-the-badge&logo=r&logoColor=white)
![RStan](https://img.shields.io/badge/RStan-B81515.svg?style=for-the-badge&logo=R&logoColor=white)

> Progetto universitario di approfondimento e valutazione critica del lavoro:
> *Finite mixtures of matrix variate Poisson-log normal distributions for three-way count data* (Silva, Qin, Rothstein, McNicholas, Subedi, Bioinformatics, 39(5), 2023, doi:10.1093/bioinformatics/btad167).

## Descrizione

Il progetto analizza un metodo di classificazione basata su modello per dati di conteggio RNA-seq organizzati in **matrici a tre vie** (tempi × variabili × geni). Il metodo si basa su una mistura finita di distribuzioni MVPLN (*Matrix Variate Poisson-log Normal*), un'estensione matriciale della MPLN, che permette di:

- Modellare dati di conteggio sovradispersi, evitando il vincolo media-varianza della distribuzione di Poisson.
- Tenere conto delle correlazioni sia tra i tempi sia tra le variabili/condizioni.
- Ridurre il numero di parametri da stimare rispetto alla mistura MPLN, grazie alla struttura di covarianza separabile ($\Phi$ per i tempi, $\Omega$ per le variabili).

## Contenuti

- **Introduzione all'RNA-seq**: Da *reads* a *counts*, matrice di espressione.
- **Classificazione basata su modello**: Vantaggi (minore arbitrarietà, confronto tra modelli, invarianza di scala).
- **Modello MVPLN e sua mistura finita**: Assunzioni, funzione di densità, calcolo del numero di parametri.
- **Framework di stima**:
  - MCMC-EM (campionamento tramite Rstan)
  - Variational Gaussian Approximation (VGA), con massimizzazione dell'ELBO
  - Approccio ibrido (VGA + MCMC)
- **Identificabilità e Model Selection**: Criteri di informazione (BIC, ICL, AIC, AIC3) e Adjusted Rand Index.
- **Studio di simulazione**: Analisi di sei scenari e confronto con *HTSCluster* e *fuzzy k-means*.
- **Applicazioni in R (mixMVPLN)** su tre casi di studio:
  1. **Lievito a fissione (S. pombe)**: RNA non codificanti, ceppo standard vs mutante.
  2. **Ratti SHR vs WKY**: Nuclei autonomi del tronco encefalico e ipertensione.
  3. **Cavie µMT vs WT**: Tessuto polmonare e BPCO.

## Analisi Critica

L'aspetto distintivo di questo lavoro è la verifica pratica del pacchetto `mixMVPLN`, che ha portato alla luce diverse criticità implementative:

-  **Errori di esecuzione:** Interruzione dell'algoritmo a causa della presenza di matrici singolari.
-  **Limiti di scalabilità e convergenza:** Mancata convergenza con più di 2 cluster (riscontrata nel caso dei ratti) e forte instabilità elaborando oltre 200 geni.
-  **Incoerenza nei dati:** Utilizzo di dati microarray normalizzati (con metodo non specificato) trattati erroneamente come conteggi RNA-seq nell'articolo originale (caso delle cavie).

**Conclusione:** Sebbene l'intuizione metodologica sia innovativa e matematicamente ben formalizzata, l'implementazione software attuale non risulta sufficientemente robusta per applicazioni reali. L'ottimizzazione del codice e del metodo rappresenta un chiaro obiettivo per la ricerca futura.

## Strumenti

- **R**
- **mixMVPLN** (Pacchetto R)
- **Rstan**

## Autori

- Delia Anamaria Bogdan
- Maria Laura Fazio
- Samuele Confalone

## Riferimenti

- Silva A., Qin X., Rothstein S.J., McNicholas P.D., Subedi S. (2023). *Finite mixtures of matrix variate Poisson-log normal distributions for three-way count data*. Bioinformatics, 39(5), btad167.
- Silva A., Rothstein S.J., McNicholas P.D. et al. (2019). *A multivariate Poisson-log normal mixture model for clustering transcriptome sequencing data*. BMC Bioinformatics.
