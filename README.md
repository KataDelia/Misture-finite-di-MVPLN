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
- **Framework di stima**: MCMC-EM, Variational Gaussian Approximation (VGA) e Approccio ibrido (VGA + MCMC).
- **Identificabilità e Model Selection**: Criteri di informazione (BIC, ICL, AIC, AIC3) e Adjusted Rand Index (ARI).
- **Applicazione Pratica (Script R): Caso di studio Ratti SHR vs WKY**
  - Strutturazione dei dati trascrittomici in matrici 5x3 (5 tempi: settimane 8, 10, 12, 16, 24 × 3 nuclei del tronco encefalico: CVLM, RVLM, NTS).
  - Preprocessing dei dati e sottocampionamento strategico per far fronte ai limiti computazionali.
  - Valutazione delle performance di clustering tramite **Adjusted Rand Index** (pacchetto `mclust`) per neutralizzare il problema del *label switching* tipico degli algoritmi non supervisionati.

## Analisi Critica

L'aspetto distintivo di questo lavoro è la verifica pratica del pacchetto `mixMVPLN` su dataset reali, che ha portato alla luce diverse criticità implementative rilevanti:

-  **Gestione degli zeri (Filtro forzato):** L'algoritmo va in errore in presenza di valori nulli nelle matrici (a causa del logaritmo insito nella componente log-normale). Ciò costringe a eliminare a priori i geni con conteggi pari a zero, una forzatura tecnica problematica dato che i dati RNA-seq sono fisiologicamente ricchi di zeri (geni non espressi).
-  **Costo computazionale estremo:** L'approccio MCMC richiede tempi di elaborazione proibitivi. Per testare l'algoritmo localmente in tempi ragionevoli è necessario campionare un numero ristrettissimo di geni e impostare iterazioni irrisorie, rendendolo poco praticabile per pipeline trascrittomiche complete.
-  **Errori di esecuzione:** Interruzione frequente dell'algoritmo a causa della formazione di matrici singolari durante la stima.
-  **Limiti di scalabilità e convergenza:** Mancata convergenza con più di 2 cluster (riscontrata direttamente nel caso dei ratti) e forte instabilità elaborando oltre 200 geni.
-  **Incoerenza nei dati (Paper originale):** Utilizzo di dati microarray normalizzati (con metodo non specificato) trattati erroneamente come conteggi interi positivi RNA-seq nell'articolo di riferimento (caso delle cavie).

**Conclusione:** Sebbene l'intuizione metodologica sia innovativa e matematicamente ben formalizzata, l'implementazione software attuale non risulta sufficientemente robusta né scalabile per applicazioni bioinformatiche reali. 

## Riproducibilità dello Script

Il repository include lo script R ottimizzato (`script_analisi.R`) per replicare il clustering sui ratti SHR e WKY. Le fasi principali includono:
1. Combinazione dei dataset dei tre nuclei (NTS, RVLM, CVLM) per creare matrici a tre vie dinamiche.
2. Pulizia dei dati (filtro sui positivi per i limiti del pacchetto) e campionamento casuale stratificato.
3. Esecuzione dei tre framework di stima (`mvplnVGAclus`, `mvplnMCMCclus`, `mvplnHybriDclus`).
4. Valutazione finale dell'accuratezza della partizione.

## Strumenti e Pacchetti

- **R**
- **mixMVPLN**, **Rstan** (Modellazione)
- **edgeR**, **MASS**, **clusterGeneration** (Manipolazione dati)
- **mclust** (Valutazione clustering)

## Autori

- Delia Anamaria Bogdan
- Maria Laura Fazio
- Samuele Confalone

## Riferimenti

- Silva A., Qin X., Rothstein S.J., McNicholas P.D., Subedi S. (2023). *Finite mixtures of matrix variate Poisson-log normal distributions for three-way count data*. Bioinformatics, 39(5), btad167.
- Silva A., Rothstein S.J., McNicholas P.D. et al. (2019). *A multivariate Poisson-log normal mixture model for clustering transcriptome sequencing data*. BMC Bioinformatics.
