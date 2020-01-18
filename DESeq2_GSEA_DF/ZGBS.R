##---------------------------------------------
# Análise em Nível de Transcrito Kallisto
# Utilizando Octuplicatas
# alpha = 0.05 (FDR, padj < 0.05) + IWH
# Data: 18/01/2020
##---------------------------------------------
library(tximport)
library(DESeq2)
library(readr)
library(dplyr)
library(rhdf5)
library(IHW)