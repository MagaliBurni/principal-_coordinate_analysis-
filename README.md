###############################################
# Script: Análisis de Coordenadas Principales (PCoA)
# Autora: Magali Burni
# Fecha: diciembre 2024
###############################################

# Análisis de Coordenadas Principales (PCoA) para caracteres Funcionales

Este repositorio contiene un script en R para realizar un **Análisis de Coordenadas Principales (PCoA)** con datos de caracteres funcionales de especies leñosas.  
El análisis incluye chequeo de colinealidad, cálculo de distancias de Gower, ordenación por PCoA, modelos GLMM, diagnóstico estadístico y visualización de resultados mediante gráficos de PCoA y boxplots.

---

## Contenido del repositorio

- `pcoa_final.R`: script principal con todo el análisis.  
- Este `README.md` con la descripción del script.  

**Nota:** La base de datos (`pca_traits.csv`) no está incluida en este repositorio por motivos de confidencialidad. Para ejecutar el script, es necesario contar con un archivo de datos con estructura similar (columnas de especies, caracteres funcionales y categoría de status (invasora, no invasora).

---

## Requisitos

El script fue desarrollado en R y requiere los siguientes paquetes:

```r
library(performance)
library(vegan)
library(patchwork)
library(DHARMa)
library(glmmTMB)
library(car)
library(ggplot2)
