## ----setup, include=FALSE--------------------------------------------------------------------------------------------
options(htmltools.dir.version = FALSE)


## ----xaringan-themer, include=FALSE----------------------------------------------------------------------------------
library(xaringanthemer)
solarized_dark(
  code_font_family = "Fira Code",
  code_font_url    = "https://cdn.rawgit.com/tonsky/FiraCode/1.204/distr/fira_code.css"
)


## /* From https://github.com/yihui/xaringan/issues/147  */

## .scroll-output {

##   height: 80%;

##   overflow-y: scroll;

## }

## 
## /* https://stackoverflow.com/questions/50919104/horizontally-scrollable-output-on-xaringan-slides */

## pre {

##   max-width: 100%;

##   overflow-x: scroll;

## }

## 
## /* From https://github.com/yihui/xaringan/wiki/Font-Size */

## .tiny{

##   font-size: 40%

## }

## 
## /* From https://github.com/yihui/xaringan/wiki/Title-slide */

## .title-slide {

##   background-image: url(https://raw.githubusercontent.com/Bioconductor/OrchestratingSingleCellAnalysis/master/images/Workflow.png);

##   background-size: 33%;

##   background-position: 0% 100%

## }


## ----all_code, cache=TRUE--------------------------------------------------------------------------------------------
library('scRNAseq')
sce.zeisel <- ZeiselBrainData(ensembl = TRUE)

# Quality control
library('scater')
is.mito <- which(rowData(sce.zeisel)$featureType == "mito")
stats <- perCellQCMetrics(sce.zeisel, subsets = list(Mt = is.mito))
qc <-
    quickPerCellQC(stats,
        percent_subsets = c("altexps_ERCC_percent", "subsets_Mt_percent"))
sce.zeisel <- sce.zeisel[, !qc$discard]


## ----all_code2, cache=TRUE, dependson='all_code'---------------------------------------------------------------------
# Library size factors
#size factor: suma cuentas a lo largo de todos los genes para cada celula

lib.sf.zeisel <- librarySizeFactors(sce.zeisel)

# Examine distribution of size factors
summary(lib.sf.zeisel)
hist(log10(lib.sf.zeisel), xlab = "Log10[Size factor]", col = "grey80")
ls.zeisel <- colSums(counts(sce.zeisel))
plot(
    ls.zeisel,
    lib.sf.zeisel,
    log = "xy",
    xlab = "Library size",
    ylab = "Size factor"
)

##########################################################################################
                                   #QUESTIONS
##########################################################################################
#Are ls.zeisel and lib.sf.zeisel identical?
identical(lib.sf.zeisel, ls.zeisel)
#No

#Are they proporcional?
#vemos que son proporcionales solo con ver la imagen del plot, pero si b tiene el mismo valor
#para todo, es proporcional
b=(lib.sf.zeisel/ls.zeisel)
#todos nos dan lo mismo, entonces son proporcionales

## ----exercise_solution, cache=TRUE, dependson='all_code'-------------------------------------------------------------
## First compute the sums
zeisel_sums <- colSums(counts(sce.zeisel))
identical(zeisel_sums, ls.zeisel)

## Next, make them have unity mean
zeisel_size_factors <- zeisel_sums/mean(zeisel_sums)
identical(zeisel_size_factors, lib.sf.zeisel)


## ----all_code3, cache=TRUE, dependson='all_code2'--------------------------------------------------------------------
# Normalization by convolution

library('scran')
# Pre-clustering
set.seed(100)
clust.zeisel <- quickCluster(sce.zeisel)

table(clust.zeisel)
# Compute deconvolution size factors
deconv.sf.zeisel <-
    calculateSumFactors(sce.zeisel, clusters = clust.zeisel, min.mean = 0.1)

# Examine distribution of size factors
summary(deconv.sf.zeisel)
hist(log10(deconv.sf.zeisel), xlab = "Log10[Size factor]",
    col = "grey80")
plot(
    ls.zeisel,
    deconv.sf.zeisel,
    log = "xy",
    xlab = "Library size",
    ylab = "Size factor"
)

sce.zeisel <-logNormCounts(sce.zeisel)
assayNames(sce.zeisel)

##########################################################################################
                                         #QUESTIONS
##########################################################################################
#how many clusters?
#12
#how many cells per cluster?
summary(clust.zeisel)
#o usando table

#How many quick clusters will we get if we set the minimum size to 200? Use 100 as the seed.
set.seed(100)
clust.zeisel_2 <- quickCluster(sce.zeisel,min.size=200)
#Ahora hay menos clusters pero con mas celulas en cada uno

#Hay mas de una linea?
#Parece que hay dos, despues graficamos los puntos con diferentes colores por tipo celular (abajo)

## ----all_code4, cache=TRUE, dependson='all_code3'--------------------------------------------------------------------
# Library size factors vs. convolution size factors

# Colouring points using the supplied cell-types
plot(
    lib.sf.zeisel,
    deconv.sf.zeisel,
    xlab = "Library size factor",
    ylab = "Deconvolution size factor",
    log = 'xy',
    pch = 16,
    col = as.integer(factor(sce.zeisel$level1class))
)
abline(a = 0, b = 1, col = "red")


## ----'reproducibility', cache = TRUE, dependson=knitr::all_labels()--------------------------------------------------
options(width = 120)
sessioninfo::session_info()


#############################################################################################################
#Correccion de batch vs normalizacion 
#la normalizacion solo considera factores tecnicos, batch tambien considera factores biologicos correccion 
#de batch vs normalizacion 

#A veces el composition bias es muy pequeño y puede que solo sea necesario el library size normalization

#En expresion diferencial que es mas interensante? (genes son rows y columnas cells)
matrix(c(50,1000,10,1100),2,2)

#hacemos la transformacion a logaritmo para ver las diferencias relativas (10-50 vs 1000-11000) y asi podemos 
#usar diferencias euclidianas sumamos uno para no tener problemas con inf  y confirmamos que en este contexto 
#es mas interesante el gen1 