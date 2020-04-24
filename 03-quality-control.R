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
## Data
library('scRNAseq')
sce.416b <- LunSpikeInData(which = "416b")
sce.416b$block <- factor(sce.416b$block)

# Download the relevant Ensembl annotation database
# using AnnotationHub resources
library('AnnotationHub')
ah <- AnnotationHub()
query(ah, c("Mus musculus", "Ensembl", "v97"))
# Annotate each gene with its chromosome location
ens.mm.v97 <- ah[["AH73905"]]
location <- mapIds(
    ens.mm.v97,
    keys = rownames(sce.416b),
    keytype = "GENEID",
    column = "SEQNAME"
)
# Identify the mitochondrial genes
is.mito <- which(location == "MT")

#Subset by mitochondrial genes
library('scater')
sce.416b <- addPerCellQC(sce.416b,
    subsets = list(Mito = is.mito))

##########################################################################################
                                      #QUESTIONS
##########################################################################################

#What changed in our sce object after addPerCellQC?
#El numero de colData names era 9 y ahora 25, esto cambia en Metadata y ahora inclye datos
#de la calidad de cada celula

## ----qc_metrics, cache=TRUE, dependson='all_code'--------------------------------------------------------------------
plotColData(sce.416b, x = "block", y = "detected")

plotColData(sce.416b, x = "block", y = "detected") +
    scale_y_log10()

plotColData(sce.416b,
    x = "block",
    y = "detected",
    other_fields = "phenotype") +
    scale_y_log10() +
    facet_wrap( ~ phenotype)


## ----all_code_part2, cache = TRUE, dependson='all_code'--------------------------------------------------------------
# Example thresholds
#The cell quality is evaluated under these conditions-> fixed
qc.lib <- sce.416b$sum < 100000
qc.nexprs <- sce.416b$detected < 5000
qc.spike <- sce.416b$altexps_ERCC_percent > 10
qc.mito <- sce.416b$subsets_Mito_percent > 10
discard <- qc.lib | qc.nexprs | qc.spike | qc.mito

# Summarize the number of cells removed for each reason
DataFrame(
    LibSize = sum(qc.lib),
    NExprs = sum(qc.nexprs),
    SpikeProp = sum(qc.spike),
    MitoProp = sum(qc.mito),
    Total = sum(discard)
)

#The cell quality is evaluated under these conditions-> adaptive, data outliers
#Assumes mst of them to be high-quality
qc.lib2 <- isOutlier(sce.416b$sum, log = TRUE, type = "lower")
qc.nexprs2 <- isOutlier(sce.416b$detected, log = TRUE,
    type = "lower")
qc.spike2 <- isOutlier(sce.416b$altexps_ERCC_percent,
    type = "higher")
qc.mito2 <- isOutlier(sce.416b$subsets_Mito_percent,
    type = "higher")
discard2 <- qc.lib2 | qc.nexprs2 | qc.spike2 | qc.mito2

# Extract the thresholds
attr(qc.lib2, "thresholds")
attr(qc.nexprs2, "thresholds")
# Summarize the number of cells removed for each reason.
DataFrame(
    LibSize = sum(qc.lib2),
    NExprs = sum(qc.nexprs2),
    SpikeProp = sum(qc.spike2),
    MitoProp = sum(qc.mito2),
    Total = sum(discard2)
)
##########################################################################################
                                    #QUESTIONS
##########################################################################################
#Was qc.lib necessary for creating discord?
y=seq(1,length(qc.lib))
y=y[(qc.lib)&(!qc.spike)&(!qc.mito)&(!qc.nexprs)]
#This tells us that 2 cells were discarder under this condition, meaning it was necessary 
#Another way:
which(qc.lib&(!qc.spike&!qc.mito&!qc.nexprs))

#Which filter discarded more cells? discard or discard2?
length(which(discard))
length(which(discard2))
#There are more TRUE in discard than in discard2, discard discarded more cells

#By considering the sample batch, did we discard more cells using automatic threshold detection?
length(which(discard3))

## More checks
plotColData(sce.416b,
    x = "block",
    y = "detected",
    other_fields = "phenotype") +
    scale_y_log10() +
    facet_wrap( ~ phenotype)

batch <- paste0(sce.416b$phenotype, "-", sce.416b$block)
qc.lib3 <- isOutlier(sce.416b$sum,
    log = TRUE,
    type = "lower",
    batch = batch)
qc.nexprs3 <- isOutlier(sce.416b$detected,
    log = TRUE,
    type = "lower",
    batch = batch)
qc.spike3 <- isOutlier(sce.416b$altexps_ERCC_percent,
    type = "higher",
    batch = batch)
qc.mito3 <- isOutlier(sce.416b$subsets_Mito_percent,
    type = "higher",
    batch = batch)
discard3 <- qc.lib3 | qc.nexprs3 | qc.spike3 | qc.mito3

# Extract the thresholds
attr(qc.lib3, "thresholds")
attr(qc.nexprs3, "thresholds")

# Summarize the number of cells removed for each reason
DataFrame(
    LibSize = sum(qc.lib3),
    NExprs = sum(qc.nexprs3),
    SpikeProp = sum(qc.spike3),
    MitoProp = sum(qc.mito3),
    Total = sum(discard3)
)


## ----use_case, cache=TRUE, dependson= c('all_code', 'all_code_part2')------------------------------------------------
sce.grun <- GrunPancreasData()
sce.grun <- addPerCellQC(sce.grun)

#el eje y es el porcentaje de lecturas de mapearon contra secuencias de ERCC, en este caso queremos que estén cerca de cero 
#n el droplet medimos el mensajero de la celula y los ERCC (mensajeros preparados, son 92 secuencias diseñadas para que funcionen bien con PCR)
#Paso de control para ver cuantas lecturas se usaron para ERCC y para los mensajeros, los que tienen la cola tan cargada no nos sirven
#Haciendo un histograma del donador 10 vemos que no es simetrico y la mayoria de las muestras salieron mal
#vs. histograma de D2, sigue teniendo cola pero es menos cargada y la podemos eliminar
plotColData(sce.grun, x = "donor", y = "altexps_ERCC_percent")

discard.ercc <- isOutlier(sce.grun$altexps_ERCC_percent,
    type = "higher",
    batch = sce.grun$donor)

#fijamos haciendo el discard respecto a los que salieron bien (escogiendolas manualmente),
#podemos descartar las células raras de D10 y D3
discard.ercc2 <- isOutlier(
    sce.grun$altexps_ERCC_percent,
    type = "higher",
    batch = sce.grun$donor,
    subset = sce.grun$donor %in% c("D17", "D2", "D7")
)

plotColData(
    sce.grun,
    x = "donor",
    y = "altexps_ERCC_percent",
    colour_by = data.frame(discard = discard.ercc)
)
plotColData(
    sce.grun,
    x = "donor",
    y = "altexps_ERCC_percent",
    colour_by = data.frame(discard = discard.ercc2)
)

# Add info about which cells are outliers
sce.416b$discard <- discard2

# Look at this plot for each QC metric
plotColData(
    sce.416b,
    x = "block",
    y = "sum",
    colour_by = "discard",
    other_fields = "phenotype"
) +
    facet_wrap( ~ phenotype) +
    scale_y_log10()

# Another useful diagnostic plot
plotColData(
    sce.416b,
    x = "sum",
    y = "subsets_Mito_percent",
    colour_by = "discard",
    other_fields = c("block", "phenotype")
) +
    facet_grid(block ~ phenotype)


## ----use_case_pbmc, cache=TRUE, dependson='all_code'-----------------------------------------------------------------
library('BiocFileCache')
bfc <- BiocFileCache()
raw.path <-
    bfcrpath(
        bfc,
        file.path(
            "http://cf.10xgenomics.com/samples",
            "cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices.tar.gz"
        )
    )
untar(raw.path, exdir = file.path(tempdir(), "pbmc4k"))

library('DropletUtils')
library('Matrix')
fname <- file.path(tempdir(), "pbmc4k/raw_gene_bc_matrices/GRCh38")
sce.pbmc <- read10xCounts(fname, col.names = TRUE)

bcrank <- barcodeRanks(counts(sce.pbmc))

# Only showing unique points for plotting speed.
uniq <- !duplicated(bcrank$rank)
plot(
    bcrank$rank[uniq],
    bcrank$total[uniq],
    log = "xy",
    xlab = "Rank",
    ylab = "Total UMI count",
    cex.lab = 1.2
)
abline(h = metadata(bcrank)$inflection,
    col = "darkgreen",
    lty = 2)
abline(h = metadata(bcrank)$knee,
    col = "dodgerblue",
    lty = 2)
legend(
    "bottomleft",
    legend = c("Inflection", "Knee"),
    col = c("darkgreen", "dodgerblue"),
    lty = 2,
    cex = 1.2
)

#lo de abajo de la linea verde los consideramos droplets vacios
#porque tienen muy pocas leturas para considerar que realmente hay celulas 

#empty drops: usa simluaciones de MonteCarlo
#pefil de expresion basado en los driplets ambientales y va a preguntar si las celulas de en medio son diferentes o no
#para determinar la "knee" de la curva buscan el punto de inflecion con las derivadas
set.seed(100)
e.out <- emptyDrops(counts(sce.pbmc))

# See ?emptyDrops for an explanation of why there are NA # values.
summary(e.out$FDR <= 0.001)

set.seed(100)
limit <- 100
all.out <-
    emptyDrops(counts(sce.pbmc), lower = limit, test.ambient = TRUE)
# Ideally, this histogram should look close to uniform.
# Large peaks near zero indicate that barcodes with total
# counts below 'lower' are not ambient in origin.
hist(all.out$PValue[all.out$Total <= limit &
        all.out$Total > 0],
    xlab = "P-value",
    main = "",
    col = "grey80")

#se usa 100 y más o menos se calcula a ojo, oara determinar cual es tu conjuto de empty droplets
#Histograma para ver si el valor de lower que escogimos es correcto o no (tomamos las que no tenian valor cero)
#Hipotesis nula: las de en medio son iguales a las lower
#los p-values de los droplets tienen que tener una distribucion uniforme si la hipotesis nula se cumple

sce.pbmc <- sce.pbmc[, which(e.out$FDR <= 0.001)]

is.mito <- grep("^MT-", rowData(sce.pbmc)$Symbol)
sce.pmbc <- addPerCellQC(sce.pbmc, subsets = list(MT = is.mito))
discard.mito <-
    isOutlier(sce.pmbc$subsets_MT_percent, type = "higher")
plot(
    sce.pmbc$sum,
    sce.pmbc$subsets_MT_percent,
    log = "x",
    xlab = "Total count",
    ylab = "Mitochondrial %"
)
abline(h = attr(discard.mito, "thresholds")["higher"], col = "red")

##########################################################################################
                                     #EXCERCISES
##########################################################################################

#Why does emptyDrops() return NA values?
#la matriz debe ser del mismo tamaño, las rows que no se usAn para el set de empty droplets 
#se les da valores de NA 
#Tambien puede ser porque no tienen datos

#Are the p-values the same for e.out and all.out?
summary(e.out$PValue)
summary(all.out$PValue)
#Con esas instrucciones podemos ver que no son identicos
  
#What if you subset to the non-NA entries?
#mientras usemos la misma semilla forzamos a que nos de el mismo resultado porque las simulaciones 
#estan forzadas a ser las mismas 

## ----marking, cache=TRUE, dependson='use_case'-----------------------------------------------------------------------
# Removing low-quality cells
# Keeping the columns we DON'T want to discard
filtered <- sce.416b[,!discard2]
# Marking low-quality cells
marked <- sce.416b
marked$discard <- discard2



## ----'reproducibility', cache = TRUE, dependson=knitr::all_labels()--------------------------------------------------
options(width = 120)
sessioninfo::session_info()








