library(mixOmics) # import the mixOmics library
library(BiocParallel)
library(ggplot2)

set.seed(5249)
source("network_test.R")


# Dataset loading

# r x c
data(liver.toxicity) # 64x3116
X <- liver.toxicity$gene # X
Y <- liver.toxicity$clinic
Z <- liver.toxicity$treatment

dim(X)
dim(Y)



# Dimensionality reduction - exploration of data


pca.gene <- pca(X, ncomp = 10, center = TRUE, scale = TRUE)

gene_scores <- as.data.frame(pca.gene$variates$X)
gene_scores$dose_group <- as.factor(Z[, "Dose.Group"])
gene_scores$time_group  <- as.factor(Z[, "Time.Group"]) #time necropsy  

barplot(
  pca.gene$prop_expl_var$X * 100,
  names.arg = paste0("PC", 1:10),
  col = "steelblue",
  las = 2,
  xlab = "Principal Component",
  ylab = "Variance Explained (%)",
  main = "Variance Explained by top 10 Components - Mouse liver gene data"
)


## PC1 variablity within the treatment times among the different groups
## PC2 subdivide majorly the treamtment at 24 hours between group groups of 50&150 to those higher dose 1500 2000mg

plot_gene <- ggplot(gene_scores, aes(x = PC1, y = PC2, colour = time_group)) +
  theme_minimal(base_size = 14) +
  geom_point(size = 3) +
  geom_text(aes(label = paste0(" ",dose_group)), hjust = -0.2, vjust = 0.5, size = 3.5) +
  labs (
    title = "Comp. 1 vs 2 on Mouse liver gene",
    x = paste0("PC1 - (", round(pca.gene$prop_expl_var$X["PC1"] * 100, 2)," %)"),
    y = paste0("PC2 - (", round(pca.gene$prop_expl_var$X["PC2"] * 100, 2)," %)"),
    colour = "Treatment time"
  )


plot_gene



clinic_scores <- as.data.frame(pca.clinic$variates$X)
clinic_scores$dose_group <- as.factor(Z[, "Dose.Group"])
clinic_scores$time_group  <- as.factor(Z[, "Time.Group"])

## PC1 variablity within the treatment times among the different groups
## PC2 subdivide majorly the treamtment at 24 hours between group groups of 50&150 to those higher dose 1500 2000mg

barplot(
  pca.clinic$prop_expl_var$X * 100,
  names.arg = paste0("PC", 1:10),
  col = "steelblue",
  las = 2,
  xlab = "Principal Component",
  ylab = "Variance Explained (%)",
  main = "Variance Explained by top 10 Components - Mouse liver gene data"
)

plot_clinic <- ggplot(clinic_scores, aes(x = PC1, y = PC2, colour = time_group)) +
  theme_minimal(base_size = 14) +
  geom_point(size = 3) +
  geom_text(aes(label = paste0(" ", dose_group)), hjust = -0.2, vjust = 0.5, size = 3.5) +
  labs (
    title = "Comp. 1 vs 2 on Mouse liver clinic",
    x = paste0("PC1 - (", round(pca.clinic$prop_expl_var$X["PC1"] * 100, 2)," %)"),
    y = paste0("PC2 - (", round(pca.clinic$prop_expl_var$X["PC2"] * 100, 2)," %)"),
    colour = "Treatment time"
  )

plot_clinic


###################################################################################################################


# sPLS & model tuning

spls.liver <- spls(X = X, Y = Y, ncomp = 5, mode = 'regression')
cores <- MulticoreParam(workers = 24)
performance.spls.liver <- perf(object = spls.liver, validation = "Mfold", progressBar = FALSE, folds = 5, nrepeat = 50, BPPARAM = cores)

plot(performance.spls.liver, criterion = "Q2.total", title = "Q2 performance assessment for the PLS model", cex = 2)


#list.keepX <- c(seq(20, 50, 5))
list.keepX <- c(seq(5, 30, 5))
list.keepY <- c(3:10) 

tune.spls.liver <- tune.spls(X, Y, ncomp = 2,
                             test.keepX = list.keepX,
                             test.keepY = list.keepY,
                             nrepeat = 50, folds = 5,
                             mode = 'regression', measure = 'cor', BPPARAM = cores) 
plot(tune.spls.liver)  # 20 ~ genes constant in Comp 1 and 3 clinical measurement, which is a strong indication that this is a meaningful sparsity level.

tune.spls.liver$choice.keepX
tune.spls.liver$choice.keepY

optimal.keepX <- tune.spls.liver$choice.keepX 
optimal.keepY <- tune.spls.liver$choice.keepY
optimal.ncomp <-  length(optimal.keepX); print(optimal.ncomp)

###################################################################################################################


final.spls.liver <- spls(X, Y, ncomp = optimal.ncomp, 
                         keepX = optimal.keepX,
                         keepY = optimal.keepY,
                         mode = "regression")

plotIndiv(final.spls.liver, ind.names = FALSE, 
          rep.space = "X-variate", # plot in X-variate subspace
          group = liver.toxicity$treatment$Time.Group, # colour by time group
          pch = as.factor(liver.toxicity$treatment$Dose.Group), 
          col.per.group = color.mixo(1:4), 
          legend = TRUE, legend.title = 'Time', legend.title.pch = 'Dose')

plotIndiv(final.spls.liver, ind.names = FALSE,
          rep.space = "Y-variate", # plot in Y-variate subspace
          group = liver.toxicity$treatment$Time.Group, # colour by time group
          pch = as.factor(liver.toxicity$treatment$Dose.Group), 
          col.per.group = color.mixo(1:4), 
          legend = TRUE, legend.title = 'Time', legend.title.pch = 'Dose')

plotIndiv(final.spls.liver, ind.names = FALSE, 
          rep.space = "XY-variate", # plot in averaged subspace
          group = liver.toxicity$treatment$Time.Group, # colour by time group
          pch = as.factor(liver.toxicity$treatment$Dose.Group), # select symbol
          col.per.group = color.mixo(1:4),                      # by dose group
          legend = TRUE, legend.title = 'Time', legend.title.pch = 'Dose')

plotArrow(final.spls.liver, ind.names = FALSE,
          group = liver.toxicity$treatment$Time.Group, # colour by time group
          col.per.group = color.mixo(1:4), arrow.alpha = 1,
          legend.title = 'Time.Group')

###################################################################################################################

# form new perf() object which utilises the final model
perf.spls.liver <- perf(final.spls.liver, 
                        folds = 5, nrepeat = 50, BPPARAM = cores, # use repeated cross-validation
                        validation = "Mfold", 
                        dist = "max.dist",  # use max.dist measure
                        progressBar = TRUE)
par(mfrow=c(1,2)) 
plot(perf.spls.liver$features$stability.X$comp1, type = 'h',
     ylab = 'Stability',
     xlab = 'Features',
     main = '(a) Comp 1', las =2,
     xlim = c(0, 150))
plot(perf.spls.liver$features$stability.X$comp2, type = 'h',
     ylab = 'Stability',
     xlab = 'Features',
     main = '(b) Comp 2', las =2,
     xlim = c(0, 300))


plotVar(final.spls.liver, cex = c(3,4), var.names = c(FALSE, TRUE))


###################################################################################################################


color.edge <- color.GreenRed(50)  # set the colours of the connecting lines

# X11() # To open a new window for Rstudio
mixOmics::network(final.spls.liver, comp = 1:2, alpha.node = 1, cutoff = 0, # only show connections with a correlation above 0.7
        shape.node = c("rectangle", "circle"), graph.scale = 0.75, keysize = c(2,1), show.edge.labels = T, cex.edge.label = 0.45, 
        color.node = c("cyan", "pink"),
        color.edge = color.edge, size.node = 0.25, cex.node.name = 0.5,
        save = 'png', # save as a png to the current working directory
        name.save = 'sPLS Liver Toxicity Case Study Network Plot')

cim(final.spls.liver, comp = 1:2, xlab = "\nclinic", ylab = "genes", title = "Cluster Image Map of Liver Genes with Liver clinical data", row.cex = 1.2, col.cex = 1.2, center = TRUE, scale = TRUE, save = "pdf")

cim(
  final.spls.liver, color = color.edge,
  comp = 1:2,
  title = "Cluster Image Map of Liver Genes with Liver Clinical Data",
  xlab = "\n\nClinical variables",
  ylab = "Genes",
  margins = c(9, 9),      # <- much larger
  row.cex = 1.2,
  col.cex = 1.1
)
    
