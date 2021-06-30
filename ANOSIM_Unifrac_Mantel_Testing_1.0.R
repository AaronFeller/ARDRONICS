## Matrix Comparison with Analysis of Similarities, Unifrac, and Mantel Test ###


### ANOSIM ###

#install.packages("vegan")
library(vegan)

data(dune)
data(dune.env)
dune.dist <- vegdist(dune)
attach(dune.env)
dune.ano <- anosim(dune.dist, Management)
summary(dune.ano)
plot(dune.ano)



### UNIFRAC ###

#install.packages("GUniFrac")
library(GUniFrac)

GUniFrac

data(throat.otu.tab)
data(throat.tree)
data(throat.meta)

groups <- throat.meta$SmokingStatus
# Rarefaction
otu.tab.rff <- Rarefy(throat.otu.tab)$otu.tab.rff
# Calculate the UniFracs
unifracs <- GUniFrac(otu.tab.rff, throat.tree, alpha=c(0, 0.5, 1))$unifracs
dw <- unifracs[, , "d_1"] 
# Weighted UniFrac

du <- unifracs[, , "d_UW"] # Unweighted UniFrac
dv <- unifracs[, , "d_VAW"] # Variance adjusted weighted UniFrac
d0 <- unifracs[, , "d_0"]      # GUniFrac with alpha 0
d5 <- unifracs[, , "d_0.5"]    # GUniFrac with alpha 0.5


# Permanova - Distance based multivariate analysis of variance
adonis(as.dist(d5) ~ groups)

# Combine d(0), d(0.5), d(1) for testing
PermanovaG(unifracs[, , c("d_0", "d_0.5", "d_1")] ~ groups)




### MANTEL TEST ###


# NOT RUN {
q1 <- matrix(runif(36), nrow = 6)
q2 <- matrix(runif(36), nrow = 6)
diag(q1) <- diag(q2) <- 0
mantel.test(q1, q2, graph = TRUE,
            main = "Mantel test: a random example with 6 X 6 matrices
representing asymmetric relationships",
            xlab = "z-statistic", ylab = "Density",
            sub = "The vertical line shows the observed z-statistic")
# }



mantel.rtest()
dist1 <- as.dist(dissm1)
dist2 <- as.dist(dissm2)

plot(r1 <- mantel.randtest(dist1,dist2), main = "Mantel's test")
r1


mantel.test(dissm1, dissm2, graph = FALSE, alternative = "two.sided", ...)


### MANTEL 2 ###
install.packages("ape")
library(ape)

mantel.test(dissm1, dissm2, nperm = 999, graph = T, alternative = "two.sided")
