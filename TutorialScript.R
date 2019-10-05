######################################################
### Exploring trait-based analysis using ade4 package
### Based on Kleyer et al. 2012 / Dray et al. 2014
### tutorials

library(ade4)
data(aravo)

#Check dimensions of the aravo data set
dim(aravo$spe)
dim(aravo$traits)

#Step 1: Correspondence analysis (ordination) for species
afcL.aravo <- dudi.coa(aravo$spe, scannf = FALSE)

#Ordination of environment variables, weighted from species ordination above
acpR.aravo <- dudi.hillsmith(aravo$env, row.w = afcL.aravo$lw, #lw refers to weights 
                                                                #from species ordination
                             scannf = FALSE)
#Ordination on trait, also weighted
acpQ.aravo <- dudi.pca(aravo$traits, row.w = afcL.aravo$cw,
                       scannf = FALSE)

#Now the rlq analysis that links the three ordinations
rlq.aravo <- rlq(acpR.aravo, afcL.aravo, acpQ.aravo,
                 scannf = FALSE)
plot(rlq.aravo)
#This plots all the elements, R scores, Q scores, eigenvalues

#To plot each plot separately:
par(mfrow = c(1, 3))
s.arrow(rlq.aravo$l1)
s.arrow(rlq.aravo$c1)
s.label(rlq.aravo$lQ, boxes = FALSE)

#To compare the structure of the trait (traits PCA), 
#the structure of the environment (Hill-Smith analysis of the environmental variables) 
#and the correlation (correspondence analysis of the sites-species table),
#look at summary:
summary(rlq.aravo)

#Next part of tutorial demonstrates Fourth corner analysis on raw data
#which shows strength of association between all env var - trait pairs
#I'm skipping straight to combining fourth corner with RLQ

#Look at significance of the RLQ
#Set number of permutations
nrepet <- 49999
testrlq.aravo <- randtest(rlq.aravo, modeltype = 6, nrepet = nrepet)

#Check significance
testrlq.aravo
#Is significant

plot(testrlq.aravo)

#Total inertia of RLQ analysis = SRLQ multivariate statistic defined in 
#Dray and Legendre (2008). 
#Use fourthcorner2 function:
Srlq <- fourthcorner2(aravo$env, aravo$spe, aravo$traits,
                        modeltype = 6, p.adjust.method.G = "fdr", nrepet = nrepet)
Srlq$trRLQ

#Run Fourth corner analysis alone. Will be used later to plot signif
#env - trait associations
nrepet <- 49999
four.comb.aravo <- fourthcorner(aravo$env, aravo$spe,
                                aravo$traits, modeltype = 6, p.adjust.method.G = "none",
                                p.adjust.method.D = "none", nrepet = nrepet)

#Produces statistics D2, D, and G that test different relationships
#D2 = correlation coefficients  between the quantitative variable and each category separately
#G = association between the quantitative variable and the whole categorical variable
#D = association is estimated between the quantitative variable and each category separately
#by a measure of the within-group homogeneity.
#Tutorial uses the D2 stat

plot(four.comb.aravo, alpha = 0.05, stat = "D2")

#Adjusted p-values for multiple comparisons 
four.comb.aravo.adj <- p.adjust.4thcorner(four.comb.aravo,
                                          p.adjust.method.G = "fdr", p.adjust.method.D = "fdr")

#NOTE: adjusted p-values can be obtained DIRECTLY using the fourthcorner function:
fourthcorner(aravo$env, aravo$spe, aravo$traits, modeltype = 6,
               p.adjust.method.G = "fdr", p.adjust.method.D = "fdr",
               nrepet = nrepet)

#Plot with adjusted p values, this reduces the number of significant associations
plot(four.comb.aravo.adj, alpha = 0.05, stat = "D2")

#Combine RLQ and fourth corner by using RLQ scores to represent traits and env vars
#in biplot. Represent significant associations with line segments
#Here blue = negative, red = positive. Only shows traits with at least one signif association
#P values are adjusted for multiple comparisons, signif = 0.05
plot(four.comb.aravo.adj, x.rlq = rlq.aravo, alpha = 0.05,
     stat = "D2", type = "biplot")


#Directly test link between RLQ axes and traits 
#(typetest="Q.axes") or environmental variables (typetest="R.axes").

testQaxes.comb.aravo <- fourthcorner.rlq(rlq.aravo, modeltype = 6,
                                         typetest = "Q.axes", nrepet = nrepet, p.adjust.method.G = "fdr",
                                         p.adjust.method.D = "fdr")
testRaxes.comb.aravo <- fourthcorner.rlq(rlq.aravo, modeltype = 6,
                                         typetest = "R.axes", nrepet = nrepet, p.adjust.method.G = "fdr",
                                         p.adjust.method.D = "fdr")

print(testQaxes.comb.aravo, stat = "D")

par(mfrow = c(1, 2))
plot(testQaxes.comb.aravo, alpha = 0.05, type = "table",
     stat = "D2")
plot(testRaxes.comb.aravo, alpha = 0.05, type = "table",
     stat = "D2")

#Plot factorial map of RLQ analysis. Here, significant
#associations with the first axis in blue, second axis in orange, 
# both axes in green, non-significant in black
par(mfrow = c(1, 2))

plot(testQaxes.comb.aravo, alpha = 0.05, type = "biplot",
     stat = "D2", col = c("black", "blue", "orange", "green"))
plot(testRaxes.comb.aravo, alpha = 0.05, type = "biplot",
     stat = "D2", col = c("black", "blue", "orange", "green"))
