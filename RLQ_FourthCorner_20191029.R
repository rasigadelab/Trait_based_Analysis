# Code for RLQ and Fourth Corner Analysis
# Based on Kleyer et al. 2012 and Dray et al. 2014

#Load libraries
library(ade4)

#Data sets:
#These data sets are URINE SAMPLES with simplified names
#Other data options can be seen in the DataProcessingManipCode_20191025 script
# L table: Species x Sites Table
sp.site.urine9

# R table: Environmental variables
ward.test.urine2

# Q table: Species X Trait table
traits.urine6

###################################################################################
#RLQ and Fourth corner analysis
#Start with the ordinations
#Step 1: Correspondence analysis (ordination) for species
afcL.lyon <- dudi.coa(sp.site.urine9, scannf = FALSE)

#Ordination of environment variables, weighted from species ordination above
acpR.lyon <- dudi.hillsmith(ward.test.urine2, row.w = afcL.lyon$lw, #lw refers to weights from species ordination
                            scannf = FALSE)

#Ordination on trait, also weighted
acpQ.lyon <- dudi.acm(traits.urine6, row.w = afcL.lyon$cw,
                      scannf = FALSE)


#Now the rlq analysis that links the three ordinations
rlq.lyon <- rlq(acpR.lyon, afcL.lyon , acpQ.lyon ,
                scannf = FALSE)
plot(rlq.lyon)
svg("rlq.testRtraitsOnly.svg",width = 8, height = 8)
plot(rlq.lyon)
dev.off()


#Separate plots
par(mfrow = c(1, 3))
s.arrow(rlq.lyon$l1)
s.arrow(rlq.lyon$c1)
s.label(rlq.lyon$lQ, boxes = TRUE)

svg("envplotRtraitsOnly.svg",width=6, height=6)
s.arrow(rlq.lyon$l1)
dev.off()

svg("qplot2RtraitsOnly.svg",width = 6, height=6)
s.arrow(rlq.lyon$c1)
dev.off()

svg("species2RtraitsOnly.svg", width=8, height=8)
s.label(rlq.lyon$lQ, boxes = TRUE)
dev.off()

#To compare the structure of the trait (traits PCA), 
#the structure of the environment (Hill-Smith analysis of the environmental variables) 
#and the correlation (correspondence analysis of the sites-species table),
#look at summary of RLQ:
summary(rlq.lyon)

#Test the significance of the RLQ by permutating
#Number of repetitions for permutations:
#NOTE: If permutation number is too high, may crash R / computer.
#On laptops, recommended to test with smaller number of permutations.
nrepet <- 49999
nrepet <- 9999
testrlq.lyon <- randtest(rlq.lyon, modeltype = 6, nrepet = nrepet)

#Check significance
testrlq.lyon
#Is significant

plot(testrlq.lyon)

#Fourthcorner2 function offers a multivariate statistic 
#(equal to the sum of eigenvalues of RLQ analysis).
#Measures the link between two variables by a square correlation coefficient (quant/quant), a Chi2/sum(L) (qual/qual) and a correlation ratio (quant/qual).
Srlq <- fourthcorner2(ward.test.urine2, sp.site.urine9, traits.urine6,
                      modeltype = 6, p.adjust.method.G = "fdr", nrepet = nrepet)
Srlq$trRLQ

#Run Fourth corner analysis alone. Will be used further down in script to 
#plot significant environment - trait associations
nrepet <- 200 #Lowering the permutations to avoid crashing laptop
four.comb.lyon <- fourthcorner(ward.test.urine2, sp.site.urine9, traits.urine6, 
                               modeltype = 6, p.adjust.method.G = "none",
                               p.adjust.method.D = "none", nrepet = nrepet)

#Produces statistics D2, D, and G that test different relationships
#D2 = correlation coefficients  between the quantitative variable and each category separately
#G = association between the quantitative variable and the whole categorical variable
#D = association is estimated between the quantitative variable and each category separately
#by a measure of the within-group homogeneity.
#Tutorial uses the D2 stat

svg("fc.unadjRtraitsOnly.svg", width=6, height=6)
plot(four.comb.lyon, alpha = 0.05, stat = "D2")
dev.off()

#Adjusted p-values for multiple comparisons 
four.comb.lyon.adj <- p.adjust.4thcorner(four.comb.lyon,
                                         p.adjust.method.G = "fdr", p.adjust.method.D = "fdr")


#NOTE: adjusted p-values can be obtained DIRECTLY using the fourthcorner function:
fourthcorner(ward.test.urine2, sp.site.urine3, traits.urine6, modeltype = 6,
             p.adjust.method.G = "fdr", p.adjust.method.D = "fdr",
             nrepet = nrepet)

#Plot with adjusted p values, this reduces the number of significant associations
plot(four.comb.lyon.adj, alpha = 0.05, stat = "D2")

svg("adjfourthcornerRtraitsOnly.svg",width=6, height=6)
plot(four.comb.lyon.adj, alpha = 0.05, stat = "D2")
dev.off()

#Combine RLQ and fourth corner by using RLQ scores to represent traits and env vars
#in biplot. Represent significant associations with line segments
#Here blue = negative, red = positive. Only shows traits with at least one signif association
#P values are adjusted for multiple comparisons, signif = 0.05
plot(four.comb.lyon.adj, x.rlq = rlq.lyon, alpha = 0.05,
     stat = "D2", type = "biplot")


#Directly test link between RLQ axes and traits 
#(typetest="Q.axes") or environmental variables (typetest="R.axes").
testQaxes.comb.lyon <- fourthcorner.rlq(rlq.lyon, modeltype = 6,
                                        typetest = "Q.axes", nrepet = nrepet, p.adjust.method.G = "fdr",
                                        p.adjust.method.D = "fdr")
testRaxes.comb.lyon <- fourthcorner.rlq(rlq.lyon, modeltype = 6,
                                        typetest = "R.axes", nrepet = nrepet, p.adjust.method.G = "fdr",
                                        p.adjust.method.D = "fdr")

print(testQaxes.comb.lyon, stat = "D")


par(mfrow = c(1, 2))
par(mfrow = c(1, 1))
plot(testQaxes.comb.lyon, alpha = 0.05, type = "table",
     stat = "D2")
plot(testRaxes.comb.lyon, alpha = 0.05, type = "table",
     stat = "D2")


svg("QtableRtraitsOnly.svg",width=3, height=6)
plot(testQaxes.comb.lyon, alpha = 0.05, type = "table",
     stat = "D2")
dev.off()

svg("RtableRtraitsOnly.svg",width=10, height=4)
plot(testRaxes.comb.lyon, alpha = 0.05, type = "table",
     stat = "D2")
dev.off()

#Plot factorial map of RLQ analysis. Here, significant
#associations with the first axis in blue, second axis in orange, 
# both axes in green, non-significant in black
par(mfrow = c(1, 2))

plot(testQaxes.comb.lyon, alpha = 0.05, type = "biplot",
     stat = "D2", col = c("black", "blue", "orange", "green"))
plot(testRaxes.comb.lyon, alpha = 0.05, type = "biplot",
     stat = "D2", col = c("black", "blue", "orange", "green"))

svg("testQaxesRtraitsOnly.biplot.svg",width=5, height=5)
plot(testQaxes.comb.lyon, alpha = 0.05, type = "biplot",
     stat = "D2", col = c("black", "blue", "orange", "green"))
dev.off()

svg("testRaxesRtraitsOnly.biplot.svg",width=7, height=7)
plot(testRaxes.comb.lyon, alpha = 0.05, type = "biplot",
     stat = "D2", col = c("black", "blue", "orange", "green"))
dev.off()
