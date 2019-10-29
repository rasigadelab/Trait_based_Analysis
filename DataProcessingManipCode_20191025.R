######################################################################
### October 25, 2019 
### Clean code for the manipulation of data for trait-based analysis
### Includes all code up until October 25, 2019
### This code manipulates the data in preparation for RLQ and
### Fourth corner analyses
######################################################################

#Load libraries
library(dplyr)
library(tibble)
library(tidyverse)

################################################################
##### Load the data
################################################################
##Load data extracted by Melissa
load('19-10-11_ecoevo_matrices_18ATBList.Rdata')

#'Dataset' contains two lists. 
#List 1: Site, Trait, and Trait x Site
#List 2: Taxonomy, Ward, and Antibiotic

#Site ID's and floras
View(dataset[[1]]$site)

#Full data set (VITEK_ATBList has many NA's due to non testing)
View(dataset[[1]]$trait)#Full data set (VITEK_ATBList has many NA's due to non testing)

#Species (taxon + profile) x Wards
View(dataset[[1]]$trait_by_site) 

#List 2
View(dataset[[2]]$taxonomy)
View(dataset[[2]]$ward)
View(dataset[[2]]$antibiotic)

################################################################
### Manipulation of the data
################################################################

################################################################
### Begin with the manipulation of the Species X Site table
################################################################

#Name the table
sp.site <- (dataset[[1]]$trait_by_site)
View(sp.site)

#Convert first column, containing the Species + Resistance Profile to the rownames
sp.site2 <- column_to_rownames(sp.site, var = "id_trait")
View(sp.site2)

#Simplify the column names (wards) using substring command
names(sp.site2) <- substring(names(sp.site2), 51)

#Transpose the table: Species are in columns and sites are in rows
sp.site3 <- as.data.frame(t(sp.site2))
View(sp.site3)

#Filter the rows for data from urine samples
sp.site.urine <- sp.site3 %>%
  rownames_to_column(var="rowname")%>%
  filter(grepl("Urine",rowname))%>%
  column_to_rownames(var="rowname")

View(sp.site.urine)
#NOTE: this includes all species from original data set, including those that
#have ZERO occurrences in urine, since we're only filtering by rowname

#Remove species with zero occurrences in urine samples
#In order to do so, sum each column then create a list
#of all taxa with a column sum > 0.
#The full data set is then filtered using this list of species
sp.sums.urine <- as.data.frame(colSums(sp.site.urine))
colnames(sp.sums.urine) <- "sum"

#Identify species present in urine
sp.site.urine.pres <- sp.sums.urine %>%
  rownames_to_column(var="rowname")%>%
  filter(sum > 0)%>%
  column_to_rownames(var="rowname")

#Make a vector of the species that have  > 0 occurences in urine
urine.sp.list <- as.vector(rownames(sp.site.urine.pres))

#Keep only species with > 0 occurrences in urine
sp.site.urine2 <- sp.site.urine[colnames(sp.site.urine) %in% urine.sp.list]

#################################################################
####### Ward data        ########################################
#################################################################

#For initial / preliminary data analysis, use ward data from
#metapopulation project.
load("abundance190507.Rdata")
head(abundance)

#Select relevant variables from metapopulation data and 
#convert ward_id to an integer
ward.dat <- abundance%>% select(ward, starts_with("ddd_"))%>%
  rename(ward_id=ward) %>%
  mutate(ward_id = as.integer(ward_id))

#Manipulate the ward data extracted by Mélissa
ward.mel <- dataset$annex_data$ward

#Join the metapopulation ward data to the data extracted by Mélissa
ward.test <- right_join(ward.mel, ward.dat, by="ward_id")

head(ward.test)

#Eliminate NA's from the joined data
ward.test <- na.omit(ward.test)

#Dataset ward.test has several duplicate rows.
#For example: Ward 21124 has two de_labels: MEDECINE INTERNE and 
#SPECIALITES MEDICALES INDIFFERENCIEES. 
#We will need to decide how to address this.

#To test for duplicates:
any(duplicated(ward.test))

#Eliminate de_label and filter to leave only unique wards.
ward.test <- ward.test %>% 
  select(-de_label)%>%
  distinct()

NROW(ward.test)
#168 wards left

#Format ward labels to match labels in Species X Site table
ward.test.urine <- ward.test

#Make a new ward ID column by pasting together "UF", 
#the ward number in the column ward_id, and "|Urine"
ward.test.urine$ward_id_flore <- paste0("UF",ward.test.urine$ward_id, "|Urine")

#Select the relevant ward variables
#Eliminate (using select function + "-" in front of columnn name)
#now irrelevant columns (ward_uri and the original ward_id column)
ward.test.urine2 <- ward.test.urine %>%
  select(bed_number,hg_label,ddd_total, ddd_carba, 
         ddd_c3g_classic,ddd_c3g_pyo, ddd_nsp,ward_id_flore,
         -ward_uri, -ward_id)%>%
  mutate_if(is.character, as.factor)%>%
  column_to_rownames(var="ward_id_flore")


#######################################################################
#### Trait data          ##############################################
#######################################################################

#Define the trait table, starting with the matrix extracted by Mélissa
traits <- dataset[[1]]$trait

#Change the order of the columns in the trait table for readability
traits.ordered <- select(traits, id_trait, species, genus, group, 
                         morpho, everything())
head(traits.ordered)
View(traits.ordered)

#Filter the traits to include only the species found in urine samples
#List of species found in urine samples:
urine.sp.list
#This list was created from filtering the species X site table

traits.urine <- traits.ordered %>%
  filter(id_trait %in% urine.sp.list)

#Select traits.
#In this case I am selecting traits that match well with the 
#current ward data that is available.
#After initial work, it is clear that taxonomic variables
#(species, genus, group, etc.) should NOT be included
#as a trait. Traits should be characteristics across taxa.
#IMPORTANT NOTE: Here we keep "species" in order to simplify the 
#variant names below.
traits.urine2 <- select(traits.urine, id_trait,species,
                        IPM, CTX, FEP, AMC)

#Keep only the data for which there are no NA's
traits.urine3 <- na.omit(traits.urine2)

#Trait table should be a dataframe in order to conduct further manipulations
traits.urine3.df <- as.data.frame(traits.urine3)

#The column id_trait should become the rownames
traits.urine4 <- traits.urine3.df %>%
  mutate_if(is.character,as.factor)%>%
  column_to_rownames(var="id_trait")


#############################################################
### Re-Filtering the Species and Trait Data       ###########
#############################################################
#The trait data has been filtered to only include traits
#that are present in urine samples.
#However, we have also now filtered the species in the trait data (e.g. no NA's)
#Therefore it is now necessary to filter the species data 
#based on the remaining traits

#Begin with the species data:
View(sp.site.urine2)

#Make a vector of the species remaining in the trait data (traits.urine4)
#This will be used to filter the species data
sp.list.urine <- as.vector(rownames(traits.urine4))

#Filter the species data using the above vector
sp.site.urine3 <- sp.site.urine2[colnames(sp.site.urine2) %in% sp.list.urine]

#Reduce the number of wards in the species x site table to match the ward data
#Make a vector of the wards for which we currently have data
ward.list <- as.vector(rownames(ward.test.urine2))

#Keep only species with > 0 occurrences in urine
sp.site.urine4 <- sp.site.urine3[rownames(sp.site.urine3) %in% ward.list,]


#################################################################################
#################################################################################
#################################################################################

#################################################################################
#"Final" data sets for RLQ and Fourth corner analysis for urine samples
#################################################################################
#Species X Site table:
sp.site.urine4

#Environmental data
ward.test.urine2

#Trait data
traits.urine4

#################################################################################
##### Editing of the variant names to improve readability of figures
#################################################################################
#Using the full trait and resistance profile as the species names
#makes the figures unreadable.
#The following code edits the variant names to include only species + a number
View(traits.urine4)

#The following code groups the variants by species and then gives each
#variant within a species a unique number to replace the full resistance profile
#currently in the name.
traits.urine5 <- traits.urine4 %>% 
  rownames_to_column(var="rowname") %>%
  group_by(species) %>% 
  mutate(id = row_number())%>%
  select(species, id, everything())%>%
  unite(new.id, c(species, id), remove=FALSE)%>%
  ungroup()
  
#The new variant names
new.urine.names <- select(traits.urine5, rowname, new.id,species)

#Now look at the species data
View(sp.site.urine4)

#Several manipulations of the species x site table are needed
#First transpose dataframe and make the rownames into a column
#Join the new.urine.names to the dataframe and make the "new.id" into the rownames. 
sp.site.urine5 <- as.data.frame(t(sp.site.urine4))

sp.site.urine6 <- rownames_to_column(sp.site.urine5,var="rowname")

sp.site.urine7 <- full_join(sp.site.urine6, new.urine.names, by="rowname")

#Eliminate the rowname, species, and set "new.id" to the rowname
sp.site.urine8 <- sp.site.urine7%>%
  select(-rowname,-species)%>%
  column_to_rownames(var="new.id")

#Re-transpose the species x site data
sp.site.urine9 <- as.data.frame(t(sp.site.urine8))

#Repeat procedure for the trait data - eliminate rowname, species, and make
#new.id the rowname
traits.urine6 <- traits.urine5 %>%
  select(-rowname,-id,-species)%>%
  mutate_if(is.character, as.factor)%>%
  column_to_rownames(var="new.id")

###########################################################################
#### Final data with simplified variant names
sp.site.urine9

traits.urine6

ward.test.urine2
