## 10000 Replications
## Import traits data
## Test Phylogenetic signal
## PGLM Model
## Table and Graphies


setwd("~/Desktop/TTU/Johnson Lab/Seed Traits and Phylogeny/R script")
# Part I: Data Preparation

## DataSet1: Germination data
library(readxl)  ## Learn from stackoverflow.com
library(dplyr)
library(data.table) ## Learn from stackoverflow.com 
library(stringr)
library(ggplot2)

## read in species code
GlobalSpeciesCode <- read_excel("Species(Name and Code).xlsx")

CoxLabSPCODE <- read_xlsx("Cox_species numbers and identities.xlsx")
CoxLabSPCODE$spp <-str_glue_data(CoxLabSPCODE, "spp.{`spp number`}")
Species <- c("Prosopis_glandulosa","Astragalus_crassicarpus", 
             "Digitaria_ciliaris","Coreopsis_tinctoria", "Salvia_azurea",
             "Prosopis?", "Monarda_citriodora","Salvia_reflexa",
             "??","Gutierrezia_sarothrae")
CoxLabSPCode <-CoxLabSPCODE %>% mutate(`spp number` = NULL, 
                                       SpeciesName = Species)
## add species function
AddSpeciesCode <- function(SpeciesCode, ControlData) {
  return(merge(SpeciesCode, ControlData, all = TRUE) %>% 
           na.omit())
}

## function extract data from Cox Lab
ExtraControl.Cox <- function(path) {
  species <- excel_sheets(path)
  read_excel_allsheets <- function(filename, tibble = FALSE) {
    sheets <- excel_sheets(filename)
    x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
    if(!tibble) x <- lapply(x, as.data.frame)
    names(x) <- sheets
    x
  }
  CoxLabG <- read_excel_allsheets(path)
  ControlFilter <- function(x){
    return(filter(CoxLabG[[x]], trt == "S0NH") %>%
             select(trt, SM, heat, rep, GC)) %>%
      mutate(spp = species[x])
  }
  i <- seq_along(species)
  ControlData <- lapply(i, ControlFilter) %>%
    rbindlist(fill = TRUE) %>%
    mutate(source = path)
}

CoxYiFang <- ExtraControl.Cox("Cox_Original data_ex1_YiFang_reformat.xlsx")
Coxdata <- ExtraControl.Cox("Cox_data_ex2.xlsx")
GList1 <- rbind(CoxYiFang, Coxdata) 
GList2 <- AddSpeciesCode(CoxLabSPCode, GList1)
GListCox <- AddSpeciesCode(GlobalSpeciesCode, GList2) %>%
  mutate(trt = NULL, SM = NULL, heat = NULL, rep = NULL, 
         Identity = NULL, source = NULL, Family = NULL,
         spp = NULL)

## Schwilk Lab Germination Data
ExtractControl.Schwilk <- function(path) {
  a <- read.csv(path, skip = 1, header = TRUE, stringsAsFactors = FALSE) %>%
    filter(treat == "c", treat2 == "c") %>%
    mutate(source = path)
  names(a) <- c("treat", "treat2", "spcode", "GC", "source")
  return(a)
}

SchwilkLab1 <- ExtractControl.Schwilk("Schwilk_smoke-dc.csv") 
SchwilkLab2 <- ExtractControl.Schwilk("Schwilk_smoke-wc.csv")
SchwilkLab3 <- rbind(SchwilkLab1,SchwilkLab2)
GListSchwilk <- AddSpeciesCode(GlobalSpeciesCode, SchwilkLab3) %>%
  mutate(treat = NULL, treat2 = NULL, source = NULL, Family = NULL)

## Mean Germination for each species
GListAll <- rbind(GListCox, GListSchwilk) %>%
  group_by(SpeciesName, spcode) %>%
  summarise(MeanGC = mean(GC)) 

#========================
## Seed Mass
## DataSet2: Seed Mass
SeedMass <- read_xlsx("Seed Mass5.xlsx") %>%
  mutate(`30Seeds + dish` = NULL, `100seeds + paper` = NULL, Paper = NULL)
## Seed Shape
SeedSurface <- read_xlsx("Seed surface area2.xlsx") %>%
  group_by(spcode) %>%
  summarise(MeanSurfaceArea = mean(Area)) %>%
  merge(GlobalSpeciesCode, all = TRUE) %>%
  na.omit()

#================
SeedTraits.Mass <- merge(GListAll, SeedMass, all = TRUE)
SeedTraits.Mass.NA <- merge(GListAll, SeedMass)
SeedTraits.Surface <- merge(GListAll, SeedSurface, all = TRUE)
SeedTraits.Surface.NA <- merge(GListAll, SeedSurface)
SeedTraits <- merge(SeedTraits.Mass, SeedTraits.Surface, all = TRUE)
SeedTraits.NA <- merge(SeedTraits.Mass.NA, SeedTraits.Surface.NA)
#========================

## Phylogenetic Tree
library(ape)
library(maps)
library(phytools)
library(readxl)
library(geiger)
library(nlme)

library(vegan)
library(permute)
library(lattice)
library(picante)
library(ade4)
library(adephylo)
library(phylobase)


## Read in a master tree from a publication (Zannel et al. 2014)
treeVascularPlants <- read.tree("Vascular_Plants_rooted.tre")
DesiredSpecies <- c(SeedTraits$SpeciesName)

## Prune the tree according to previous Publications
treeDesiredSpecies <- drop.tip(treeVascularPlants, 
                               setdiff(treeVascularPlants$tip.label, 
                                       DesiredSpecies))

plotTree(treeDesiredSpecies)

## Adding "Astragalus_crassicarpus"
tip1 <- "Astragalus_crassicarpus"
sister1 <- "Prosopis_glandulosa"
tree1 <- bind.tip(treeDesiredSpecies,tip1,
                  edge.length = 0.5*treeDesiredSpecies$edge.length[which(treeDesiredSpecies$edge[,2]==
                                                                           which(treeDesiredSpecies$tip.label==sister1))],
                  where=which(treeDesiredSpecies$tip.label==sister1),
                  position=0.5*treeDesiredSpecies$edge.length[which(treeDesiredSpecies$edge[,2]==
                                                                      which(treeDesiredSpecies$tip.label==sister1))])
plotTree(tree1)

## Adding "Liatris_mucronata"
tip2 <- "Liatris_mucronata"
sister2 <- "Liatris_pycnostachya"
tree2 <- bind.tip(tree1,tip2,
                  edge.length = 0.5*tree1$edge.length[which(tree1$edge[,2]==
                                                              which(tree1$tip.label==sister2))],
                  where=which(tree1$tip.label==sister2),
                  position=0.5*tree1$edge.length[which(tree1$edge[,2]==
                                                         which(tree1$tip.label==sister2))])
plotTree(tree2)


## Adding "Salvia penstemonoide"
## Adding species to genus
tree3 <- add.species.to.genus(tree2, "Salvia_penstemonoide", genus = NULL, where = "random")
tree1$edge.length
plotTree(tree3)

tree4 <- add.species.to.genus(tree3, "Salvia_reflexa", genus = NULL, where = "random")
plotTree(tree4)

tree5 <- add.species.to.genus(tree4, "Bouteloua_eriopoda", genus = NULL, where = "random")
plotTree(tree5)

#======================
## Phylogenetic Independent Contrast(PIC) of germination and seed mass 

Germination.mass <- SeedTraits.Mass.NA$MeanGC
names(Germination.mass) <- SeedTraits.Mass.NA$SpeciesName
treeGerminationSpecies.mass <- drop.tip(tree5, 
                                        setdiff(tree5$tip.label, SeedTraits.Mass.NA$SpeciesName))
pic.germination.mass <- pic(Germination.mass, treeGerminationSpecies.mass)

Mass <- SeedTraits.Mass.NA$`100Seeds`
names(Mass) <- SeedTraits.Mass.NA$SpeciesName
treeMassSpecies <- drop.tip(tree5, 
                            setdiff(tree5$tip.label, SeedTraits.Mass.NA$SpeciesName))
pic.mass <- pic(Mass, treeMassSpecies)

#----------

Germination.surface <-SeedTraits.Surface.NA$MeanGC
names(Germination.surface) <- SeedTraits.Surface.NA$SpeciesName
treeGerminationSpecies.surface <-drop.tip(tree5, 
                                          setdiff(tree5$tip.label, SeedTraits.Surface.NA$SpeciesName))
pic.germination.surface <- pic(Germination.surface, treeGerminationSpecies.surface)

Surface <- SeedTraits.Surface.NA$MeanSurfaceArea
names(Surface) <- SeedTraits.Surface.NA$SpeciesName
treeSurfaceSpecies <- drop.tip(tree5, 
                               setdiff(tree5$tip.label, SeedTraits.Surface.NA$SpeciesName))
pic.surface <- pic(Surface, treeSurfaceSpecies)

#------------

par(mar=c(5,5,5,1))
plot(SeedTraits.Mass.NA$`100Seeds` ~ SeedTraits.Mass.NA$MeanGC,
     xlab = "Seed Mass per 100 / g",
     ylab = "Seed Germination Rate (%)")
abline(lm(SeedTraits.Mass.NA$`100Seeds` ~ SeedTraits.Mass.NA$MeanGC))

plot(pic.mass ~ pic.germination.mass, 
     xlab = "PIC Value of Seed Mass (g)",
     ylab = "PIC Value of Seed Germination Rate (%)")
abline(lm(pic.mass ~ pic.germination.mass))


plot(SeedTraits.Surface.NA$MeanSurfaceArea ~ SeedTraits.Surface.NA$MeanGC, 
     xlab = "Seed Surface Area (pixel)",
     ylab = "Seed Germination Rate (%)")
abline(lm(SeedTraits.Surface.NA$MeanSurfaceArea ~ SeedTraits.Surface.NA$MeanGC))

plot(pic.surface ~ pic.germination.surface, 
     xlab = "PIC Value of Seed Surface Area",
     ylab = "PIC Value of Seed Germination Rate")
abline(lm(pic.surface ~ pic.germination.surface))

#===================
##Phylogenetic Signal

phylosig(treeGerminationSpecies.mass, Mass, 
         method="lambda", test=TRUE, nsim=999)
phylosig(treeGerminationSpecies.mass, Mass, 
         method="K", test=TRUE, nsim=999)


phylosig(treeGerminationSpecies.surface, Surface, 
         method="lambda", test=TRUE, nsim=999)
phylosig(treeGerminationSpecies.surface, Surface, 
         method="K", test=TRUE, nsim=999)


#------------

phylosig(treeGerminationSpecies.surface, Surface, 
         method="lambda", test=TRUE, nsim=999)
phylosig(treeGerminationSpecies.surface, Surface, 
         method="K", test=TRUE, nsim=999)