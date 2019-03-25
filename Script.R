## Set working directory
setwd("~/Desktop/TTU/Johnson Lab/Seed Traits and Phylogeny/R script")

## Preparing Phylogenetic Tree
library(ape)
library(maps)
library(phytools)
library(readxl)

## Prune the tree with Existing Species

## Read in a master tree from a publication (Zannel et al. 2014)
treeVascularPlants <- read.tree("Vascular_Plants_rooted.tre")

## Prune the tree according to previous Publications
GlobalSpeciesCode <- read_excel("Species(Name and Code).xlsx")
DesiredSpecies <- c(GlobalSpeciesCode$SpeciesName)
### WHY it doesn't work
treeDesiredSpecies <- keep.tip(treeVascularPlants, DesiredSpecies) 

treeDesiredSpecies <- drop.tip(treeVascularPlants, 
                               setdiff(treeVascularPlants$tip.label, 
                                       DesiredSpecies))

## Adding New Species_5 species







## 10000 Replications
## Import traits data
## Test Phylogenetic signal
## PGLM Model
## Table and Graphies
