## ImmGen miRNA qPCR normalization and batch correction ##

# load libraries
library(dplyr)
library(mice)
library(readxl)
library(sva)

#### load functions ####
## arithmetic mean normalization
dNormAMbulk <- function(exp){
  
  am <- apply(exp, 2, function(x){ mean(x, na.rm = T)})
  med <- median(am, na.rm = T)
  nf <- med/am
  print(nf)
  exp.n <- t(t(exp) * nf)
  colnames(exp.n) <- colnames(exp)
  return(exp.n)
}
## arithmetic mean normalization by plate
dNormAMplateBulk <- function(exp, plate_index){
  sp <- split(names(plate_index), factor(plate_index))
  d <- lapply(sp, function(p) {
    dNormAMbulk(exp[p,])
  })
  r <- do.call(rbind, d)
  return(r)
}


#### data processing ####
d <- read_xls("~/Dropbox/Data/geo/Immgen_miRNA_GEO_meta.xls",  sheet = 2)
qpcr <- as.matrix(d[,-1])
# convert columns to numeric
qpcr <- apply(qpcr, 2, as.numeric)
# set rownames
rownames(qpcr) <- d$...1


# fix names temporarily to not indexing
colnames(qpcr) <- gsub("\\+", "plus", colnames(qpcr))
colnames(qpcr) <- gsub("\\-", "minus", colnames(qpcr))

# specify outlier assays (Supplementary Table 5)
outlierAssays <- read_xlsx("~/Dropbox/Immgen miRNA Atlas Manuscript/final_for_journal/supp_info/SupplementaryTables.xlsx", 
                           sheet = 5, skip = 1, 
                           col_types = c("guess", "guess", "skip", "skip", "skip"), 
                           n_max = 32)
outlierAssays$Sample <- gsub("\\+", "plus", outlierAssays$Sample)
outlierAssays$Sample <- gsub("\\-", "minus", outlierAssays$Sample)

# set outlier assay measurements to NA
for(i in seq(from = 1, to = nrow(outlierAssays))){
  print(i)
  print(unlist(outlierAssays[i,'miRNA']))
  print(unlist(outlierAssays[i,'Sample']))
  qpcr[unlist(outlierAssays[i,'miRNA']), unlist(outlierAssays[i,'Sample'])] <- NA
}

# load sample metadata from GEO spreadsheet
pDat <- read_xls("~/Dropbox/Data/geo/Immgen_miRNA_GEO_meta_20210413.xls", sheet = 1, 
                        skip = 10, n_max = 141) %>%
  mutate(sample_name.mod = gsub("\\+", "plus", `Sample name`)) %>%
  mutate(sample_name.mod = gsub("\\-", "minus", sample_name.mod)) %>%
  mutate(cell_type = gsub("_[0-9]$", "", sample_name.mod)) 


# load assay metadata from GEO spreadsheet
assayDat <- read_xls("~/Dropbox/Data/geo/Immgen_miRNA_GEO_meta.xls", 
                     sheet = 3, 
                     col_types = c(rep("guess", 8), 
                                   rep("skip", 141)))
# arithmetic mean normalization by plate
plate_idx <- assayDat$plate_no
names(plate_idx) <- rownames(qpcr)

set.seed(1234)
# arithmetic mean normalization by plate
exp.am <- dNormAMplateBulk(qpcr, plate_index = plate_idx)
# input missing values using mice
## need to specify use.matcher = T to replicate processing done previously
exp.am <- mice::complete(mice::mice(exp.am, m = 10, use.matcher = T))
rownames(exp.am) <- rownames(qpcr)

#### batch correction ####
# load batch associated miRs from supplemental table
batch_mirs <- read_xlsx("~/Dropbox/Immgen miRNA Atlas Manuscript/final_for_journal/supp_info/SupplementaryTables.xlsx", 
                        sheet = 5, skip = 2, 
                        col_types = c("skip", "skip", "skip", "skip", "guess"), 
                        col_names = F) %>% unlist()

# filter expression for only miRNAs that will be batch corrected
exp.am.b <- exp.am[batch_mirs,pDat$sample_name.mod]
# create cell type model matrix to preserve this information
mod <- model.matrix(~ pDat$cell_type)
# ComBat batch correction
cleandat.am <- ComBat(as.matrix(exp.am.b), pDat$`characteristics:processing batch`, mod)
# combine uncorrected and corrected data
cleandat.am <- rbind(cleandat.am, exp.am[!(rownames(exp.am) %in% batch_mirs),])[rownames(exp.am),]


#### scale transform ####
# a Ct of 28 as 1 and then scaling from there.
AQtransform <- function(exp, cutoff = 28){
  exp.a <- (exp - cutoff) * -1
  exp.aq <- 2^exp.a
  return(exp.aq)
}
# convert to linear scale
exp.l <- AQtransform(cleandat.am)
# add a pseudocount and convert to log2 scale
exp.l <- exp.l + .1
exp.l2 <- log2(exp.l)

# following this, names can be converted back to their original form

