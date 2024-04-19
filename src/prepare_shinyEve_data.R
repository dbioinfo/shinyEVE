### combine all the files into one dataframe and save it as a csv file
library(tidyverse)

#set working directory
setwd("~/Postdoc/Hindle/broad/")

#import EVE results
thetafiles <- list.files('data/shinyEveData', pattern='twoTheta')
thetas <- tibble(filename=thetafiles) %>% 
  mutate(file = map(filename, ~read_delim(paste0('data/shinyEveData/',.x),'\t'))) %>%  #add in local path here
  unnest(file) %>%
  separate(filename, sep = '-', into = c('trash','datalevel','clade','temp')) %>% #split the filename to determine params
  mutate(clade = case_when(
    clade == 'allflexible' ~ 'flexible',
    clade == 'hibernators' ~ 'hibernators',
  )) %>% #fix clade names
  mutate(temp = gsub(pattern = "\\.tsv$", "", temp)) %>%
  select(-trash) %>% 
  rename(gene = 'OG')
  


#import expression data, prune out non temp experiments
banned_samps <- c('human_HS2096_37C_2-5mM',
                   'human_HS2096_37C_30mM',
                   'human_HS3099_37C_2-5mM',
                   'human_HS3099_37C_30mM',
                   'human_HS3896_37C_2-5mM',
                   'human_HS3896_37C_30mM',
                   'human_HS4057_37C_2-5mM',
                   'human_HS4057_37C_30mM',
                   'human_HS967_37C_2-5mM',
                   'human_HS967_37C_30mM')
expr <- read_delim('data/shinyEveData/TPM_9-species.tsv','\t') %>% #add in local path here
  rename(gene='Gene') %>% 
  select(-any_of(banned_samps) )


#import lfc data
lfcfiles <- list.files('data/shinyEveData', pattern='Gene-score')
lfc <- tibble(filename=lfcfiles) %>% 
  mutate(file = map(filename, ~read_delim(paste0('data/shinyEveData/',.x),'\t'))) %>%  #add in local path here
  unnest(file) %>%
  separate(filename, sep = '-', into = c('trash','trash1','trash2','temp')) %>% #split the filename to determine params
  select(-c(trash, trash1, trash2)) %>% 
  mutate(temp = gsub(pattern = "\\.tsv$", "", temp)) %>%
  relocate(gene)

#filter out genes that were not in theta
theta_genes <- thetas %>% pull(gene) %>% unique()
expr_genes <- expr %>% pull(gene) %>% unique()
lfc_genes <- lfc %>% pull(gene) %>% unique()
final_genes <- theta_genes  #intersect(intersect(theta_genes, expr_genes), lfc_genes) #only genes in all 3

thetas <- thetas %>% filter(gene %in% final_genes)
expr <- expr %>% filter(gene %in% final_genes)
lfc <- lfc %>% filter(gene %in% final_genes)

#see how many genes are common among all 3 datasets, note how few fall in the middle under final genes
library(ggvenn)
ggvenn(list(theta=theta_genes, expr=expr_genes, lfc=lfc_genes))

#get gene symbols from ensembl symbols
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes <- getBM(attributes = c('ensembl_gene_id','hgnc_symbol'), mart = ensembl) %>% 
  mutate(clean_gname = case_when(
    hgnc_symbol == "" ~ ensembl_gene_id,
    T ~ hgnc_symbol
  )) %>% 
  dplyr::select(-hgnc_symbol)


#store gene names in theta
thetas <- thetas %>% left_join(genes, by = c('gene' = 'ensembl_gene_id'))


#save the data
write_csv(thetas, 'data/shinyEveData/cleaned_EVE_thetas.csv')
write_csv(expr, 'data/shinyEveData/cleaned_EVE_expr.csv')
write_csv(lfc, 'data/shinyEveData/cleaned_EVE_lfc.csv')

