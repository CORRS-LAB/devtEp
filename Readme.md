# TODO

1. This project aims to identify cell type specific DNA methylation regions associated with phenotypes using RRBS data.
   To derive a robust estimation of immune cell fractions, we first perform a benchmarking study to compare 19 deconvolution algorithms using RRBS data.
   Four key steps are included in this project: marker selection, benchmarking of deconvolution algorithms, estimating the immune cell fractions from real RRBS data, identifing cell type specific DNA methylation regions associated with phenotypes.

2. Immune Cell Type Specific DNA Methylation Regions associate with phenotypes
   
3. HELP Manual:
    - marker selection and preparing the reference signatures
    - benchmarking of deconvolution algorithms
    - estimating the immune cell fractions from real RRBS data
    - identifying cell type specific DNA methylation regions associated with phenotypes


# Demo: marker selection and preparing the reference signatures.

Regions were selected as deconvolution makers using the algorithm described by Luo et al.[PMID: 22110609; PMID: 23284283]. 
To increase the computational efficiency, we optimized the original marker selection algorithm using C/C++ language and incorporated it into our R package “devtEp”. 
For each region, two-sample t-tests were applied to assess methylation difference between the target cell type and all other cell types. 
The top 200 regions with lowest P values were selected for each cell type, and from these, the 100 regions with the highest mean methylation differences were retained. 
After removing duplicated regions across cell types, the final reference signature was derived.




```{r}
library(devtEp)
library(plyr)

## A function to select top N values in a data matrix, return the row ids.
select_top_N <- function(dt.input,N, dec){

   # Find the top N values in each column
   top_N_indices <- apply(dt.input, 2, function(x) order(x, decreasing = dec)[1:N])

   # Select the rows with the top N values in each column
   selected_rows <- data.matrix(unique(unlist(top_N_indices)))
   return(selected_rows)
}


data <- demo_devEpt_1k

names <- colnames(data)  
colnames(data)[-c(1:3)] <- paste(names[-c(1:3)],".Z00R", sep = "") #### TO DO: Can we remove this?
uu <- selectSite(df = data, skip.col = 3)

n <- 100
row.ids <- c(unlist(select_top_N(uu$pvalue,2*n, FALSE)))
row.ids <- row.ids[order(row.ids)]

Dif <- uu$abs_diff[row.ids,]
rownames(Dif) <- row.ids

top.dif <- apply(Dif, 2, function(x) order(x, decreasing = TRUE)[1:n])
ids <- unique(as.numeric(rownames(Dif[top.dif,])))
row.out <- ids[order(ids)]

out <- data[row.out,]
colnames(out) <- names
write.table(out, "reference.mker.ref",row.names = F, sep = "\t", quote = F)

```


# Demo: benchmarking of deconvolution algorithms


# Demo: estimating the immune cell fractions from real RRBS data


# Demo: identifying cell type specific DNA methylation regions associated with phenotypes

# Performance

# Logs


- 2024-11-03: (TODO) **Manual Update**
- 2024-10-29: **Project initialization**. The first vesiion was finished in July.
