setwd("D:/CV/Krumsiek_task/")
load("Diseasome_preprocessed.RData")

#explore the diseasome data frame
str(diseasome)

# number of unique diseases 
num_unique_diseases <- length(unique(diseasome$Disease_ID))

# number of unique genes 

num_unique_genes <- length(unique(diseasome$Gene_ID))

#number of edges 
length(diseasome$Edge)

#matrix of the diseasome 
diseasome_matrix <- matrix(nrow = length(diseasome$Disease_ID) , ncol =length(diseasome$Gene_ID) )

dim(diseasome_matrix)

diseasome_matrix <- matrix(0,nrow = 100 , ncol =100)

#fill the matrix by the edges 
for(i in 1:nrow(diseasome[c(1:100),])){
  diseasome_matrix[diseasome[i,1],diseasome[i,2]] <- 1
}


#dimensions of the matrix 
dim(diseasome_matrix)


# Transpose the diseasome_matrix
gene_disease_matrix <- t(diseasome_matrix)
gene_disease_matrix
# Multiply the transposed matrix by the original matrix
gene_gene_matrix <- gene_disease_matrix %*% diseasome_matrix
gene_gene_matrix
# Convert the gene-by-gene matrix to a disease-by-disease matrix
disease_disease_matrix <- t(gene_gene_matrix) %*% gene_gene_matrix

#give the rows and columns names 
rownames(disease_disease_matrix) <- colnames(disease_disease_matrix) <- paste0("D", 1:nrow(disease_disease_matrix))

disease_disease_matrix


# Count the number of edges THEN 
#subtract the number of diagonal values which represented by the number of rows 
num_edges <- sum(disease_disease_matrix > 0) - nrow(disease_disease_matrix)

# Multiply by 2 to account for symmetry
total_edges <- 2 * num_edges

total_edges
