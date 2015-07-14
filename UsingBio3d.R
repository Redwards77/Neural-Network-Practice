#---
#title: "Downloading Protein Data Bank (PDB) Info"
#author: "Paul Pearson"
#date: "July 2, 2015"
#output: html_document
#---
  
# Download and read amino acid sequences from PDB using bio3d library
library(bio3d)

# Parsing the mutation string

parse_mutation <- function(mutation) {
  # search: do a greedy match for a digit [0-9] 
  # followed by any number of characters .* 
  # up to the end of the line $.  
  # replace: empty
  before_mutation <- gsub("[0-9].*$","", mutation)
  # do a global substitution for any alpha character [A-Z] 
  # and replace it by the empty string.
  mutation_character_position <- gsub("[A-Z]", "", mutation)
  # search: do a greedy match for the beginning of the line ^ 
  # followed by any number of characters .* up to a digit [0-9].  
  # replace: empty
  after_mutation <- gsub("^.*[0-9]","", mutation)
  
  return( list(before = as.character(before_mutation), 
               position = as.integer(mutation_character_position), 
               after = as.character(after_mutation) ) )
}

# Making a feature vector from a mutation

mutate_feature <- function(mutation) {
  
  mut <- parse_mutation(mutation)
  # character vector of 20 amino acids
  #AA <- sort(c("A","C","G","T","U","R","Y","K","M","S","W","B","D","H","V","N","X"))
  AA <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
  feature_vector <- rep(0, times=length(AA))
  feature_vector[ which( AA == mut$before ) ] <- -1
  feature_vector[ which( AA == mut$after  ) ] <-  1
  return(feature_vector)
}

# Get the sequence and mutant strings

mutant_seq_string <- function(PDB_ID, mutation) {
  
  mut <- parse_mutation(mutation)
  aa_seq <- pdbseq(read.pdb(PDB_ID))
  # Note: by using mut_index <- which()[1] below, 
  # we're just grabbing the first match for convenience
  # and this choice could be a potential source for error.
  mut_index <- which( names(aa_seq) == mut$position )[1]
  
  aa_seq <- paste0( as.vector(aa_seq), collapse="")
  
  mut_string <- paste0(
    substr(aa_seq, 1, mut_index-1), 
    mut$after, 
    substr(aa_seq, mut_index+1, nchar(aa_seq)) )
  
  return(mut_string)
  
}

# Get solvent accessibility

solventAcc <- function(PDB_ID, mutation) {
  
  pdb <- read.pdb(PDB_ID)
  seq <- pdbseq(pdb)
  acc <- bio3d::dssp(pdb, exepath = dsspExePath )$acc
  mut <- parse_mutation(mutation)
  mut_location <- as.integer(mut[2])
  
  return( acc[ which( names(seq) == mut_location ) ] )
  
}
