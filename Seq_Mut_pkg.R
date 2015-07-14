AASeq <- function(PDB_ID){
  library(RCurl)
  #PDB_ID <- "1A23"
  protein_seq <- getURL( paste0("http://www.rcsb.org/pdb/files/fasta.txt?structureIdList=", PDB_ID) )
  protein_seq
  
  # Remove the first line of text and get just the sequence
  # use sub to match only once (instead of gsub)
  # ^ matches the beginning of a line
  # . matches any character and .* matches any character any number of times
  # SEQUENCE\\n matches SEQUENCE followed by a newline
  protein_seq <- sub("^.*SEQUENCE\\n", "", protein_seq)
  return(protein_seq)
}

Mut <- function(AASeq,mutation){
  #mutation <- "H32L"
  # search: do a greedy match for a digit [0-9] followed by any number of characters .* up to the end of the line     $.  replace: empty
  before_mutation <- gsub("[0-9].*$","", mutation)
  # do a global substitution for any alpha character [A-Z] and replace it by the empty string.
  mutation_character_position <- gsub("[A-Z]", "", mutation)
  # search: do a greedy match for the beginning of the line ^ followed by any number of characters .* up to a         digit [0-9].  replace: empty
  after_mutation <- gsub("^.*[0-9]","", mutation)
  
  # In position 32 of the protein sequence, replace H by L
  # Use the substring command to get the first 1:(n-1) positions of protein_seq, followed by the replacement          character after_mutation, followed by positions (n+1):end in the original protein_seq
  n <- as.integer(mutation_character_position)
  mut <- paste(substr(protein_seq, 1, n-1), after_mutation, substr(protein_seq, n+1, nchar(protein_seq)), 
               sep = "")
  return(mut)
}