#---
#  title: "Downloading Protein Data Bank (PDB) Info"
#author: "Paul Pearson"
#date: "July 2, 2015"
#output: html_document
#---
  
  # Download and read amino acid sequences from PDB using bio3d library
  
library(bio3d)

# download protein 1AAR from PDB
pdb_1AAR <- read.pdb("1AAR")
seq_1AAR <- pdbseq(pdb_1AAR)
seq_1AAR
seq_1AAR[77]
class(seq_1AAR) # character vector with names
names(seq_1AAR)
as.vector(seq_1AAR)

pdb_1BEN <- read.pdb("1BEN")
seq_1BEN <- pdbseq(pdb_1BEN)
seq_1BEN

pdb_1AG2 <- read.pdb("1AG2")
seq_1AG2 <- pdbseq(pdb_1AG2)
seq_1AG2
# mutation: E200K
mut_1AG2 <- seq_1AG2
mut_1AG2[ which( names(seq_1AG2) == 200 ) ] <- "K"
mut_1AG2

# Get the frequency of sequence neighbors

sequence_neighbors <- function(PDB_ID,seq_index) {
  
  pdb <- read.pdb(PDB_ID)
  seq_vect <- pdbseq(pdb)   
  # seq_index: indexing of proteins in pdb file, does not necessarily start with 1
  # sequence position: indexing that starts with 1
  seq_pos <- which( names(seq_vect) == seq_index )
  
  if (seq_pos < 16) {
    seq_start <- 1
  } else {
    seq_start <- seq_pos - 15
  }
  
  if (length(seq_vect) - 16 < seq_pos) {
    seq_end <- length(seq_vect)
  } else {
    seq_end <- seq_pos + 15
  }
  
  #return(seq_vect[seq_start:seq_end])
  
  AA <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
  
  #freq_list <- rle(sort( seq_vect[seq_start:seq_end]))
  
  feature_vector <- rep(0,length(AA))
  
  for (i in 1:length(AA)) {
    feature_vector[i] <- sum( seq_vect[seq_start:seq_end] == AA[i] )
  }
  
  return(feature_vector)
}

sequence_neighbors("1AG2",131)

# Get DSSP working with bio3d

# Note: in RStudio you may need to choose "Session -> Set working directory -> To source file location" in order for DSSP to be found.
fileURL <- "ftp://ftp.cmbi.ru.nl/pub/software/dssp/dssp-2.0.4-win32.exe"
fileName <- "dssp-2.0.4-win32.exe"
if (!file.exists(fileName)) {
  download.file( url = fileURL, destfile = fileName, mode = "wb" )
  system("cp dssp-2.0.4-win32.exe dssp.exe")
}

dsspExePath <- "C:/research/Summer_Research_2015/code/Using_bio3d/dssp.exe"
#dsspExePath <- "C:/research/Summer_Research_2015/code/Using_bio3d/dssp-2.0.4-win32.exe"

dssp_1AG2 <- bio3d::dssp(pdb_1AG2, exepath = dsspExePath )
dssp_1AG2
dssp_1AG2$acc


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

#mutate_feature("E200K")

# Get the sequence and mutant strings

native_seq_string <- function(PDB_ID) {
  return( paste0( as.vector(pdbseq(read.pdb(PDB_ID))), collapse="") )
}
native_seq_string("1AG2")

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

mutant_seq_string("1AG2","E200K")

# Get solvent accessibility

solventAcc <- function(PDB_ID, mutation) {
  
  pdb <- read.pdb(PDB_ID)
  seq <- pdbseq(pdb)
  acc <- bio3d::dssp(pdb, exepath = dsspExePath )$acc
  mut <- parse_mutation(mutation)
  mut_location <- as.integer(mut[2])
  
  return( acc[ which( names(seq) == mut_location ) ] )
  
}

#solventAcc("1AG2","E200K") # should return 154
#solventAcc("1AG2","F198S") # should return 10

# Using the Peptides package

#**aacomp** The output is a matrix with the number and percentage of amino acids of a particular class

#- Tiny (A + C + G + S + T)
#- Small (A + B + C + D + G + N + P + S + T + V)
#- Aliphatic (A + I + L + V)
#- Aromatic (F + H + W + Y)
#- Non-polar (A + C + F + G + I + L + M + P + V + W + Y)
#- Polar (D + E + H + K + N + Q + R + S + T + Z)
#- Charged (B + D + E + H + K + R + Z)
#- Basic (H + K + R)
#- Acidic (B + D + E + Z)

#library(Peptides)

#native_seq_string_1AG2 <- native_seq_string("1AG2")
#mutant_seq_string_1AG2 <- mutant_seq_string("1AG2","E200K")

#aacomp_native_1AG2 <- aacomp(native_seq_string_1AG2)
#aacomp_mutant_1AG2 <- aacomp(mutant_seq_string_1AG2)

#as.vector(aacomp_mutant_1AG2[,1] - aacomp_native_1AG2[,1])
#as.vector(aacomp_mutant_1AG2[,2] - aacomp_native_1AG2[,2])

#**instaindex** This function calculates the instability index proposed by Guruprasad (1990). A protein whose instability index is smaller than 40 is predicted as stable, a value above 40 predicts that the protein may be unstable.

#instaindex(native_seq_string_1AG2)
#instaindex(mutant_seq_string_1AG2)

#**mw** This function calculates the molecular weight of a protein sequence. It is calculated as the sum of the mass of each amino acid using the scale available on Compute pI/Mw tool.

#mw(native_seq_string_1AG2)
#mw(mutant_seq_string_1AG2)


