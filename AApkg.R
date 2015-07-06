#Amino Acids                                         NH2
#                                                   |
#All Amino Acids have the general composition of   HC--C(O)OH with the 'R' group
#                                                   |         specifying the 
#                                                   R         amino acid residue

#Defining Aromatic 'R' Groups
Phenyl <- c("C",rep("CH",5))
Imidazole <- c("C","CH","NH","CH","N+H")
Indole <- c("C","CH","NH","C","C",rep("CH",4))
Phenol <- c("C","CH","CH","C_OH","CH","CH")
#The following will relate the 1-letter abreviation of 
#the 20 AA residues with their respecitve 'R' groups

A <- Ala <- c("CH3")                              #Alanine (Ala)
C <- Cys <- c("CH2", "SH")                        #Cysteine (Cys)
D <- Asp <- c("CH2","CO","O-")                       #Aspartate (Asp)
E <- Glu <- c("CH2","CH2","CO","O-")                 #Glutamate (Glu)
F <- Phe <- c("CH2",Phenyl)                       #Phenylalanine (Phe)
G <- Gly <- c("H")                                #Glycine (Gly)
H <- His <- c("CH2",Imidazole)                    #Histidine (His)
I <- Ile <- c("Methyl","CH","CH2","CH3")            #IsoLeucine (Ile)
K <- Lys <- c(rep("CH2",4), "N+H3")               #Lysine (Lys)
L <- Leu <- c("CH2","Methyl","CH","CH3")            #Leucine (Leu)
M <- Met <- c("CH2","CH2","S","CH3")              #Methionine (Met)
N <- Asn <- c("CH2","CO","NH2")                   #Asparagine (Asn)
P <- Pro <- c("CH2","CH2","CH2")                  #Proline (Pro)
Q <- Gln <- c("CH2","CH2","CO","NH2")             #Glutamine (Gln)
R <- Arg <- c(rep("CH2",3), "NH", "CN+H2","NH2")  #Arginine (Arg)
S <- Ser <- c("CH2","OH")                         #Serine (Ser)
T <- Thr <- c("CH","OH","CH3")                   #Threonine (Thr)
V <- Val <- c("Methyl","CH","CH3")                  #Valine (Val)
W <- Trp <- c("CH2",Indole)                       #Tryptophan (Trp)
Y <- Tyr <- c("CH2", Phenol)                      #Tyrosine (Tyr)

#Convert from subgroups to individual atoms
Mol2Atm = function(AA){
  AAcid <- gsub("_","",AA); AAcid <- gsub("Methyl","CH3",AAcid)
  AAcid <- gsub("2","H",AAcid); AAcid <- gsub("3","HH",AAcid)
  AAcid <- gsub("N+","X",AAcid,fixed=TRUE)
  AAcid <- gsub("O-","Z",AAcid,fixed=TRUE)
  AAcid <- strsplit(AAcid, split="")
  AAcid <- unlist(AAcid)
  AAcid <- gsub("X","N+",AAcid); AAcid <- gsub("Z","O-",AAcid)
  return(AAcid)
}

#Define path for conecting bonds in the Amino Acid
fpath <- function(AA){
  
  if(length(P) == length(AA) && all(P == AA) == TRUE){
    
    pAA <- c(c(1,2,2,3,3,2,2,4,8),c(1,1,5,5,6),
             c(1,9,9,12,12,15,15,5,5,7),c(10,9,9,11),
             c(13,12,12,14),c(17,15,15,16))
  }
  
  else{
    if(length(A) == length(AA) && all(A == AA) == TRUE){
      pAA <- c(c(1,10,10,11),c(12,10,10,13))
    }
    
    if(length(C) == length(AA) && all(C == AA) == TRUE){
      pAA <- c(c(1,10,10,13,13,14),c(11,10,10,12))
    }
    
    if(length(D) == length(AA) && all(D == AA) == TRUE){
      pAA <- c(c(1,10,10,13,13,14,14,13,13,15),
               c(11,10,10,12))
    }
    
    if(length(E) == length(AA) && all(E == AA) == TRUE){
      pAA <- c(c(1,10,10,13,13,16,16,17,17,16,16,18),
               c(11,10,10,12),c(14,13,13,15))
    }
    
    if(length(F) == length(AA) && all(F == AA) == TRUE){
      pAA <- c(c(1,10,10,13,13,14,14,16,16,18,18,20,20,22,22,13,13,22,22,23),
               c(11,10,10,12),c(15,14,14,16,16,17),c(19,18,18,20,20,21))
    }
    
    if(length(G) == length(AA) && all(G == AA) == TRUE){
      pAA <- c(1,10)
    }
    
    if(length(H) == length(AA) && all(H == AA) == TRUE){
      pAA <- c(c(1,10,10,13,13,14,14,16,16,18,18,20,20,13,13,14,14,15),
               c(11,10,10,12),c(16,17),c(19,18,18,20,20,21))
    }
    
    if(length(I) == length(AA) && all(I == AA) == TRUE){
      pAA <- c(c(1,14,14,15),c(11,10,10,14,14,16,16,19,19,20),
               c(12,10,10,13),c(17,16,16,18),c(21,19,19,22))
      
    }
    if(length(K) == length(AA) && all(K == AA) == TRUE){
      pAA <- c(c(1,10,10,13,13,16,16,19,19,22,22,25),c(11,10,10,12), 
               c(14,13,13,15),c(20,19,19,21),c(23,22,22,24),c(17,16,16,18))
    }
    
    if(length(L) == length(AA) && all(L == AA) == TRUE){
      pAA <- c(c(1,10,10,17,17,13,13,14),c(11,10,10,12),c(15,13,13,16),
               c(18,17,17,19,19,21),c(20,19,19,22))
    }
    
    if(length(M) == length(AA) && all(M == AA) == TRUE){
      pAA <- c(c(1,10,10,13,13,16,16,17,17,18),c(11,10,10,12), 
               c(14,13,13,15),c(20,17,17,19))
    }
    
    if(length(N) == length(AA) && all(N == AA) == TRUE){
      pAA <- c(c(1,10,10,13,13,14,14,13,13,15,15,16),
               c(11,10,10,12),c(17,15))
    }  
    
    if(length(Q) == length(AA) && all(Q == AA) == TRUE){
      pAA <- c(c(1,10,10,13,13,16,16,17,17,16,16,18,18,20),
               c(11,10,10,12),c(14,13,13,15),c(19,18))
    }
    
    if(length(R) == length(AA) && all(R == AA) == TRUE){
      pAA <- c(c(1,10,10,13,13,16,16,19,19,20),c(23,22,22,21,21,19), 
               c(24,22,22,21,21,25),c(26,25,25,27),c(11,10,10,12),
               c(14,13,13,15),c(17,16,16,18))
    }
    
    if(length(S) == length(AA) && all(S == AA) == TRUE){
      pAA <- c(c(1,10,10,13,13,14),c(11,10,10,12))
    }
    
    if(length(T) == length(AA) && all(T == AA) == TRUE){
      pAA <- c(c(1,10,10,14,14,15),c(11,10,10,12,12,13),
               c(16,14,14,17))
    }
    
    if(length(V) == length(AA) && all(V == AA) == TRUE){
      pAA <- c(c(1,14,14,10,10,13),c(11,10,10,12),
               c(15,14,14,16,16,19),c(17,16,16,18))
    }
    
    if(length(W) == length(AA) && all(W == AA) == TRUE){
      pAA <- c(c(1,10,10,13,13,14,14,13,13,19,19,20),c(20,22,22,24),
               c(24,26,26,18,18,19),c(11,10,10,12),c(19,18,18,16,16,17),
               c(15,14,14,16),c(21,20,20,22,22,23),c(25,24,24,26,26,27))
    }
    
    if(length(Y) == length(AA) && all(Y == AA) == TRUE){
      pAA <- c(c(1,10,10,13,13,14,14,16,16,18,18,21),c(21,23,23,13),
               c(13,14,14,15),c(17,16,16,18,18,19,19,20),
               c(22,21,21,23,23,24),c(11,10,10,12))
    }
    AABbPath <- c(1,2,2,3,3,2,2,4,9,1,1,5,5,6,7,5,5,8)
    AAPath <- c(AABbPath,pAA)
    return(AAPath)
  }
}

#Define Vertices and bonds for Amino Acid Backbone (all but the 'R' sidechain)

AABackbone <- c("C","C","O","O-","N+",rep("H",4))

VizAA = function(AA){                             ##visualize the amino acid
  
  if(length(P) == length(AA) && all(P %in% AA) == TRUE  
     && all(AA %in% P) == TRUE){
    
    gAA <- graph.empty(directed=FALSE) + vertices(AABackbone[1:8],Mol2Atm(AA))
    gAA <- gAA + edges(fpath(P))
  }
  else{
    gAA <- graph.empty(directed=FALSE) + vertices(AABackbone,Mol2Atm(AA))
    gAA <- gAA + edges(fpath(AA))
  }
}

#Generate a graphical representation of the Amino Acid

#gAla <- graph.empty(directed=FALSE) + vertices("C","C",rep("H",3),"N+",rep("H",3),"C","O","O-","H")
#gAla <- gAla + path(5,2,4) + path(3,2,1,6,9) + path(7,6,8) + path(13,1,10,11) + path(11,10,12)
#layout <- matrix(c(0,0,0,1.5215,0,0,1.8875,1.0242,0,1.88710,-.51530,-0.88510,1.884,-.513,.8904,-.4972,.7077,1.2161,-.15040,1.6829,1.2197,-1.5319,.7112,1.2183,-.15510,.2212,2.0645,-.5073,-1.4149,-.002,-1.2091,-1.8122,-.9158,-.2057,-2.1995,.95010,-.35860,.51470,-.8867),ncol=3,byrow=TRUE)
#gAla$layout <- layout.norm(layout)
#rglplot(gAla)

AAclassifier <- function(AA){
  
  n = length(AA)
  
  AAclass= data.frame(row.names = "classes")
  
  if(n <= 3){                   #classifier for Small
    #and Tiny AA residues
    AAclass$size <- "Small"     
    
    if(n < 3){
      
      AAclass$size <-"Tiny"
    }
  }
  else{
    AAclass$size <- "Mid"
  }
  
  if(length(grep("+",AA,fixed=TRUE)) !=0 || length(grep("-",AA,fixed=TRUE)) !=0)  {      #classifier for charge
    
    
    if(!is.na(match("O-",AA)) == TRUE){ #classifier for neg. charge
      
      AAclass$charge <- "Negative"
    }
    
    else{
      AAclass$charge <- "Positive"
    }
  }
  
  else{
    AAclass$charge <- "Neutral"
  }
  
  if(length(grep("OH|SH|NH2",AA)) != 0 || length(grep("+",AA,fixed=TRUE)) !=0 ||
     length(grep("-",AA,fixed=TRUE)) !=0 || length(grep("C",AA)) == 0){                   #classifier for hydrophilicity
    
    AAclass$hydroph <- "Hydrophilic"
    
    if(length(grep("SH",AA)) != 0 || length(grep("C",AA)) == 0){        #sub-classifier for polarity
      AAclass$polarity <- "Nonpolar"
    }
    
    else{
      AAclass$polarity <- "Polar"
    }
  }
  
  else{
    AAclass$hydroph <- "Hydrophobic"
    AAclass$polarity <- "Nonpolar"
  }
  
  if(n > 6 && all(AA[2:7] == Phenyl) == TRUE || 
     n > 5 && all(AA[2:6] == Imidazole) == TRUE ||
     n > 9 && all(AA[2:10] == Indole) == TRUE || 
     n > 6 && all(AA[2:7] == Phenol) == TRUE) {      #classifier for Aromaticity
    
    AAclass$Aromatic <- "Yes"
  }
  
  else{
    AAclass$Aromatic <- "No"
  }
  
  if(!is.na(match("Methyl",AA)) == TRUE){    #classifierf for Aliphatic
    
    AAclass$Aliphatic <- "Yes"    
  }
  
  else{
    AAclass$Aliphatic <- "No"
  }
  
  if(length(grep("O|N",AA)) != 0){
    AAclass$HBonds <- "Yes"
  }
  
  else{
    AAclass$HBonds <- "No"
  }
  
  if(length(grep("+",AA,fixed=TRUE)) !=0 || length(grep("-",AA,fixed=TRUE)) !=0)  {
    AAclass$Ion <- "Yes"
  }
  
  else{
    AAclass$Ion <- "No"
  }
  
  if(length(grep("SH",AA)) != 0){
    AAclass$CBond <- "Disulfide"
  }
  
  else{
    AAclass$CBond <- "No"
  }
  
  return(AAclass)
}

#Now need to convert the dataframe classifier into a numerical vector
num.AAclass <- function(AA){
  
  numClass <- as.vector(as.matrix(AAclassifier(AA)))  #convert the the dataframe to a character vector
  
  numClass <- sub("Tiny",0,numClass); numClass <- sub("Small",1,numClass); numClass <- sub("Mid",2,numClass) #numericalize the size
  numClass <- sub("Hydrophilic|Polar", 0,numClass);  numClass <- sub("Hydrophobic|Nonpolar",1,numClass)   #numericalize hydroph & polar
  numClass <- gsub("No",0,numClass);   numClass <- gsub("Yes",1,numClass)  #numericalize the 
  numClass <- sub("Disulfide",1,numClass)                                  #special attributes
  numClass <- sub("Negative",-1,numClass); numClass <- sub("Neutral",0,numClass); numClass <- sub("Positive",1,numClass) #numericalize the charge
  numClass <- as.numeric(numClass)
  
  return(numClass)
}