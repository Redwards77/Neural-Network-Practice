Haarwavelet = function(firstRow, endRow){
  
  k = endRow
  
  fRow_new = function(Row_prev,Row_next,Row_no){
    m = Row_no
    Row_next = c()
    n = length(Row_prev)
      for(i in 1:n){
         if(i <= 2^(m-2)){
            Row_next[i] = Row_prev[i]/2
            #print(Row_next)
         }
         else{
            Row_next[i] = (Row_prev[i-(2^(m-2))]+Row_prev[i])/2
            #print(Row_next)
         }
      }
    return(Row_next)
   }

  fRow_dif = function(Row_prev, Row_next, Row_no){
    m = Row_no
    Row_next = c()
    n = length(Row_prev)
    for(i in 1:n){
      if(i <= 2^(m-1)){
        Row_next[i] = (-Row_prev[i])/2
        #print(Row_next)
      }
      else{
        Row_next[i] = (Row_prev[i-2^(m-1)]-Row_prev[i])/2
        #print(Row_next)
      }
    }
    return(Row_next)
  }
  i = 1:8
  j=c(0:3)
  n= mapply(function(j){i*2^j},j)
  
  
  Averages = matrix(c(-12, 0, 12, 16, 12, 0, -12, -16),1,8)
  Row <- c()
  fAvg = function(Averages){
    (Averages[nrow(Averages),n[i-1,1]]+Averages[nrow(Averages),n[(i),1]])/2
  }
  Rows <- function(Averages){for(i in 1:8){
    if(i == 1){
      Row[i] <- (Averages[nrow(Averages),n[i,1]])/2
    }
    else{
      Row[i] <- (Averages[nrow(Averages),n[i-1,1]]+Averages[nrow(Averages),n[(i),1]])/2
  }
  }
    rbind(Averages,Row)}
}

