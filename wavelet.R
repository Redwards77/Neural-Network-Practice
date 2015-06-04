Haarwavelet = function(firstRow, endRow){
  n=length(Row_1)
  
  fRow_new = function(Row_prev, Row_next, Row_no){
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
  }

