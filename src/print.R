# 3) function for printing a matrix as a table in laTex
mat.print=function(mat,rnames=rownames(mat),cnames=colnames(mat),d1=3,d2=1){
  
  nr=nrow(mat)
  nc=ncol(mat)
  
  # the start of the table
  if(is.null(rnames)){
    tb.start=paste("\\begin{tabular}{*{",nc,"}{l}}",sep="")
  }else{
    tb.start=paste("\\begin{tabular}{l|*{",nc,"}{l}}",sep="")
  }
  
  # set header
  header=""
  if(is.null(cnames)==F){
    if(is.null(rnames)==F){
      header=" &"
    }
    for (k in 1:(nc-1)){
      header=paste(header,cnames[k],"&",sep="")
    }
    header=paste(header,cnames[nc],"\\\\ \\hline ",sep="")  
  }
  
  # string to print (in progress)
  string=paste(tb.start,header,sep="")
  
  for (l in 1:nr){
    # the name of the row
    if(is.null(rnames)){
      row=""  
    }else {
      row=paste(rnames[l],"&")
    }   
    for (k in 1:(nc-1)){
      row=paste(row,my.print(mat[l,k],showDollar=T,d1=d1,d2=d2),"&",sep="")
    }
    row=paste(row,my.print(mat[l,nc],showDollar=T,d1=d1,d2=d2),sep="")
    # include row to string to print (in progress)
    string=paste(string,row,"\\\\",sep="")
  }
  # end of table
  tb.end=paste("\\end{tabular}",sep="")
  
  #paste the final string
  paste(string,tb.end,sep="")
}

my.print=function(x,d1=3,d2=1,showDollar=F){
  xx<-0
  if(x!=0){
    xx<-round(x,d1)
    if(xx==0){
      xx<-scinot(x,digits=d2,showDollar=F,showZero=T)
    }
    else{
      string=paste("%.",d1,"f",sep="")
      xx<-sprintf(string,xx)
    }
  }
  if (showDollar)paste("$",xx,"$",sep="")
  else paste(xx)
}