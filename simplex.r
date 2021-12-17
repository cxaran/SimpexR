library("rjson")

#leemos la instancia
inst <- fromJSON(file = "instance.json")
inst

#leemos la instancia
#t <- read.table("instancia.txt",sep = ",");

#guardamos el numero de variables
#vars <- length(t[1,])-1;vars
vars <- inst$vars; vars

#leer restricciones
t<-matrix(t(unlist(inst$rest)),ncol=length(inst$rest),nrow=vars+2)
sig<-t[vars+1,];sig[sig=="<="] <- 1;sig[sig==">="] <- -1;
t<-matrix(as.numeric(t(t[-(vars+1),])), ncol = ncol(t(t[-(vars+1),])));t

t<-t*as.numeric(sig);t


#leer funcion objetivo
obj<-append(unlist(inst$obj),0)*-1;
#if(inst$type =="max")obj<-obj*-1
t<- rbind(obj,t);t

#agregamos variables de exceso


##agregamos variables de olgura

#colnames(t)[vars+1] <- "R"
#colnames(t)[vars+2] <- "Z"
rownames(t)[1] <- "Z"
t

#guardar
res<-replicate(length(t[1,]), 0)
res[vars+1]<-1
res[vars+1]<-1

while(length(which(t[1,]<0))){
  #encontramos pivote
  pa <- which(t[1,]==min(t[1,])); #t[1,pa]
  if(inst$type =="min"){
    pa <- t[1,]
    pa <- which(t[1,]==max(pa[pa!=0])); #t[1,pa]
  }
  #pb <- t[2:length(t[,1]),pa]/(t[1,pa]); pb
  pb <-t[2:length(t[,1]),pa];
  pc <- t[2:length(t[,1]),vars+1];
  pb[pb>0] <-pc[pb>0]/pb[pb>0];
  if(inst$type =="max"){
    pb[pb<=0] <- Inf;
  }
  pb <- which(pb==min(pb))+1; #t[pb,pa]
  
  print(t[pb,pa])
  
  #guardar
  res[which(res==pb)]=0
  res[pa]<-pb
  
  #guardamos renglon
  tp <- as.numeric(t[pb,]/t[pb,pa]); #tp
  #renglon pivote por fila pivote
  tx <- sweep(t+t[,pa]-t,MARGIN=2,tp,`*`); #tx
  #restamos a tabla original
  t  <- t-tx; #t
  t[pb,]<- tp; #t
  print(t)
}
t
res[which(res>0)]<-t[res,vars+1];
res<-t(as.matrix(res))
colnames(res)<-colnames(t)
res


