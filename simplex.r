library("rjson")

#leemos la instancia
inst <- fromJSON(file = "instance.json")

#guardamos el numero de variables
vars <- inst$vars;

#leer restricciones
t<-matrix(t(unlist(inst$rest)),ncol=length(inst$rest),nrow=vars+2);
sig<-t[vars+1,];
t<-matrix(as.numeric(t(t[-(vars+1),])), ncol = ncol(t(t[-(vars+1),])));
t<-cbind(t[,vars+1],t[,-(vars+1)]);


restricciones<-t[,-1]
derecho<-t[,1]

#agregamos variables >= ==
va<- sig
va[va==">="]=1;
va[va=="=="]=1;
va[va=="<="]=0;

#agregamos variables <=
vo<- sig
vo[vo==">="]=-1;
vo[vo=="=="]=0;
vo[vo=="<="]=1;

#juantamos las variables
vall<-cbind(diag(vo),diag(va));
vall<-vall[, colSums(vall != 0) > 0];

t<-cbind(t,vall);


#leer funcion objetivo
objetivo<-unlist(inst$obj)
t<-rbind(append(as.numeric(append(0,unlist(inst$obj))*-1),replicate(length(vall[1,]),0)),t);

rbind(t[-1,],t[1,])

#guardar
res<-replicate(length(t[1,]), 0)
res[1]<-1

#condicon de parada
f<-FALSE;
fac<-TRUE

#Metodo de dos fases
if(length(va[va==1])>0  ){
  tx<-t;
  tx[1,]<-append(replicate(length(t[1,])-length(va[va==1]),0),replicate(length(va[va==1]),-1));
  cs<- which(va==1);
  ts<- tx[cs+1,];
  ts<-rbind(ts,tx[1,])
  tx[1,]<-colSums(ts);
  tx[1,]<-tx[1,]*-1;
  res<-replicate(length(tx[1,]), 0)
  res[1]<-1
  
  print(rbind(tx[-1,],tx[1,]))
  
  #condicion de parada max
  f<-length(which(tx[1,-1]<0)) > 0
  while(f){
    #obtener pivote
    pa<-tx[1,-1]
    pa <- (which(pa==min(pa)[1])[1])+1;
    pb <-tx[2:length(tx[,1]),pa];
    pc <- tx[2:length(tx[,1]),1];
    pb[pb<=0]<-0;
    pc[pb==0]<-1;
    pb <-pc/pb;
    pb <- (which(pb==min(pb)[1])[1])+1;
    
    print(tx[pb,pa])
    
    #guardar
    res[which(res==pb)]=0;
    res[pa]<-pb;
    
    #Algortimo
    tp <- as.numeric(tx[pb,]/tx[pb,pa]);
    ty <- sweep(tx+tx[,pa]-tx,MARGIN=2,tp,`*`);
    tx <- tx-ty; 
    tx[pb,]<- tp; 
    print(rbind(tx[-1,],tx[1,]))
    
    f<-length(which(tx[1,-1]<0)) > 0
  }
  #fase dos
  to<-t[1,]
  t<-tx[, 1 : (length(tx[1,])-length(va[va==1])) ]
  t[1,]<-to[1:length(t[1,])];
  tm<-(-t[1,]);tm
  tm[res==0]<-0
  tm[tm<0]<-0
  tm<-which(tm>0)
  for (i in tm) {
    min<-which(t[-1,i]>0)[1]
    t[1,]<-t[1,]+(t[min+1,]*((-t[1,i])/t[min+1,i]))
    t[1,]
    tm<-(-t[1,1:(vars+1)]);
    tm[1]<-0
    tm[tm<0]<-0
    tm<-which(tm>0)
  }
  #if(t[1,1]>0) fac<-FALSE
}
rbind(t[-1,],t[1,])

if(fac & inst$type =="max")f<-length(which(t[1,-1]<0)) > 0
if(fac & inst$type =="min")f<-length(which(t[1,-1]>0)) > 0

#Metodo simplex
while(f){
  #encontramos pivote
  pa<-t[1,]
  if(inst$type =="max") pa <- which(pa==min(pa))[1];
  if(inst$type =="min") pa <- which(pa==max(pa))[1];
  pb <-t[2:length(tx[,1]),pa];
  pc <- t[2:length(tx[,1]),1];
  pb[pb<=0]<-0;
  pc[pb<=0]<-1;
  pb <-pc/pb;pb
  pb <- (which(pb==min(pb)[1])[1])+1;
  
  print(t[pb,pa])
  
  #guardar
  res[which(res==pb)]=0
  res[pa]<-pb
  
  #guardamos renglon
  tp <- as.numeric(t[pb,]/t[pb,pa]);
  tx <- sweep(t+t[,pa]-t,MARGIN=2,tp,`*`);
  t  <- t-tx;
  t[pb,]<- tp;
  print(rbind(t[-1,],t[1,]))
  
  if(inst$type =="max")f<-length(which(t[1,-1]<0)) > 0
  if(inst$type =="min")f<-length(which(t[1,-1]>0)) > 0
}
res[which(res>0)]<-t[res,1];
if(!fac) res<-0
res[2:(vars+1)]
res[1]

library(lpSolve)
o <- lp(inst$type,objetivo,restricciones,sig,derecho,all.int = T)
o$solution
o$objval