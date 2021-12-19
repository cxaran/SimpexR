library("rjson")

#leemos la instancia
inst <- fromJSON(file = "instance.json")

#guardamos el numero de variables
vars <- inst$vars; vars

#leer restricciones
t<-matrix(t(unlist(inst$rest)),ncol=length(inst$rest),nrow=vars+2);t
sig<-t[vars+1,];
t<-matrix(as.numeric(t(t[-(vars+1),])), ncol = ncol(t(t[-(vars+1),])));
t<-cbind(t[,vars+1],t[,-(vars+1)]);

#sing<-t[vars+1,];sing[sig=="<="] <- 1;sing[sig==">="] <- -1;
#t<-t*as.numeric(sing);t;

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
#obj<-append(0,unlist(inst$obj))*-1;
#obj<-append(as.numeric(obj),replicate(length(vall[1,]),0)) 
#t<- rbind(obj,t);t
t<-rbind(append(as.numeric(append(0,unlist(inst$obj))*-1),replicate(length(vall[1,]),0)),t);
rbind(t[-1,],t[1,])

#guardar
res<-replicate(length(t[1,]), 0)
res[1]<-1

#condicon de parada
f<-FALSE;
fac<-TRUE

library(lpSolve)
o <- lp(inst$type,objetivo,restricciones,sig,derecho,all.int = T)

#Metodo de dos fases
if(length(va[va==1])>0  ){
  tx<-t;
  tx[1,]<-append(replicate(length(t[1,])-length(va[va==1]),0),replicate(length(va[va==1]),-1));
  #tx<- cbind(tx, append(1,replicate( length(t[,1])-1,0 ))); tx
  
  cs<- which(va==1);
  ts<- tx[cs+1,];
  ts<-rbind(ts,tx[1,])
  tx[1,]<-colSums(ts);
  tx[1,]<-tx[1,]*-1;
  
  #rbind(tx[-1,],tx[1,])
  tx
  
  res<-replicate(length(tx[1,]), 0)
  res[1]<-1
  
  #condicion de parada max
  f<-length(which(tx[1,-1]<0)) > 0
  while(f){
    #obtener pivote
    pa<-tx[1,-1]
    pa <- (which(pa==min(pa)[1])[1])+1; tx[1,pa]; pa
    pb <-tx[2:length(tx[,1]),pa];
    pc <- tx[2:length(tx[,1]),1];
    pb[pb<=0]<-0;pb
    pc[pb==0]<-1;pc
    pb <-pc/pb;pb
    pb <- which(pb==min(pb)[1])+1; tx[pb,pa];pb
    
    #guardar
    res[which(res==pb)]=0
    res[pa]<-pb
    
    #Algortimo
    tp <- as.numeric(tx[pb,]/tx[pb,pa]);tp
    ty <- sweep(tx+tx[,pa]-tx,MARGIN=2,tp,`*`);
    tx <- tx-ty; 
    tx[pb,]<- tp; 
    rbind(tx[-1,],tx[1,])
    
    f<-length(which(tx[1,-1]<0)) > 0
  }
  #tx
  #no se si meter la z
  to<-t[1,]
  t<-tx[, 1 : (length(tx[1,])-length(va[va==1])) ]
  t[1,]<-to[1:length(t[1,])];
  #t<-cbind(t,append(0,matrix(0, nrow =  length(t[,1])))); 
  #t[1,length(t[1,])]<-1;
  t
  
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
f

#Metodo simplex
while(f){
  #encontramos pivote
  pa<-t[1,]
  if(inst$type =="max") pa <- which(pa==min(pa))[1];
  if(inst$type =="min") pa <- which(pa==max(pa))[1];
  t[1,pa]
  pb <-t[2:length(tx[,1]),pa];pb
  pc <- t[2:length(tx[,1]),1];pc
  pb[pb<=0]<-0;
  pc[pb<=0]<-1;
  pb <-pc/pb;pb
  pb <- (which(pb==min(pb))[1])+1; t[pb,pa]
  
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
  print(rbind(t[-1,],t[1,]))
  
  if(inst$type =="max")f<-length(which(t[1,-1]<0)) > 0
  if(inst$type =="min")f<-length(which(t[1,-1]>0)) > 0
  f
}
t
res[which(res>0)]<-t[res,1];
if(!fac) res<-0
res[2:(vars+1)]
res[1]

o$solution
o$objval