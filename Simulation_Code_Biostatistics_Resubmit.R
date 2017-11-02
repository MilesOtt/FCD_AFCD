library(amen)
set.seed(1)

#------------#

q.n<-function(probs){
	p.o<-probs[1:(length(probs)-1)]
	m<-probs[length(probs)]
	n<-length(p.o)
	q<-matrix(rep(0,((m+1)*n)), nrow=n)

p<-p.o
#n.rows by m+1 columns
q[1,1]<-1*(1-p[1])
q[1,2]<-1*p[1]
	for(i in 2:n){
		q[i,1]<-q[(i-1),1]*(1-p[i])
		for(j in 2:min((i+1),(m+1))){
			q[i,j]<-q[(i-1),j]*(1-p[i])+q[(i-1),(j-1)]*p[i]
		}
	}
real.q<-q[n,(m+1)]	
return(real.q)
}

#--------------#


q.n.no.N.i<-function(probs){
  p.o<-probs[1:(length(probs)-1)]
  m<-probs[length(probs)]
  n<-length(p.o)
  
  p<-p.o
  if(m>0){
    q<-matrix(rep(0,((m+1)*n)), nrow=n)
    
    
    #n.rows by m+1 columns
    q[1,1]<-1*(1-p[1])
    q[1,2]<-1*p[1]
    for(i in 2:n){
      q[i,1]<-q[(i-1),1]*(1-p[i])
      for(j in 2:min((i+1),(m+1))){
        q[i,j]<-q[(i-1),j]*(1-p[i])+q[(i-1),(j-1)]*p[i]
      }
    }
    real.q<-q[n,]	
  }
  if (m==0){
    real.q<-prod((1-p))	
  }
  return(real.q)
}
#--------------#

loglike.probit<-function(theta){
  
  #for observed data
  obs<-(miss.vec==0)  	
  y.vec<-y.vecs.full
  theta0<-theta[1]
  theta1<-theta[2]
  theta2<-theta[3]
  
  
  part.1<-sum(y.vec[obs]*log(pnorm((theta0+theta1*x.vecs[obs,2] +theta2*x.vecs[obs,3]))))
  part.2<-sum((1-y.vec[obs])*log(1-pnorm((theta0+theta1*x.vecs[obs,2] +theta2*x.vecs[obs,3]))))
  
  part.3<-0
  #for unobserveddata
  index.no.obs<-index.vec[miss.vec==1]
  prob.no.obs<-pnorm(x.vecs[!obs,]%*%(theta))
  
  missing.row<-(N.i>max.here)
  for(q in index[missing.row]){
    n.missing<-N.i[q]-max.here
    prob.no.obs[index.no.obs==q]<-prob.no.obs[index.no.obs==q]			
    passer<-c(prob.no.obs[index.no.obs==q],(N.i[q]-max.here))
    part.3<-part.3+log(q.n(passer)*(1/choose((max.here+n.missing), max.here)))
  }
  this<-part.1+part.2+part.3
  return(-this)
}

loglike.no.N.i.probit<-function(theta){
  #global variables
  #max.here
  #miss.vec 
  obs<-(miss.vec==0)  	
  y.vec<-y.vecs.full
  theta0<-theta[1]
  theta1<-theta[2]
  theta2<-theta[3]
  
  
  part.1<-sum(y.vec[obs]*log(pnorm((theta0+theta1*x.vecs[obs,2] +theta2*x.vecs[obs,3]))))
  part.2<-sum((1-y.vec[obs])*log(1-pnorm((theta0+theta1*x.vecs[obs,2] +theta2*x.vecs[obs,3]))))
  
  part.3<-0
  #for unobserveddata
  index.no.obs<-index.vec[miss.vec==1]
  prob.no.obs<-pnorm(x.vecs[!obs,]%*%(theta))
  
  
  #---- retry
  missing.row<-(N.i>=max.here)
  for (q in index[missing.row]){
    n.missing<-0:(n.nodes-max.here-1)
    prob.no.obs[index.no.obs==q]<-prob.no.obs[index.no.obs==q]	
    passer<-c(prob.no.obs[index.no.obs==q],(n.nodes-max.here-1))
    part.3.q<-t(1/choose((max.here+n.missing), max.here))%*%q.n.no.N.i(passer)
    part.3<-part.3+log(part.3.q)
  }	
  #---- end retry
  
  this<-part.1+part.2+part.3
  return(-this)
}

#--------#

loglike<-function(theta){
	#global variables
	#max.here
	#miss.vec 
	part.1<-0
	part.2<-0
	obs<-(miss.vec==0)		
	y.vec<-y.vecs.full
	theta0<-theta[1]
	theta1<-theta[2]
	theta2<-theta[3]
	part.1<-sum(y.vec[obs]*(theta0+theta1*x.vecs[obs,2] +theta2*x.vecs[obs,3]))
	part.2<-(exp(theta0+theta1*x.vecs[obs,2]+theta2*x.vecs[obs,3]))
	part.2<-1+part.2
	part.2<--log(part.2)
	part.2<-sum(part.2)
	part.3<-0
	
	index.no.obs<-index.vec[miss.vec==1]
	exp.v.y.no.obs.logit<-(x.vecs[!obs,]%*%(theta))
	exp.v.y.no.obs<-exp(exp.v.y.no.obs.logit)/(1+exp(exp.v.y.no.obs.logit))

	missing.row<-(N.i>max.here)
	for(q in index[missing.row]){
		n.missing<-N.i[q]-max.here
		exp.v.y.no.obs[index.no.obs==q]<-exp.v.y.no.obs[index.no.obs==q]			
		passer<-c(exp.v.y.no.obs[index.no.obs==q],(N.i[q]-max.here))
		part.3<-part.3+log(q.n(passer))
	}
	this<-part.1+part.2+part.3
	return(-this)
}


#-------------------------------- End Functions-----------------------------------------------#



keep.track.N.i.only.0.0.0<-NULL
keep.track.N.i.only.n1.n1.n1<-NULL
keep.track.N.i.only.n5.0.n5<-NULL

keep.track.N.i.1<-NULL
keep.track.N.i.2<-NULL
keep.track.N.i.3<-NULL
keep.track.N.i.4<-NULL
keep.track.N.i.5<-NULL
keep.track.N.i.6<-NULL
keep.track.N.i.7<-NULL
keep.track.N.i.8<-NULL
keep.track.N.i.9<-NULL
keep.track.full<-NULL

keep.track.straw.1<-NULL
keep.track.straw.2<-NULL
keep.track.straw.3<-NULL
keep.track.straw.4<-NULL
keep.track.straw.5<-NULL
keep.track.straw.6<-NULL
keep.track.straw.7<-NULL
keep.track.straw.8<-NULL
keep.track.straw.9<-NULL


keep.track.CB.1<-NULL
keep.track.CB.2<-NULL
keep.track.CB.3<-NULL
keep.track.CB.4<-NULL
keep.track.CB.5<-NULL
keep.track.CB.6<-NULL
keep.track.CB.7<-NULL
keep.track.CB.8<-NULL
keep.track.CB.9<-NULL


keep.track.no.N.i.1<-NULL
keep.track.no.N.i.2<-NULL
keep.track.no.N.i.3<-NULL
keep.track.no.N.i.4<-NULL
keep.track.no.N.i.5<-NULL
keep.track.no.N.i.6<-NULL
keep.track.no.N.i.7<-NULL
keep.track.no.N.i.8<-NULL
keep.track.no.N.i.9<-NULL


n.nodes<-100
N.i.record<-rep(0, (n.nodes+1))

x.0.matrix.data<-matrix(rep(1, (n.nodes*n.nodes)), nrow=n.nodes)

young.age<-c(seq(18,20), seq(18.5,27.5), seq(18.25,28.25), seq(18.75,34))
middle.age<-c(seq(34,50), seq(34.25, 40.25), seq(38.5, 53.5))
older.age<-c(seq(50,63), seq(60,65))

age<-c(young.age, middle.age, older.age)

index<-1:n.nodes
index.matrix<-matrix(index, ncol=n.nodes,nrow=n.nodes)
x.1.matrix.data<-matrix(age, ncol=n.nodes,nrow=n.nodes)
x.2.matrix.data<-x.1.matrix.data

for(i in 1:n.nodes){
	for(j in 1:n.nodes){
		x.2.matrix.data[i,j]<-abs(age[i]-age[j])
	}
}


#--------- setting beta parameters---------------#
beta.0<-(-1); beta.1<-(.002); beta.2<-(-0.025)



n.loop<-101

#---------------------------------------------#
for(iii in 67:n.loop){




#---- setting up the data--------#

y.matrix.data.full<-matrix(rep(0, (n.nodes*n.nodes)), nrow=n.nodes)
for(i in 1:n.nodes){
	for(j in 1:n.nodes){
		probit<-pnorm(beta.0+beta.1*x.1.matrix.data[i,j] +beta.2*x.2.matrix.data[i,j] )
		y.matrix.data.full[i,j]<-rbinom(1, 1, probit)
	}
}
diag(y.matrix.data.full)<-0
N.i<-rowSums(y.matrix.data.full)
table(N.i)
hist(N.i)
for(nn in 1:n.nodes)
	N.i.record[(N.i[nn]+1)]<-N.i.record[(N.i[nn]+1)]+1
dim.y<-dim(y.matrix.data.full)[1]


y.matrix.data.9<-y.matrix.data.full
y.matrix.data.8<-y.matrix.data.full
y.matrix.data.7<-y.matrix.data.full
y.matrix.data.6<-y.matrix.data.full
y.matrix.data.5<-y.matrix.data.full
y.matrix.data.4<-y.matrix.data.full
y.matrix.data.3<-y.matrix.data.full
y.matrix.data.2<-y.matrix.data.full
y.matrix.data.1<-y.matrix.data.full

y.miss.9<-matrix(rep(0, (n.nodes*n.nodes)), nrow=n.nodes)
y.miss.8<-matrix(rep(0, (n.nodes*n.nodes)), nrow=n.nodes)
y.miss.7<-matrix(rep(0, (n.nodes*n.nodes)), nrow=n.nodes)
y.miss.6<-matrix(rep(0, (n.nodes*n.nodes)), nrow=n.nodes)
y.miss.5<-matrix(rep(0, (n.nodes*n.nodes)), nrow=n.nodes)
y.miss.4<-matrix(rep(0, (n.nodes*n.nodes)), nrow=n.nodes)
y.miss.3<-matrix(rep(0, (n.nodes*n.nodes)), nrow=n.nodes)
y.miss.2<-matrix(rep(0, (n.nodes*n.nodes)), nrow=n.nodes)
y.miss.1<-matrix(rep(0, (n.nodes*n.nodes)), nrow=n.nodes)

y.miss.no.N.i.9<-matrix(rep(0, (n.nodes*n.nodes)), nrow=n.nodes)
y.miss.no.N.i.8<-matrix(rep(0, (n.nodes*n.nodes)), nrow=n.nodes)
y.miss.no.N.i.7<-matrix(rep(0, (n.nodes*n.nodes)), nrow=n.nodes)
y.miss.no.N.i.6<-matrix(rep(0, (n.nodes*n.nodes)), nrow=n.nodes)
y.miss.no.N.i.5<-matrix(rep(0, (n.nodes*n.nodes)), nrow=n.nodes)
y.miss.no.N.i.4<-matrix(rep(0, (n.nodes*n.nodes)), nrow=n.nodes)
y.miss.no.N.i.3<-matrix(rep(0, (n.nodes*n.nodes)), nrow=n.nodes)
y.miss.no.N.i.2<-matrix(rep(0, (n.nodes*n.nodes)), nrow=n.nodes)
y.miss.no.N.i.1<-matrix(rep(0, (n.nodes*n.nodes)), nrow=n.nodes)


index<-1:dim.y


index<-1:dim.y
for(i in 1:dim.y){
  t.c<-index[y.matrix.data.full[i,]==1]
  t.c.o<-sample(t.c, length(t.c))
  n.c<-length(t.c.o)
  
  if(n.c>1){
    y.matrix.data.1[i,t.c.o[2:n.c]]<-0
    y.miss.1[i,y.matrix.data.1[i,]==0]<-1
    y.miss.no.N.i.1[i,y.matrix.data.1[i,]==0]<-1
  }
  if (n.c==1){
    y.miss.no.N.i.1[i,y.matrix.data.1[i,]==0]<-1
  }
  if(n.c>2){
    y.matrix.data.2[i,t.c.o[3:n.c]]<-0
    y.miss.2[i,y.matrix.data.2[i,]==0]<-1
    y.miss.no.N.i.2[i,y.matrix.data.2[i,]==0]<-1
  }
  if (n.c==2){
    y.miss.no.N.i.2[i,y.matrix.data.2[i,]==0]<-1
  }
  
  
  
  if(n.c>4){
    y.matrix.data.4[i,t.c.o[5:n.c]]<-0
    y.miss.4[i,y.matrix.data.4[i,]==0]<-1
    y.miss.no.N.i.4[i,y.matrix.data.4[i,]==0]<-1
  }
  if (n.c==4){
    y.miss.no.N.i.4[i,y.matrix.data.4[i,]==0]<-1
  }
  
 
  
  if(n.c>6){
    y.matrix.data.6[i,t.c.o[7:n.c]]<-0
    y.miss.6[i,y.matrix.data.6[i,]==0]<-1
    y.miss.no.N.i.6[i,y.matrix.data.6[i,]==0]<-1
  }
  if (n.c==6){
    y.miss.no.N.i.6[i,y.matrix.data.6[i,]==0]<-1
  }
  
 
  
  if(n.c>8){
    y.matrix.data.8[i,t.c.o[9:n.c]]<-0
    y.miss.8[i,y.matrix.data.8[i,]==0]<-1
    y.miss.no.N.i.8[i,y.matrix.data.8[i,]==0]<-1
  }
  if (n.c==8){
    y.miss.no.N.i.8[i,y.matrix.data.8[i,]==0]<-1
  }
  
  
  
	}






vectorized.x.0<-c(x.0.matrix.data[upper.tri(x.0.matrix.data)], x.0.matrix.data[lower.tri(x.0.matrix.data)])
vectorized.x.1<-c(x.1.matrix.data[upper.tri(x.1.matrix.data)], x.1.matrix.data[lower.tri(x.1.matrix.data)])
vectorized.x.2<-c(x.2.matrix.data[upper.tri(x.2.matrix.data)], x.2.matrix.data[lower.tri(x.2.matrix.data)])

x.vecs.a<-cbind((t(vectorized.x.0)),(t(vectorized.x.1)), (t(vectorized.x.2)))
x.vecs<-matrix(x.vecs.a, ncol=3, byrow=FALSE)
y.vecs.full<-c(y.matrix.data.full[upper.tri(y.matrix.data.full)], y.matrix.data.full[lower.tri(y.matrix.data.full)])
y.miss.vec.1<-c(y.miss.1[upper.tri(y.miss.5)], y.miss.1[lower.tri(y.miss.5)])
y.miss.vec.2<-c(y.miss.2[upper.tri(y.miss.5)], y.miss.2[lower.tri(y.miss.5)])
y.miss.vec.3<-c(y.miss.3[upper.tri(y.miss.5)], y.miss.3[lower.tri(y.miss.5)])
y.miss.vec.4<-c(y.miss.4[upper.tri(y.miss.5)], y.miss.4[lower.tri(y.miss.5)])
y.miss.vec.5<-c(y.miss.5[upper.tri(y.miss.5)], y.miss.5[lower.tri(y.miss.5)])
y.miss.vec.6<-c(y.miss.6[upper.tri(y.miss.5)], y.miss.6[lower.tri(y.miss.5)])
y.miss.vec.7<-c(y.miss.7[upper.tri(y.miss.5)], y.miss.7[lower.tri(y.miss.5)])
y.miss.vec.8<-c(y.miss.8[upper.tri(y.miss.5)], y.miss.8[lower.tri(y.miss.5)])
y.miss.vec.9<-c(y.miss.9[upper.tri(y.miss.5)], y.miss.9[lower.tri(y.miss.5)])

y.vec.1<-c(y.matrix.data.1[upper.tri(y.miss.5)], y.matrix.data.1[lower.tri(y.miss.5)])
y.vec.2<-c(y.matrix.data.2[upper.tri(y.miss.5)], y.matrix.data.2[lower.tri(y.miss.5)])
y.vec.3<-c(y.matrix.data.3[upper.tri(y.miss.5)], y.matrix.data.3[lower.tri(y.miss.5)])
y.vec.4<-c(y.matrix.data.4[upper.tri(y.miss.5)], y.matrix.data.4[lower.tri(y.miss.5)])
y.vec.5<-c(y.matrix.data.5[upper.tri(y.miss.5)], y.matrix.data.5[lower.tri(y.miss.5)])
y.vec.6<-c(y.matrix.data.6[upper.tri(y.miss.5)], y.matrix.data.6[lower.tri(y.miss.5)])
y.vec.7<-c(y.matrix.data.7[upper.tri(y.miss.5)], y.matrix.data.7[lower.tri(y.miss.5)])
y.vec.8<-c(y.matrix.data.8[upper.tri(y.miss.5)], y.matrix.data.8[lower.tri(y.miss.5)])
y.vec.9<-c(y.matrix.data.9[upper.tri(y.miss.5)], y.matrix.data.9[lower.tri(y.miss.5)])

y.miss.no.N.i.vec.1<-c(y.miss.no.N.i.1[upper.tri(y.miss.5)], y.miss.no.N.i.1[lower.tri(y.miss.5)])
y.miss.no.N.i.vec.2<-c(y.miss.no.N.i.2[upper.tri(y.miss.5)], y.miss.no.N.i.2[lower.tri(y.miss.5)])
y.miss.no.N.i.vec.3<-c(y.miss.no.N.i.3[upper.tri(y.miss.5)], y.miss.no.N.i.3[lower.tri(y.miss.5)])
y.miss.no.N.i.vec.4<-c(y.miss.no.N.i.4[upper.tri(y.miss.5)], y.miss.no.N.i.4[lower.tri(y.miss.5)])
y.miss.no.N.i.vec.5<-c(y.miss.no.N.i.5[upper.tri(y.miss.5)], y.miss.no.N.i.5[lower.tri(y.miss.5)])
y.miss.no.N.i.vec.6<-c(y.miss.no.N.i.6[upper.tri(y.miss.5)], y.miss.no.N.i.6[lower.tri(y.miss.5)])
y.miss.no.N.i.vec.7<-c(y.miss.no.N.i.7[upper.tri(y.miss.5)], y.miss.no.N.i.7[lower.tri(y.miss.5)])
y.miss.no.N.i.vec.8<-c(y.miss.no.N.i.8[upper.tri(y.miss.5)], y.miss.no.N.i.8[lower.tri(y.miss.5)])
y.miss.no.N.i.vec.9<-c(y.miss.no.N.i.9[upper.tri(y.miss.5)], y.miss.no.N.i.9[lower.tri(y.miss.5)])



index.vec<-c(index.matrix[upper.tri(index.matrix)], index.matrix[lower.tri(index.matrix)])

#--------------- End Creating the Data--------------#


#------------N.i only -----------------------------#

max.here<-0
miss.vec<-y.miss.vec.5
miss.vec[miss.vec==0]<-1
theta.try.1<-c(0,0,0)
params.N.i.only.0.0.0<-optim(theta.try.1, loglike)$par
keep.track.N.i.only.0.0.0<-rbind(keep.track.N.i.only.0.0.0,params.N.i.only.0.0.0)


#------------ N.i 2--------------------------------#
max.here<-2
miss.vec<-y.miss.vec.2
theta.try<-c(0,0,0)
params2<-optim(theta.try, loglike.probit)$par
keep.track.N.i.2<-rbind(keep.track.N.i.2,params2)



#------------ N.i 4--------------------------------#
max.here<-4
miss.vec<-y.miss.vec.4
theta.try<-c(0,0,0)
params<-optim(theta.try, loglike.probit)$par
keep.track.N.i.4<-rbind(keep.track.N.i.4,params)


#------------ N.i 6--------------------------------#
max.here<-6
miss.vec<-y.miss.vec.6
theta.try<-c(0,0,0)
params<-optim(theta.try, loglike.probit)$par
keep.track.N.i.6<-rbind(keep.track.N.i.6,params)


##------------ N.i 8--------------------------------#
max.here<-8
miss.vec<-y.miss.vec.8
theta.try<-c(0,0,0)
params<-optim(theta.try, loglike.probit)$par
keep.track.N.i.8<-rbind(keep.track.N.i.8,params)

#--------- Full Data ------------------------------#

if(iii==1){
fit.full<-glm(y.vecs.full~x.vecs[,2]+x.vecs[,3],family=binomial(link="probit"))
params.full.real<-fit.full$coef
}
keep.track.full<-rbind(keep.track.full,params.full.real)



#--------- Straw 2 ------------------------------#
fit.full<-glm(y.vec.2~x.vecs[,2]+x.vecs[,3], family=binomial(link="probit"))
params.full<-fit.full$coef
keep.track.straw.2<-rbind(keep.track.straw.2,params.full)


#--------- Straw 4 ------------------------------#
fit.full<-glm(y.vec.4~x.vecs[,2]+x.vecs[,3], family=binomial(link="probit"))
params.full<-fit.full$coef
keep.track.straw.4<-rbind(keep.track.straw.4,params.full)



#--------- Straw 6 ------------------------------#
fit.full<-glm(y.vec.6~x.vecs[,2]+x.vecs[,3], family=binomial(link="probit"))
params.full<-fit.full$coef
keep.track.straw.6<-rbind(keep.track.straw.6,params.full)


#--------- Straw 8 ------------------------------#
fit.full<-glm(y.vec.8~x.vecs[,2]+x.vecs[,3], family=binomial(link="probit"))
params.full<-fit.full$coef
keep.track.straw.8<-rbind(keep.track.straw.8,params.full)



#------------ no N.i 2--------------------------------#
max.here<-2
miss.vec<-y.miss.no.N.i.vec.2
theta.try<-c(0,0,0)
params1<-optim(theta.try, loglike.no.N.i.probit, method="BFGS")$par
keep.track.no.N.i.2<-rbind(keep.track.no.N.i.2,params1)

#------------ no N.i 4--------------------------------#
max.here<-4
miss.vec<-y.miss.no.N.i.vec.4
theta.try<-c(0,0,0)
params1<-optim(theta.try, loglike.no.N.i.probit, method="BFGS")$par
keep.track.no.N.i.4<-rbind(keep.track.no.N.i.4,params1)


#------------ no N.i 6--------------------------------#
max.here<-6
miss.vec<-y.miss.no.N.i.vec.6
theta.try<-c(0,0,0)
params1<-optim(theta.try, loglike.no.N.i.probit, method="BFGS")$par
keep.track.no.N.i.6<-rbind(keep.track.no.N.i.6,params1)


#------------ no N.i 8--------------------------------#
max.here<-8
miss.vec<-y.miss.no.N.i.vec.8
theta.try<-c(0,0,0)
params1<-optim(theta.try, loglike.no.N.i.probit, method="BFGS")$par
keep.track.no.N.i.8<-rbind(keep.track.no.N.i.8,params1)


#--------- CB 2 ------------------------------#

fit.full<-ame(y.matrix.data.2, Xdyad=x.2.matrix.data, Xrow=age, 
              model="cbin", intercept=TRUE, odmax=2, print=FALSE, plot=FALSE,
              burn=5000,nscan=50000, odens=25, seed=iii)


params.full<-colMeans(fit.full$BETA)
keep.track.CB.2<-rbind(keep.track.CB.2,params.full)


#--------- CB 4 ------------------------------#
fit.full<-ame(y.matrix.data.4, Xdyad=x.2.matrix.data, Xrow=age, 
              model="cbin", intercept=TRUE, odmax=4, print=FALSE, plot=FALSE,
              burn=5000,nscan=50000, odens=25, seed=iii)


params.full<-colMeans(fit.full$BETA)
keep.track.CB.4<-rbind(keep.track.CB.4,params.full)

#--------- CB 6 ------------------------------#
fit.full<-ame(y.matrix.data.6, Xdyad=x.2.matrix.data, Xrow=age, 
              model="cbin", intercept=TRUE, odmax=6, print=FALSE, plot=FALSE,
              burn=5000,nscan=50000, odens=25, seed=iii)


params.full<-colMeans(fit.full$BETA)
keep.track.CB.6<-rbind(keep.track.CB.6,params.full)


#--------- CB 8 ------------------------------#
fit.full<-ame(y.matrix.data.8, Xdyad=x.2.matrix.data, Xrow=age, 
              model="cbin", intercept=TRUE, odmax=8, print=FALSE, plot=FALSE,
              burn=5000,nscan=50000, odens=25, seed=iii)


params.full<-colMeans(fit.full$BETA)
keep.track.CB.8<-rbind(keep.track.CB.8,params.full)



print(iii)
if(iii %in% c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100))
  save(list = ls(all = TRUE),file=paste("AFCD_CB_Likelihood_Sim_100_loop_probit",iii,"_5_21_2017.RData"))

}


#----------Get rid of index 67 for the ones with 101 observations

#------------- set iii to 100 

keep.track.CB.2.new<-keep.track.CB.2[1:100,]
keep.track.CB.4.new<-keep.track.CB.4[1:100,]
keep.track.CB.6.new<-keep.track.CB.6[1:100,]
keep.track.CB.8.new<-keep.track.CB.8[1:100,]

keep.track.full.new<-keep.track.full[-c(67,102),]

keep.track.N.i.2.new<-keep.track.N.i.2[-c(67,102),]
keep.track.N.i.4.new<-keep.track.N.i.4[-c(67,102),]
keep.track.N.i.6.new<-keep.track.N.i.6[-c(67,102),]
keep.track.N.i.8.new<-keep.track.N.i.8[-c(67,102),]
keep.track.N.i.only.0.0.0.new<-keep.track.N.i.only.0.0.0[-c(67,102),]

keep.track.no.N.i.2.new<-keep.track.no.N.i.2[-c(67,102),]
keep.track.no.N.i.4.new<-keep.track.no.N.i.4[-c(67,102),]
keep.track.no.N.i.6.new<-keep.track.no.N.i.6[-c(67,102),]
keep.track.no.N.i.8.new<-keep.track.no.N.i.8[-c(67,102),]

keep.track.straw.2.new<-keep.track.straw.2[-c(67,102),]
keep.track.straw.4.new<-keep.track.straw.4[-c(67,102),]
keep.track.straw.6.new<-keep.track.straw.6[-c(67,102),]
keep.track.straw.8.new<-keep.track.straw.8[-c(67,102),]



save(list = ls(all = TRUE),file=paste("AFCD_CB_Likelihood_Sim_100_loop_probit_100_fixed_5_21_2017.RData"))

iii<-100


#---------fix stuff below

labs.AFCD<-c( rep("8a",iii), rep("6a",iii), rep("4a",iii), rep("2a",iii), rep("0a",iii))
labs.CC<-c( rep("8c",iii), rep("6c",iii), rep("4c",iii), rep("2c",iii))
labs.FCD<-c(  rep("8b",iii), rep("6b",iii), rep("4b",iii), rep("2b",iii))
labs.CB<-c(  rep("8d",iii), rep("6d",iii), rep("4d",iii), rep("2d",iii))


beta0.AFCD<-c(keep.track.N.i.8.new[1:iii,1],keep.track.N.i.6.new[1:iii,1],keep.track.N.i.4.new[1:iii,1] ,
              keep.track.N.i.2.new[1:iii,1],keep.track.N.i.only.0.0.0.new[1:iii,1] )

beta0.CC<-c(keep.track.straw.8.new[1:iii,1],keep.track.straw.6.new[1:iii,1],keep.track.straw.4.new[1:iii,1], 
            keep.track.straw.2.new[1:iii,1])

beta0.FCD<-c(keep.track.no.N.i.8.new[1:iii,1],keep.track.no.N.i.6.new[1:iii,1],keep.track.no.N.i.4.new[1:iii,1] ,
             keep.track.no.N.i.2.new[1:iii,1] )

beta0.CB<-c(keep.track.CB.8.new[1:iii,1],keep.track.CB.6.new[1:iii,1],keep.track.CB.4.new[1:iii,1] ,
            keep.track.CB.2.new[1:iii,1] )



beta1.AFCD<-c(keep.track.N.i.8.new[1:iii,2],keep.track.N.i.6.new[1:iii,2],keep.track.N.i.4.new[1:iii,2] ,
              keep.track.N.i.2.new[1:iii,2],keep.track.N.i.only.0.0.0.new[1:iii,2] )

beta1.CC<-c(keep.track.straw.8.new[1:iii,2],keep.track.straw.6.new[1:iii,2],keep.track.straw.4.new[1:iii,2], 
            keep.track.straw.2.new[1:iii,2])

beta1.FCD<-c(keep.track.no.N.i.8.new[1:iii,2],keep.track.no.N.i.6.new[1:iii,2],keep.track.no.N.i.4.new[1:iii,2] ,
             keep.track.no.N.i.2.new[1:iii,2] )

beta1.CB<-c(keep.track.CB.8.new[1:iii,2],keep.track.CB.6.new[1:iii,2],keep.track.CB.4.new[1:iii,2] ,
            keep.track.CB.2.new[1:iii,2] )


beta2.AFCD<-c(keep.track.N.i.8.new[1:iii,3],keep.track.N.i.6.new[1:iii,3],keep.track.N.i.4.new[1:iii,3] ,
              keep.track.N.i.2.new[1:iii,3],keep.track.N.i.only.0.0.0.new[1:iii,3] )

beta2.CC<-c(keep.track.straw.8.new[1:iii,3],keep.track.straw.6.new[1:iii,3],keep.track.straw.4.new[1:iii,3], 
            keep.track.straw.2.new[1:iii,3])

beta2.FCD<-c(keep.track.no.N.i.8.new[1:iii,3],keep.track.no.N.i.6.new[1:iii,3],keep.track.no.N.i.4.new[1:iii,3] ,
             keep.track.no.N.i.2.new[1:iii,3] )

beta2.CB<-c(keep.track.CB.8.new[1:iii,3],keep.track.CB.6.new[1:iii,3],keep.track.CB.4.new[1:iii,3] ,
            keep.track.CB.2.new[1:iii,3] )


estimator=c("0", "2","4","6", "8")

pdf("Sim100_Beta0_with_CB_New.pdf")
par(mfrow=c(1,1), mar=c(4,5,1,2))

boxplot(c(beta0.AFCD, beta0.FCD,beta0.CC, beta0.CB)~c(labs.AFCD,labs.FCD,labs.CC, labs.CB) , 
        ylab= expression(beta[0]), xlab="Maximum Number of Nominations", cex.lab=2, cex.axis=2, 
        xlim=c(0.8, 6.1),boxwex=.15,yaxt="n",
        at=c(1, 1.8,2,2.2,2.4, 3,3.2,3.4,3.6, 4.2,4.4,4.6,4.8, 5.4,5.6,5.8,6 ),  xaxt="n")

abline(h=beta.0)


boxplot(c(beta0.AFCD, beta0.FCD,beta0.CC, beta0.CB)~
          c(labs.AFCD,labs.FCD,labs.CC, labs.CB), ylab=expression(beta[0]), 
        xlab="Maximum Number of Nominations", add=TRUE,cex.lab=2, cex.axis=2,
        col=(c("blue","blue", "red", "green","yellow", 
               "blue", "red", "green","yellow",
               "blue", "red", "green","yellow",
               "blue", "red", "green", "yellow")),
        xlim=c(0.8, 6.1),,  xaxt="n",boxwex=.15,at=c(1, 1.8,2,2.2,2.4, 3,3.2,3.4,3.6, 4.2,4.4,4.6,4.8, 5.4,5.6,5.8,6 ))
axis(1, at=(c(1, 2.1, 3.3, 4.5, 5.7)),labels=estimator, cex.axis=2,cex=2, tick=FALSE)
#leg.text<-c("AFCD", "FCD", "Naive", "CB")
#legend("bottomright", legend=leg.text, fill=(c("blue", "red", "green", "yellow")), pt.cex=c(1.4,1,1, 1.4), bty="n")


dev.off()

pdf("Sim100_Beta1_with_CB_new.pdf")

par(mfrow=c(1,1), mar=c(4,5,1,2))

boxplot(c(beta1.AFCD, beta1.FCD,beta1.CC, beta1.CB)~c(labs.AFCD,labs.FCD,labs.CC, labs.CB) , 
        ylab=expression(beta[1]), xlab="Maximum Number of Nominations", cex.lab=2, cex.axis=2, 
        xlim=c(0.8, 6.1),boxwex=.15,yaxt="n", 
        at=c(1, 1.8,2,2.2,2.4, 3,3.2,3.4,3.6, 4.2,4.4,4.6,4.8, 5.4,5.6,5.8,6 ),  xaxt="n")
abline(h=beta.1)

boxplot(c(beta1.AFCD, beta1.FCD,beta1.CC, beta1.CB)~
          c(labs.AFCD,labs.FCD,labs.CC, labs.CB), ylab=expression(beta[1]), cex.lab=2, cex.axis=2,
        xlab="Maximum Number of Nominations", add=TRUE,
        col=(c("blue","blue", "red", "green","yellow", 
               "blue", "red", "green","yellow",
               "blue", "red", "green","yellow",
               "blue", "red", "green", "yellow")),
        xlim=c(0.8, 6.1),,  xaxt="n",boxwex=.15,at=c(1, 1.8,2,2.2,2.4, 3,3.2,3.4,3.6, 4.2,4.4,4.6,4.8, 5.4,5.6,5.8,6 ))
axis(1, at=(c(1, 2.1, 3.3, 4.5, 5.7)),labels=estimator, cex.axis=2,cex=2, tick=FALSE)
#leg.text<-c("AFCD", "FCD", "Naive", "CB")
#legend("bottomright", legend=leg.text, fill=(c("blue", "red", "green", "yellow")), pt.cex=c(1.4,1,1, 1.4), bty="n")

dev.off()


pdf("Sim100_Beta2_with_CB_new.pdf")
par(mfrow=c(1,1), mar=c(4,5,1,2))

boxplot(c(beta2.AFCD, beta2.FCD,beta2.CC, beta2.CB)~c(labs.AFCD,labs.FCD,labs.CC, labs.CB) , 
        ylab=expression(beta[2]), xlab="Maximum Number of Nominations", cex.lab=2, cex.axis=2, 
        xlim=c(0.8, 6.1),boxwex=.15,yaxt="n", ylim=c(-.20, -.01),
        at=c(1, 1.8,2,2.2,2.4, 3,3.2,3.4,3.6, 4.2,4.4,4.6,4.8, 5.4,5.6,5.8,6 ),  xaxt="n")

abline(h=beta.2)

boxplot(c(beta2.AFCD, beta2.FCD,beta2.CC, beta2.CB)~
          c(labs.AFCD,labs.FCD,labs.CC, labs.CB), ylab=expression(beta[2]), 
        xlab="Maximum Number of Nominations", add=TRUE,cex.lab=2, cex.axis=2,ylim=c(-.20, -.01),
        col=(c("blue","blue", "red", "green","yellow", 
               "blue", "red", "green","yellow",
               "blue", "red", "green","yellow",
               "blue", "red", "green", "yellow")),
        xlim=c(0.8, 6.1),,  xaxt="n",boxwex=.15,at=c(1, 1.8,2,2.2,2.4, 3,3.2,3.4,3.6, 4.2,4.4,4.6,4.8, 5.4,5.6,5.8,6 ))
axis(1, at=(c(1, 2.1, 3.3, 4.5, 5.7)),labels=estimator, cex.axis=2, tick=FALSE)
leg.text<-c("AFCD", "FCD", "Naive", "CB")
legend("bottomright", legend=leg.text, fill=(c("blue", "red", "green", "yellow")), 
       pt.cex=c(1.4,1,1, 1.4), bty="n", cex=2)


dev.off()


#beta 0-----MSE ----------#

mean((keep.track.full.new[,1]-beta.0)^2)*10^2
mean((keep.track.N.i.only.0.0.0.new[,1]-beta.0)^2)*10^2

mean((keep.track.N.i.8.new[,1]-beta.0)^2)*10^2
mean((keep.track.no.N.i.8.new[,1]-beta.0)^2)*10^2
mean((keep.track.straw.8.new[,1]-beta.0)^2)*10^2
mean((keep.track.CB.8.new[,1]-beta.0)^2)*10^2

mean((keep.track.N.i.6.new[,1]-beta.0)^2)*10^2
mean((keep.track.no.N.i.6.new[,1]-beta.0)^2)*10^2
mean((keep.track.straw.6.new[,1]-beta.0)^2)*10^2
mean((keep.track.CB.6.new[,1]-beta.0)^2)*10^2

mean((keep.track.N.i.4.new[,1]-beta.0)^2)*10^2
mean((keep.track.no.N.i.4.new[,1]-beta.0)^2)*10^2
mean((keep.track.straw.4.new[,1]-beta.0)^2)*10^2
mean((keep.track.CB.4.new[,1]-beta.0)^2)*10^2

mean((keep.track.N.i.2.new[,1]-beta.0)^2)*10^2
mean((keep.track.no.N.i.2.new[,1]-beta.0)^2)*10^2
mean((keep.track.straw.2.new[,1]-beta.0)^2)*10^2
mean((keep.track.CB.2.new[,1]-beta.0)^2)*10^2

#beta 0-----bias ----------#

mean((keep.track.full.new[,1]-beta.0))*10
mean((keep.track.N.i.only.0.0.0.new[,1]-beta.0))*10

mean((keep.track.N.i.8.new[,1]-beta.0))*10
mean((keep.track.no.N.i.8.new[,1]-beta.0))*10
mean((keep.track.straw.8.new[,1]-beta.0))*10
mean((keep.track.CB.8.new[,1]-beta.0))*10


mean((keep.track.N.i.6.new[,1]-beta.0))*10
mean((keep.track.no.N.i.6.new[,1]-beta.0))*10
mean((keep.track.straw.6.new[,1]-beta.0))*10
mean((keep.track.CB.6.new[,1]-beta.0))*10

mean((keep.track.N.i.4.new[,1]-beta.0))*10
mean((keep.track.no.N.i.4.new[,1]-beta.0))*10
mean((keep.track.straw.4.new[,1]-beta.0))*10
mean((keep.track.CB.4.new[,1]-beta.0))*10

mean((keep.track.N.i.2.new[,1]-beta.0))*10
mean((keep.track.no.N.i.2.new[,1]-beta.0))*10
mean((keep.track.straw.2.new[,1]-beta.0))*10
mean((keep.track.CB.2.new[,1]-beta.0))*10


#beta 0-----sd ----------#
sd((keep.track.full.new[,1]))*10
sd((keep.track.N.i.only.0.0.0.new[,1]))*10

sd((keep.track.N.i.8.new[,1]))*10
sd((keep.track.no.N.i.8.new[,1]))*10
sd((keep.track.straw.8.new[,1]))*10
sd((keep.track.CB.8.new[,1]))*10

sd((keep.track.N.i.6.new[,1]))*10
sd((keep.track.no.N.i.6.new[,1]))*10
sd((keep.track.straw.6.new[,1]))*10
sd((keep.track.CB.6.new[,1]))*10

sd((keep.track.N.i.4.new[,1]))*10
sd((keep.track.no.N.i.4.new[,1]))*10
sd((keep.track.straw.4.new[,1]))*10
sd((keep.track.CB.4.new[,1]))*10

sd((keep.track.N.i.2.new[,1]))*10
sd((keep.track.no.N.i.2.new[,1]))*10
sd((keep.track.straw.2.new[,1]))*10
sd((keep.track.CB.2.new[,1]))*10


#beta 1-----MSE ----------#

mean((keep.track.full.new[,2]-beta.1)^2)*10^5
mean((keep.track.N.i.only.0.0.0.new[,2]-beta.1)^2)*10^5

mean((keep.track.N.i.8.new[,2]-beta.1)^2)*10^5
mean((keep.track.no.N.i.8.new[,2]-beta.1)^2)*10^5
mean((keep.track.straw.8.new[,2]-beta.1)^2)*10^5
mean((keep.track.CB.8.new[,2]-beta.1)^2)*10^5

mean((keep.track.N.i.6.new[,2]-beta.1)^2)*10^5
mean((keep.track.no.N.i.6.new[,2]-beta.1)^2)*10^5
mean((keep.track.straw.6.new[,2]-beta.1)^2)*10^5
mean((keep.track.CB.6.new[,2]-beta.1)^2)*10^5

mean((keep.track.N.i.4.new[,2]-beta.1)^2)*10^5
mean((keep.track.no.N.i.4.new[,2]-beta.1)^2)*10^5
mean((keep.track.straw.4.new[,2]-beta.1)^2)*10^5
mean((keep.track.CB.4.new[,2]-beta.1)^2)*10^5

mean((keep.track.N.i.2.new[,2]-beta.1)^2)*10^5
mean((keep.track.no.N.i.2.new[,2]-beta.1)^2)*10^5
mean((keep.track.straw.2.new[,2]-beta.1)^2)*10^5
mean((keep.track.CB.2.new[,2]-beta.1)^2)*10^5

#beta 1-----bias ----------#

mean((keep.track.full.new[,2]-beta.1))*10^4
mean((keep.track.N.i.only.0.0.0.new[,2]-beta.1))*10^4

mean((keep.track.N.i.8.new[,2]-beta.1))*10^4
mean((keep.track.no.N.i.8.new[,2]-beta.1))*10^4
mean((keep.track.straw.8.new[,2]-beta.1))*10^4
mean((keep.track.CB.8.new[,2]-beta.1))*10^4

mean((keep.track.N.i.6.new[,2]-beta.1))*10^4
mean((keep.track.no.N.i.6.new[,2]-beta.1))*10^4
mean((keep.track.straw.6.new[,2]-beta.1))*10^4
mean((keep.track.CB.6.new[,2]-beta.1))*10^4

mean((keep.track.N.i.4.new[,2]-beta.1))*10^4
mean((keep.track.no.N.i.4.new[,2]-beta.1))*10^4
mean((keep.track.straw.4.new[,2]-beta.1))*10^4
mean((keep.track.CB.4.new[,2]-beta.1))*10^4

mean((keep.track.N.i.2.new[,2]-beta.1))*10^4
mean((keep.track.no.N.i.2.new[,2]-beta.1))*10^4
mean((keep.track.straw.2.new[,2]-beta.1))*10^4
mean((keep.track.CB.2.new[,2]-beta.1))*10^4

#beta 1-----sd ----------#
sd((keep.track.full.new[,2]))*10^2
sd((keep.track.N.i.only.0.0.0.new[,2]))*10^2

sd((keep.track.N.i.8.new[,2]))*10^2
sd((keep.track.no.N.i.8.new[,2]))*10^2
sd((keep.track.straw.8.new[,2]))*10^2
sd((keep.track.CB.8.new[,2]))*10^2

sd((keep.track.N.i.6.new[,2]))*10^2
sd((keep.track.no.N.i.6.new[,2]))*10^2
sd((keep.track.straw.6.new[,2]))*10^2
sd((keep.track.CB.6.new[,2]))*10^2

sd((keep.track.N.i.4.new[,2]))*10^2
sd((keep.track.no.N.i.4.new[,2]))*10^2
sd((keep.track.straw.4.new[,2]))*10^2
sd((keep.track.CB.4.new[,2]))*10^2

sd((keep.track.N.i.2.new[,2]))*10^2
sd((keep.track.no.N.i.2.new[,2]))*10^2
sd((keep.track.straw.2.new[,2]))*10^2
sd((keep.track.CB.2.new[,2]))*10^2




#beta 2-----MSE ----------#

mean((keep.track.full.new[,3]-beta.2)^2)*10^5
mean((keep.track.N.i.only.0.0.0.new[,3]-beta.2)^2)*10^5

mean((keep.track.N.i.8.new[,3]-beta.2)^2)*10^5
mean((keep.track.no.N.i.8.new[,3]-beta.2)^2)*10^5
mean((keep.track.straw.8.new[,3]-beta.2)^2)*10^5
mean((keep.track.CB.8.new[,3]-beta.2)^2)*10^5

mean((keep.track.N.i.6.new[,3]-beta.2)^2)*10^5
mean((keep.track.no.N.i.6.new[,3]-beta.2)^2)*10^5
mean((keep.track.straw.6.new[,3]-beta.2)^2)*10^5
mean((keep.track.CB.6.new[,3]-beta.2)^2)*10^5

mean((keep.track.N.i.4.new[,3]-beta.2)^2)*10^5
mean((keep.track.no.N.i.4.new[,3]-beta.2)^2)*10^5
mean((keep.track.straw.4.new[,3]-beta.2)^2)*10^5
mean((keep.track.CB.4.new[,3]-beta.2)^2)*10^5

mean((keep.track.N.i.2.new[,3]-beta.2)^2)*10^5
mean((keep.track.no.N.i.2.new[,3]-beta.2)^2)*10^5
mean((keep.track.straw.2.new[,3]-beta.2)^2)*10^5
mean((keep.track.CB.2.new[,3]-beta.2)^2)*10^5

#beta 2-----bias ----------#

mean((keep.track.full.new[,3]-beta.2))*10^4
mean((keep.track.N.i.only.0.0.0.new[,3]-beta.2))*10^4

mean((keep.track.N.i.8.new[,3]-beta.2))*10^4
mean((keep.track.no.N.i.8.new[,3]-beta.2))*10^4
mean((keep.track.straw.8.new[,3]-beta.2))*10^4
mean((keep.track.CB.8.new[,3]-beta.2))*10^4

mean((keep.track.N.i.6.new[,3]-beta.2))*10^4
mean((keep.track.no.N.i.6.new[,3]-beta.2))*10^4
mean((keep.track.straw.6.new[,3]-beta.2))*10^4
mean((keep.track.CB.6.new[,3]-beta.2))*10^4

mean((keep.track.N.i.4.new[,3]-beta.2))*10^4
mean((keep.track.no.N.i.4.new[,3]-beta.2))*10^4
mean((keep.track.straw.4.new[,3]-beta.2))*10^4
mean((keep.track.CB.4.new[,3]-beta.2))*10^4

mean((keep.track.N.i.2.new[,3]-beta.2))*10^4
mean((keep.track.no.N.i.2.new[,3]-beta.2))*10^4
mean((keep.track.straw.2.new[,3]-beta.2))*10^4
mean((keep.track.CB.2.new[,3]-beta.2))*10^4

#beta 2-----sd ----------#
sd((keep.track.full.new[,3]))*10^3
sd((keep.track.N.i.only.0.0.0.new[,3]))*10^3

sd((keep.track.N.i.8.new[,3]))*10^3
sd((keep.track.no.N.i.8.new[,3]))*10^3
sd((keep.track.straw.8.new[,3]))*10^3
sd((keep.track.CB.8.new[,3]))*10^3

sd((keep.track.N.i.6.new[,3]))*10^3
sd((keep.track.no.N.i.6.new[,3]))*10^3
sd((keep.track.straw.6.new[,3]))*10^3
sd((keep.track.CB.6.new[,3]))*10^3

sd((keep.track.N.i.4.new[,3]))*10^3
sd((keep.track.no.N.i.4.new[,3]))*10^3
sd((keep.track.straw.4.new[,3]))*10^3
sd((keep.track.CB.4.new[,3]))*10^3

sd((keep.track.N.i.2.new[,3]))*10^3
sd((keep.track.no.N.i.2.new[,3]))*10^3
sd((keep.track.straw.2.new[,3]))*10^3
sd((keep.track.CB.2.new[,3]))*10^3
