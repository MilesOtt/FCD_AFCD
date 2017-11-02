#install.packages("foreign")
library(foreign)
library(amen)
#----------- Reading in data ----------------#


pnom<-as.matrix(read.dta("nomination_matrix_129.dta"), nrow=129, ncol=130)
subj<-as.matrix(read.dta("subjectdata.dta"), nrow=129, ncol=132)
sna_dorm<-as.matrix(read.dta("sndata_long_dorm_only.dta"))
sna_try<-(read.dta("urweb129.dta"))
attach(sna_try)
sna_try.m<-as.matrix(read.dta("urweb129.dta"))

w<-matrix(rep(0,129*129), nrow=129, ncol=129)
	y<-pnom[1:129,2:130] 

	for(i in 1:129){
		for(j in 1:129){
			w[i,j]<-as.numeric(y[i,j])
				}
		}
y.matrix.data.full<-w

c.id<-subj[,1]
index.cid<-1:length(c.id)

order.nom<-matrix(rep(0, (length(c.id)*20)), ncol=20)
for (i in index.cid){
	a<-0
	if(SN2P1[customid==c.id[i]] %in% c.id)
	a<-c(a, index.cid[c.id==SN2P1[customid==c.id[i]]])
	
	if(SN2P2[customid==c.id[i]] %in% c.id)
	a<-c(a, index.cid[c.id==SN2P2[customid==c.id[i]]])
	
	if(SN2P3[customid==c.id[i]] %in% c.id)
	a<-c(a, index.cid[c.id==SN2P3[customid==c.id[i]]])

	if(SN2P4[customid==c.id[i]] %in% c.id)
	a<-c(a, index.cid[c.id==SN2P4[customid==c.id[i]]])
	
	if(SN2P5[customid==c.id[i]] %in% c.id)
	a<-c(a, index.cid[c.id==SN2P5[customid==c.id[i]]])
	
	if(SN2P6[customid==c.id[i]] %in% c.id)
	a<-c(a, index.cid[c.id==SN2P6[customid==c.id[i]]])
	
	if(SN2P7[customid==c.id[i]] %in% c.id)
	a<-c(a, index.cid[c.id==SN2P7[customid==c.id[i]]])

	if(SN2P8[customid==c.id[i]] %in% c.id)
	a<-c(a, index.cid[c.id==SN2P8[customid==c.id[i]]])
	
	if(SN2P9[customid==c.id[i]] %in% c.id)
	a<-c(a, index.cid[c.id==SN2P9[customid==c.id[i]]])

	if(SN2P10[customid==c.id[i]] %in% c.id)
	a<-c(a, index.cid[c.id==SN2P10[customid==c.id[i]]])
	
	if(SN2P11[customid==c.id[i]] %in% c.id)
	a<-c(a, index.cid[c.id==SN2P11[customid==c.id[i]]])
	
	if(SN2P12[customid==c.id[i]] %in% c.id)
	a<-c(a, index.cid[c.id==SN2P12[customid==c.id[i]]])
	
	if(SN2P13[customid==c.id[i]] %in% c.id)
	a<-c(a, index.cid[c.id==SN2P13[customid==c.id[i]]])
	
	if(SN2P14[customid==c.id[i]] %in% c.id)
	a<-c(a, index.cid[c.id==SN2P14[customid==c.id[i]]])
	
	if(SN2P15[customid==c.id[i]] %in% c.id)
	a<-c(a, index.cid[c.id==SN2P15[customid==c.id[i]]])
	
	if(SN2P16[customid==c.id[i]] %in% c.id)
	a<-c(a, index.cid[c.id==SN2P16[customid==c.id[i]]])

	if(SN2P17[customid==c.id[i]] %in% c.id)
	a<-c(a, index.cid[c.id==SN2P17[customid==c.id[i]]])

	if(SN2P18[customid==c.id[i]] %in% c.id)
	a<-c(a, index.cid[c.id==SN2P18[customid==c.id[i]]])
	
	if(SN2P19[customid==c.id[i]] %in% c.id)
	a<-c(a, index.cid[c.id==SN2P19[customid==c.id[i]]])
	
	if(SN2P20[customid==c.id[i]] %in% c.id)
	a<-c(a, index.cid[c.id==SN2P20[customid==c.id[i]]])

	if(length(a)>1){
		a<-a[2:length(a)]
		order.nom[i,1:length(a)]<-a
	}
}

iii<-1

#------------#

se.from.observed.info<-function(betas){
x.0.vec<-vectorized.x.0 #global variables
x.1.vec<-vectorized.x.1 #global variables
x.2.vec<-vectorized.x.2 #global variables

logit<-betas[1]+(betas[2]*x.1.vec)+(betas[3]*x.2.vec)
pi.vec<-exp(logit)/(1+exp(logit))
x.matrix.se<-matrix(c(x.0.vec, x.1.vec,x.2.vec), ncol=3)
v.matrix<-matrix(rep(0, (length(x.0.vec)^2)),nrow=length(x.0.vec))
diag(v.matrix)<-pi.vec

info.matrix<-t(x.matrix.se)%*%v.matrix%*%x.matrix.se
vcov<-solve(info.matrix)
se.beta.0<-sqrt(vcov[1,1])
se.beta.1<-sqrt(vcov[2,2])
se.beta.2<-sqrt(vcov[3,3])

se<-c(se.beta.0, se.beta.1, se.beta.2)
return(se)
}

#------------#

q.n<-function(probs){
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
	real.q<-q[n,(m+1)]	
	}
	if (m==0){
		real.q<-prod((1-p))	
	}
return(real.q)
}
#----------#


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
  

  missing.row<-(N.i>=max.here)
  for (q in index[missing.row]){
    n.missing<-0:(n.nodes-max.here-1)
    prob.no.obs[index.no.obs==q]<-prob.no.obs[index.no.obs==q]	
    passer<-c(prob.no.obs[index.no.obs==q],(n.nodes-max.here-1))
    part.3.q<-t(1/choose((max.here+n.missing), max.here))%*%q.n.no.N.i(passer)
    part.3<-part.3+log(part.3.q)
  }	

  
  this<-part.1+part.2+part.3
  return(-this)
}

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
		part.3<-part.3+log(q.n(passer)*(1/choose((max.here+n.missing), max.here)))
	}
	this<-part.1+part.2+part.3
	return(-this)
}

#--------#

loglike.no.N.i<-function(theta){
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


	missing.row<-(N.i>=max.here)
	for (q in index[missing.row]){
		n.missing<-0:(n.nodes-max.here-1)
		exp.v.y.no.obs[index.no.obs==q]<-exp.v.y.no.obs[index.no.obs==q]	
		passer<-c(exp.v.y.no.obs[index.no.obs==q],(n.nodes-max.here-1))
		part.3.q<-t(1/choose((max.here+n.missing), max.here))%*%q.n.no.N.i(passer)
		part.3<-part.3+log(part.3.q)
	}	

	this<-part.1+part.2+part.3
	return(-this)
}
#-------------------------------- End Functions-----------------------------------------------#

#-------------------------------- End Functions-----------------------------------------------#


keep.track.CB.1<-NULL
keep.track.CB.2<-NULL
keep.track.CB.3<-NULL
keep.track.CB.4<-NULL
keep.track.CB.5<-NULL
keep.track.CB.6<-NULL
keep.track.CB.7<-NULL
keep.track.CB.8<-NULL
keep.track.CB.9<-NULL


keep.track.N.i.only.0<-NULL
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

keep.track.no.N.i.1<-NULL
keep.track.no.N.i.2<-NULL
keep.track.no.N.i.3<-NULL
keep.track.no.N.i.4<-NULL
keep.track.no.N.i.5<-NULL
keep.track.no.N.i.6<-NULL
keep.track.no.N.i.7<-NULL
keep.track.no.N.i.8<-NULL
keep.track.no.N.i.9<-NULL



n.nodes<-129


x.0.matrix.data<-matrix(rep(1, (n.nodes*n.nodes)), nrow=n.nodes)

drink.days.month<-as.numeric(subj[,112])

index<-1:n.nodes
index.matrix<-matrix(index, ncol=n.nodes,nrow=n.nodes)
x.1.matrix.data<-matrix(drink.days.month, ncol=n.nodes,nrow=n.nodes)
x.2.matrix.data<-x.1.matrix.data

for(i in 1:n.nodes){
	for(j in 1:n.nodes){
		x.2.matrix.data[i,j]<-abs(drink.days.month[i]-drink.days.month[j])
	}
}


#for(iii in 1:1){


#---- setting up the data--------#


N.i<-rowSums(y.matrix.data.full)
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
		t.c<-order.nom[i,]
		t.c<-t.c[t.c!=0]
		t.c.o<-(t.c)
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


		if(n.c>3){
			y.matrix.data.3[i,t.c.o[4:n.c]]<-0
			y.miss.3[i,y.matrix.data.3[i,]==0]<-1
			y.miss.no.N.i.3[i,y.matrix.data.3[i,]==0]<-1
		}
		if (n.c==3){
			y.miss.no.N.i.3[i,y.matrix.data.3[i,]==0]<-1
		}

		if(n.c>4){
			y.matrix.data.4[i,t.c.o[5:n.c]]<-0
			y.miss.4[i,y.matrix.data.4[i,]==0]<-1
			y.miss.no.N.i.4[i,y.matrix.data.4[i,]==0]<-1
		}
		if (n.c==4){
			y.miss.no.N.i.4[i,y.matrix.data.4[i,]==0]<-1
		}

		if(n.c>5){
			y.matrix.data.5[i,t.c.o[6:n.c]]<-0
			y.miss.5[i,y.matrix.data.5[i,]==0]<-1
			y.miss.no.N.i.5[i,y.matrix.data.5[i,]==0]<-1
		}
		if (n.c==5){
			y.miss.no.N.i.5[i,y.matrix.data.5[i,]==0]<-1
		}

		if(n.c>6){
			y.matrix.data.6[i,t.c.o[7:n.c]]<-0
			y.miss.6[i,y.matrix.data.6[i,]==0]<-1
			y.miss.no.N.i.6[i,y.matrix.data.6[i,]==0]<-1
		}
		if (n.c==6){
			y.miss.no.N.i.6[i,y.matrix.data.6[i,]==0]<-1
		}

		if(n.c>7){
			y.matrix.data.7[i,t.c.o[8:n.c]]<-0
			y.miss.7[i,y.matrix.data.7[i,]==0]<-1
			y.miss.no.N.i.7[i,y.matrix.data.7[i,]==0]<-1
		}
		if (n.c==7){
			y.miss.no.N.i.7[i,y.matrix.data.7[i,]==0]<-1
		}

		if(n.c>8){
			y.matrix.data.8[i,t.c.o[9:n.c]]<-0
			y.miss.8[i,y.matrix.data.8[i,]==0]<-1
			y.miss.no.N.i.8[i,y.matrix.data.8[i,]==0]<-1
		}
		if (n.c==8){
			y.miss.no.N.i.8[i,y.matrix.data.8[i,]==0]<-1
		}

		if(n.c>9){
			y.matrix.data.9[i,t.c.o[10:n.c]]<-0
			y.miss.9[i,y.matrix.data.9[i,]==0]<-1
			y.miss.no.N.i.9[i,y.matrix.data.9[i,]==0]<-1
		}
		if (n.c==9){
			y.miss.no.N.i.9[i,y.matrix.data.9[i,]==0]<-1
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

y.miss.no.N.i.vec.1<-c(y.miss.no.N.i.1[upper.tri(y.miss.5)], y.miss.no.N.i.1[lower.tri(y.miss.5)])
y.miss.no.N.i.vec.2<-c(y.miss.no.N.i.2[upper.tri(y.miss.5)], y.miss.no.N.i.2[lower.tri(y.miss.5)])
y.miss.no.N.i.vec.3<-c(y.miss.no.N.i.3[upper.tri(y.miss.5)], y.miss.no.N.i.3[lower.tri(y.miss.5)])
y.miss.no.N.i.vec.4<-c(y.miss.no.N.i.4[upper.tri(y.miss.5)], y.miss.no.N.i.4[lower.tri(y.miss.5)])
y.miss.no.N.i.vec.5<-c(y.miss.no.N.i.5[upper.tri(y.miss.5)], y.miss.no.N.i.5[lower.tri(y.miss.5)])
y.miss.no.N.i.vec.6<-c(y.miss.no.N.i.6[upper.tri(y.miss.5)], y.miss.no.N.i.6[lower.tri(y.miss.5)])
y.miss.no.N.i.vec.7<-c(y.miss.no.N.i.7[upper.tri(y.miss.5)], y.miss.no.N.i.7[lower.tri(y.miss.5)])
y.miss.no.N.i.vec.8<-c(y.miss.no.N.i.8[upper.tri(y.miss.5)], y.miss.no.N.i.8[lower.tri(y.miss.5)])
y.miss.no.N.i.vec.9<-c(y.miss.no.N.i.9[upper.tri(y.miss.5)], y.miss.no.N.i.9[lower.tri(y.miss.5)])

y.vec.1<-c(y.matrix.data.1[upper.tri(y.miss.5)], y.matrix.data.1[lower.tri(y.miss.5)])
y.vec.2<-c(y.matrix.data.2[upper.tri(y.miss.5)], y.matrix.data.2[lower.tri(y.miss.5)])
y.vec.3<-c(y.matrix.data.3[upper.tri(y.miss.5)], y.matrix.data.3[lower.tri(y.miss.5)])
y.vec.4<-c(y.matrix.data.4[upper.tri(y.miss.5)], y.matrix.data.4[lower.tri(y.miss.5)])
y.vec.5<-c(y.matrix.data.5[upper.tri(y.miss.5)], y.matrix.data.5[lower.tri(y.miss.5)])
y.vec.6<-c(y.matrix.data.6[upper.tri(y.miss.5)], y.matrix.data.6[lower.tri(y.miss.5)])
y.vec.7<-c(y.matrix.data.7[upper.tri(y.miss.5)], y.matrix.data.7[lower.tri(y.miss.5)])
y.vec.8<-c(y.matrix.data.8[upper.tri(y.miss.5)], y.matrix.data.8[lower.tri(y.miss.5)])
y.vec.9<-c(y.matrix.data.9[upper.tri(y.miss.5)], y.matrix.data.9[lower.tri(y.miss.5)])


index.vec<-c(index.matrix[upper.tri(index.matrix)], index.matrix[lower.tri(index.matrix)])

#--------------- End Creating the Data--------------#
set.seed(iii)

#------------N.i only -----------------------------#
if(iii==1){
  
  max.here<-0
  miss.vec<-y.miss.vec.5
  miss.vec[miss.vec==0]<-1
  theta.try.1<-c(0,0,0)
  #theta.try.2<-c(-1,-1,-1)
  #theta.try.3<-c(-.5, 0, -.5)
  params.N.i.only.0<-optim(theta.try.1, loglike.probit, method="BFGS")$par
  
  keep.track.N.i.only.0<-rbind(keep.track.N.i.only.0,params.N.i.only.0)
}

#------------ N.i 2--------------------------------#
max.here<-2
miss.vec<-y.miss.vec.2
theta.try<-c(0,0,0)
params2<-optim(theta.try, loglike.probit, method="BFGS")$par
keep.track.N.i.2<-rbind(keep.track.N.i.2,params2)



#------------ N.i 4--------------------------------#
max.here<-4
miss.vec<-y.miss.vec.4
theta.try<-c(0,0,0)
params<-optim(theta.try, loglike.probit, method="BFGS")$par
keep.track.N.i.4<-rbind(keep.track.N.i.4,params)


#------------ N.i 6--------------------------------#
max.here<-6
miss.vec<-y.miss.vec.6
theta.try<-c(0,0,0)
params<-optim(theta.try, loglike.probit, method="BFGS")$par
keep.track.N.i.6<-rbind(keep.track.N.i.6,params)

##------------ N.i 8--------------------------------#
max.here<-8
miss.vec<-y.miss.vec.8
theta.try<-c(0,0,0)
params<-optim(theta.try, loglike.probit, method="BFGS")$par
keep.track.N.i.8<-rbind(keep.track.N.i.8,params)

#--------- Full Data ------------------------------#


fit.full<-glm(y.vecs.full~x.vecs[,2]+x.vecs[,3], family="binomial"(link="probit"))
params.full.real<-fit.full$coef
keep.track.full<-rbind(keep.track.full,params.full.real)


#--------- Straw 2 ------------------------------#
fit.full<-glm(y.vec.2~x.vecs[,2]+x.vecs[,3], family="binomial"(link="probit"))
params.full<-fit.full$coef
keep.track.straw.2<-rbind(keep.track.straw.2,params.full)


#--------- Straw 4 ------------------------------#
fit.full<-glm(y.vec.4~x.vecs[,2]+x.vecs[,3], family="binomial"(link="probit"))
params.full<-fit.full$coef
keep.track.straw.4<-rbind(keep.track.straw.4,params.full)


#--------- Straw 6 ------------------------------#
fit.full<-glm(y.vec.6~x.vecs[,2]+x.vecs[,3], family="binomial"(link="probit"))
params.full<-fit.full$coef
keep.track.straw.6<-rbind(keep.track.straw.6,params.full)


#--------- Straw 8 ------------------------------#
fit.full<-glm(y.vec.8~x.vecs[,2]+x.vecs[,3], family="binomial"(link="probit"))
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

fit.2<-ame(y.matrix.data.2, Xdyad=x.2.matrix.data, Xrow=drink.days.month, 
              model="cbin", intercept=TRUE, odmax=2, print=FALSE, plot=FALSE,
              burn=2000,nscan=10000, odens=5, seed=iii)


#--------- CB 4 ------------------------------#
fit.4<-ame(y.matrix.data.4, Xdyad=x.2.matrix.data, Xrow=drink.days.month, 
              model="cbin", intercept=TRUE, odmax=4, print=FALSE, plot=FALSE,
              burn=2000,nscan=10000, odens=5, seed=iii)


#--------- CB 6 ------------------------------#
fit.6<-ame(y.matrix.data.6, Xdyad=x.2.matrix.data, Xrow=drink.days.month, 
              model="cbin", intercept=TRUE, odmax=6, print=FALSE, plot=FALSE,
              burn=2000,nscan=10000, odens=5, seed=iii)


#--------- CB 8 ------------------------------#
fit.8<-ame(y.matrix.data.8, Xdyad=x.2.matrix.data, Xrow=drink.days.month, 
              model="cbin", intercept=TRUE, odmax=8, print=FALSE, plot=FALSE,
              burn=2000,nscan=10000, odens=5, seed=iii)

