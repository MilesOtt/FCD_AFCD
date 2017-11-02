#install.packages("foreign")
library(foreign)
#install.packages("amen")
library(amen)


##----------- Reading in data ----------------#
setwd("C:/Users/Miles/Dropbox/Miles Backup/crystal networks/updated dta files")
getwd()

pnom<-as.matrix(read.dta("nomination_matrix_129.dta"), nrow=129, ncol=130)
subj<-as.matrix(read.dta("subjectdata.dta"), nrow=129, ncol=132)
sna_dorm<-as.matrix(read.dta("sndata_long_dorm_only.dta"))

w<-matrix(rep(0,129*129), nrow=129, ncol=129)
	y<-pnom[1:129,2:130] 

	for(i in 1:129){
		for(j in 1:129){
			w[i,j]<-as.numeric(y[i,j])
				}
		}
y.matrix.data.full<-w


keep.track.CB.1<-NULL
keep.track.CB.2<-NULL
keep.track.CB.3<-NULL
keep.track.CB.4<-NULL
keep.track.CB.5<-NULL
keep.track.CB.6<-NULL
keep.track.CB.7<-NULL
keep.track.CB.8<-NULL
keep.track.CB.9<-NULL


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

n.loop<-100
for(iii in 1:n.loop){


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


#--------- CB 2 ------------------------------#

fit.full<-ame(y.matrix.data.2, Xdyad=x.2.matrix.data, Xrow=drink.days.month, 
             model="cbin", intercept=TRUE, odmax=2, print=FALSE, plot=FALSE,
             burn=5000,nscan=50000, odens=25, seed=iii)

params.full<-colMeans(fit.full$BETA)
keep.track.CB.2<-rbind(keep.track.CB.2,params.full)


#--------- CB 4 ------------------------------#
fit.full<-ame(y.matrix.data.4, Xdyad=x.2.matrix.data, Xrow=drink.days.month, 
              model="cbin", intercept=TRUE, odmax=4, print=FALSE, plot=FALSE,
              burn=5000,nscan=50000, odens=25, seed=iii)

params.full<-colMeans(fit.full$BETA)
keep.track.CB.4<-rbind(keep.track.CB.4,params.full)

#--------- CB 6 ------------------------------#
fit.full<-ame(y.matrix.data.6, Xdyad=x.2.matrix.data, Xrow=drink.days.month, 
              model="cbin", intercept=TRUE, odmax=6, print=FALSE, plot=FALSE,
              burn=5000,nscan=50000, odens=25, seed=iii)


params.full<-colMeans(fit.full$BETA)
keep.track.CB.6<-rbind(keep.track.CB.6,params.full)


#--------- CB 8 ------------------------------#
fit.full<-ame(y.matrix.data.8, Xdyad=x.2.matrix.data, Xrow=drink.days.month, 
              model="cbin", intercept=TRUE, odmax=8, print=FALSE, plot=FALSE,
              burn=5000,nscan=50000, odens=25, seed=iii)


params.full<-colMeans(fit.full$BETA)
keep.track.CB.8<-rbind(keep.track.CB.8,params.full)

print(iii)


}
