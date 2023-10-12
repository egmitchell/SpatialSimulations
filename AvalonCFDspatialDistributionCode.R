######### Code to simulate Mistaken Point Communities


library(spatstat)
library(mclust)

Mclust(marks(d.species1[[1]]))$parameters
Mclust(marks(d.species1[[2]]))$parameters
Mclust(log(marks(d.species1[[3]])))$parameters


Mclust((pect.data$width)/10)$parameters
Mclust(bra.data$Width/10)$parameters
Mclust(log(fra.data$Width/10))$parameters


#going for majority, ignoring small numbers

win.sim.d<-square(0,1)

#orientations
orients<-sample(seq(0,360,1),10,replace=TRUE)

#size
size.bra<-rnorm(10,12,3)
size.pect<-rnorm(10,67,19)
size.fra<-exp(rlnorm(10,0.181,0.341))

#
#kppm.fra<-(kppm(unmark(d.species1[[3]]), clusters = "Thomas",statistic="pcf"))

#simulate the points

res.sim.d<-list()

par(mfrow=c(5,5))
for(i in 1:25)
{
	sim.fra1<-rThomas(55.92535055,0.06704791,0.3112936,win=win.sim.d,saveparents = TRUE)
	lengths1<-rlnorm(sim.fra1$n,2.176,0.440)
	widths1<-2.6+0.25*lengths1
	marks(sim.fra1)<-(cbind(lengths1,widths1,sample(seq(0,360,1),sim.fra1$n,replace=TRUE)))
	plot(sim.fra1,which.marks="lengths1",bg="cyan",col="grey",main=paste("Simulation",i))
	
	sim.bra1<-rpoispp(3,window=win.sim.d)
	lengths1<-rnorm(sim.bra1$n,12.88,2.95)
	widths1<-1.5+0.95*lengths1
	marks(sim.bra1)<-(cbind(lengths1,widths1,sample(seq(0,360,1),sim.bra1$n,replace=TRUE)))
	plot(sim.bra1,which.marks="lengths1",bg="darkgreen",add=TRUE)

	sim.pect1<-rpoispp(2,window=win.sim.d)
	lengths1<-rnorm(sim.pect1$n,17.34,6.08)
	widths1<-4.1+0.14*lengths1
	marks(sim.pect1)<-(cbind(lengths1,widths1,sample(seq(0,360,1),sim.pect1$n,replace=TRUE)))
	plot(sim.pect1,which.marks="lengths1",bg="navy",add=TRUE)

	sim.d.1<-superimpose(fra=sim.fra1,bra=sim.bra1,pect=sim.pect1)
	res.sim.d[[i]]<-sim.d.1
	sim.d.1.table<-cbind(coords(sim.d.1),marks(sim.d.1))
	colnames(sim.d.1.table)<-c("x","y","Length","Width","Orientation","species")

	filename1<-paste("D://simDa_",paste(i),".txt",sep="")
	write.table(sim.d.1.table,filename1,sep="\t",row.names=FALSE)
	
}

sim.lmp.sp<-list()

win.size<-10#n
win.sim.lmp<-square(win.size) #100m2

#kppm(unmark(rescale(lmp.species1[[1]],100)),cluster="Thomas",rmax=0.4)
#kppm(unmark(rescale(lmp.species1[[2]],100)),cluster="Thomas",rmax=0.2)
#kppm(unmark(rescale(lmp.species1[[3]],100)),cluster="Thomas",rmax=0.2)

#nb input parameter needs to be 2x point where closers pcf=1

#rSpeciesTC(npar,nclust,r1,parMean,parStd,kidsMean,kidsStd)
#rSpeciesTC(round(win.size^2*(meanDensity)(propLarge),propSmall/propLarge,numberPerCluster,mean2,var2,mean1,varience1)
#rSpeciesRC(ppp.par1,nclust2,r2,newkidsMean,newkidsStd)

#fit tc to beo culmo and ost


#Beothukis 
sim.lmp.sp[[1]]<-rSpeciesTC(round(win.size^2*7.46/2),round(0.5544219/0.4455781),1.26,10.100947,sqrt(23.531954),4.388379,sqrt(1.887486))
sim.lmp.sp[[1]]<-subset(sim.lmp.sp[[1]],marks(sim.lmp.sp[[1]])>1)

#Culmofrons
sim.lmp.sp[[2]]<-rSpeciesTC(round(win.size^2*1.963233/2),round(0.4375468/0.5624532),0.1553098,8.979649,sqrt(18.134112),2.958462,sqrt(0.3482764))
sim.lmp.sp[[2]]<-subset(sim.lmp.sp[[2]],marks(sim.lmp.sp[[2]])>1)


#ostrich
#hpTC
hpBeo<-density(sim.lmp.sp[[1]],0.1)/100
sim.lmp.sp[[3]]<-rpoint(3.69*win.size^2,hpBeo,win=win.sim.lmp)
marks(sim.lmp.sp[[3]])<-rnorm(sim.lmp.sp[[3]]$n,3.968519,sqrt(2.960676))
sim.lmp.sp[[3]]<-subset(sim.lmp.sp[[3]],marks(sim.lmp.sp[[3]])>1)

#Charniodiscus
sim.lmp.sp[[4]]<-rpoispp(1,win=win.sim.lmp)
marks(sim.lmp.sp[[4]])<-rnorm(sim.lmp.sp[[4]]$n,9.4,5.559)
sim.lmp.sp[[4]]<-subset(sim.lmp.sp[[4]],marks(sim.lmp.sp[[4]])>1)


#fronds
sim.lmp.sp[[5]]<-rpoispp(3,win=win.sim.lmp)
marks(sim.lmp.sp[[5]])<-rnorm(sim.lmp.sp[[5]]$n,5.8,3.7)
sim.lmp.sp[[5]]<-subset(sim.lmp.sp[[5]],marks(sim.lmp.sp[[5]])>1)


########

names(sim.lmp.sp)<-c("Beothukis","Culmofrons","OstrichFeather","Charniodiscus","Fronds")

#add dimensions


sim.lmp.1<-superimpose(Beothukis=sim.lmp.sp[[1]],Culmofrons=sim.lmp.sp[[2]],OstrichFeather=sim.lmp.sp[[3]],Charniodiscus=sim.lmp.sp[[4]],Fronds=sim.lmp.sp[[5]])
	sim.lmp.1.table<-cbind(coords(sim.lmp.1),marks(sim.lmp.1))
	colnames(sim.lmp.1.table)<-c("x","y","Length","species")
	sim.lmp.1.table.pos<-sim.lmp.1.table[sim.lmp.1.table[,3]>1,]
	filename1<-paste("D://simLMP2_",paste(i),".txt",sep="")
	write.table(sim.lmp.1.table.pos,filename1,sep="\t",row.names=FALSE)


#subsample to 1m2 and ouput

ppp2<-sim.lmp.1

par(mfrow=c(5,5))
count.lmp<-c(1:25)
ppp2<-(rescale(ppp2,1.17))
for(i in 1:25)
{
	maxWindSize<-8
	coordsMin<-sample(maxWindSize-1,2)
	sim.lmp.bi.s1<-ppp(coords(ppp2)[,1],coords(ppp2)[,2],marks=marks(ppp2),window = owin(c(coordsMin[1],coordsMin[1]+1),c(coordsMin[2],coordsMin[2]+1)))
	sim.lmp.bi.table<-cbind(coords(sim.lmp.bi.s1)[,1]-coordsMin[1],coords(sim.lmp.bi.s1)[,2]-coordsMin[2],marks(sim.lmp.bi.s1))
	colnames(sim.lmp.bi.table)<-c("x","y","Length","species")
	count.lmp[i]<-sim.lmp.bi.s1$n
	res1<-add.fossil.dim.lmp(sim.lmp.bi.table)
	filename1<-paste("D://simLMPDim_",paste(i),".txt",sep="")
	write.table(res1,filename1,sep="\t",row.names=FALSE)
}



##########################################################################
##below all changed
add.fossil.dim.lmp<-function(sim.table)
{
	dim.table<-matrix(0,nrow(sim.table),4)
	colnames(dim.table)<-c("Frond Length", "Frond Width", "Stem Length", "Stem Width")

	for(i in 1:nrow(sim.table))
	{
		length1<-sim.table[i,3]

		if(sim.table[i,4]=="Beothukis")
		{
			dim.table[i,1]<-length1
			dim.table[i,2]<-0.2555*length1+1.01
			#print(length1)
		}

		else if(sim.table[i,4]=="Culmofrons")
		{

				dim.table[i,1]<-0.63*length1
				dim.table[i,2]<-0.4975*dim.table[i,1]+ 0.31892
				dim.table[i,3]<-length1-dim.table[i,1]
				dim.table[i,4]<-0.18*dim.table[i,3] #these are based off measurements off my unpublished measurements for a limited subset of data. 

					#print(dim.table[i,1])

		}
		
		else if(sim.table[i,4]=="OstrichFeather")
		{
				dim.table[i,1]<-0.72*length1
				dim.table[i,2]<-0.4552*dim.table[i,1]+0.401 ####
				dim.table[i,3]<-length1-dim.table[i,1]
				dim.table[i,4]<-0.43*dim.table[i,3] #these are based off measurements off my unpublished measurements for a limited subset of data. 
		}
		else if(sim.table[i,4]=="Charniodiscus")
		{

				dim.table[i,1]<-0.86*length1
				dim.table[i,2]<-0.483*dim.table[i,1]+0.34007
				dim.table[i,3]<-length1-dim.table[i,1]
				dim.table[i,4]<-0.66*dim.table[i,3] #these are based off measurements off my unpublished measurements for a limited subset of data. 

		}

		else if(sim.table[i,4]=="Fronds")
		{ 
				dim.table[i,1]<-0.78*length1
				dim.table[i,2]<-0.459*dim.table[i,1]+0.40376
				dim.table[i,3]<-length1-dim.table[i,1]
				dim.table[i,4]<-0.145*dim.table[i,3] #these are based off measurements off my unpublished measurements for a limited subset of data. 

		}


	}

	res.lmp.sim.dem<-cbind(sim.table,dim.table)
	
	return(res.lmp.sim.dem)
}

sim.lmp.f.ab<-add.fossil.dim.lmp(sim.lmp.1.table.pos)

sim.lmp.f.ab<-

add.fossil.dim.lmp(sim.lmp.1.table.pos[690:695,])


filename1<-paste("C://simLMPDim_10thFeb_",paste(i),".txt",sep="")
write.table(sim.lmp.f.ab,filename1,sep="\t",row.names=FALSE)

#beo
width = 0.2555*length+1.01
Mclust(marks(lmp.species1$`Charnia II`))$parameters

#Beothukis
lengths1<-rlnorm(sim.fra1$n,2.176,0.440)

lengths1<-sample()
frondWidth<-0.255*lengths1
#Culmofrons
frondlength<-0.41*lengths1
frondwidth<-0.24*lengths1
stemLength1<-0.59*length1
stemwidth1<-0.18*stemLength1 #these are based off measurements off the single holotype since stem width isn't reliably measured

#ostrichs
frondlength<-0.68*lengths1+ 1.8975
frondwidth<-0.39*lengths1+ 2.1503
stemLength1<-0.32*length1
stemwidth1<-0.43*stemLength1 #these are based off measurements off the single holotype since stem width isn't reliably measured

#charnios heights
lmp.cha.h<-c(19.8,12.4,15.3,3.8,25.6,22.4,2.8,7.7,4.6,9.3,16.5)
frondlength<-0.68*lengths1
frondwidth<-0.35*lengths1
stemLength1<-0.32*length1
stemwidth1<-0.66*stemLength1 #these are based off measurements off the single holotype since stem width isn't reliably measured

#fronds
lmp.frond.h<-c(5.8,3.6,6.9,14.6,3.4,8.1,9.4,6.8,3.7,1.7,3.2)
frondlength1<-0.78*lengths1
frondWidth<-0.76*lengths1
stemLength1<-length1-frondlength1
stemwidth1<-0.145*stemLength1 - 1.24 #these are based off measurements three holotype since stem width isn't reliably measured




##ChaP
##########
hp.e.all<-density(superimpose(sim.e.sp[[4]],sim.e.sp[[5]],sim.e.sp[[6]],sim.e.sp[[7]]),0.5)
all.pts<-superimpose(sim.e.sp[[4]],sim.e.sp[[5]],sim.e.sp[[6]],sim.e.sp[[7]])
bi.e.all.pts<-round(win.size^2*summary(all.pts)$intensity*0.5)

bi.chaP<-rpoint(bi.e.all.pts*sim.e.sp[[4]]$n/all.pts$n,hp.e.all,win=win.sim.e)
marks(bi.chaP)<-rnorm(bi.chaP$n,3.8,sqrt(0.52))
sim.e.sp.bi[[4]]<-superimpose(sim.e.sp[[4]],bi.chaP)



ppp.par1<-rSpeciesTC(round(win.size^2*3.37779*0.1723134),round(0.369313/0.1723134),1,13.5,sqrt(26.5),7.02,sqrt(3.95))
sim.e.sp[[4]]<-rSpeciesRC(ppp.par1,2.5,0.5,3.8,sqrt(0.52))




## code for E sims below

#Thectardis
sim.lmp.sp[[1]]<-rSpeciesTC(round(win.size^2*0.2216398*0.2571073),3,2,11.45,sqrt(0.1523),7.23,sqrt(3.47))

#Bradgatia
sim.lmp.sp[[2]]<-rSpeciesTC(round(win.size^2*2.003623*0.1391414),8,2,8.85,sqrt(1.525866),5.76,sqrt(1.525866))


ppp1<-sim.e.bi
#par(mfrow=c(5,5))
for(i in 1:25)
{
	maxWindSize<-10
	coordsMin<-sample(maxWindSize-1,2)
	sim.e.bi.s1<-ppp(coords(ppp1)[,1],coords(ppp1)[,2],marks=marks(ppp1),window = owin(c(coordsMin[1],coordsMin[1]+1),c(coordsMin[2],coordsMin[2]+1)))
		
	sim.e.bi.table<-cbind(coords(sim.e.bi.s1)[,1]-coordsMin[1],coords(sim.e.bi.s1)[,2]-coordsMin[2],marks(sim.e.bi.s1))
	colnames(sim.e.bi.table)<-c("x","y","Length","species")
	print(i)
	res1<-add.fossil.dim(sim.e.bi.table)
	filename1<-paste("C:/simEbivDim4_",paste(i),".txt",sep="")
	write.table(res1,filename1,sep="\t",row.names=FALSE)
}


sim.e.bi<-superimpose(Thectardis=sim.e.sp.bi[[1]],Bradgatia=sim.e.sp.bi[[2]],Beothukis=sim.e.sp.bi[[3]],CharnioP=sim.e.sp.bi[[4]],CharnioS=sim.e.sp.bi[[5]],Fractofusus=sim.e.sp.bi[[6]],Primocandelabrum=sim.e.sp.bi[[7]],Plumeropriscum=sim.e.sp.bi[[8]])
	
	sim.e.bi.table.all<-cbind(coords(sim.e.bi)[,1],coords(sim.e.bi)[,2],marks(sim.e.bi))
	colnames(sim.e.bi.table.all)<-c("x","y","Length","species")
	print(i)
	res1.all<-add.fossil.dim(sim.e.bi.table.all)
	filename2<-paste("C://simEbivDim4_ALL_",paste(i),".txt",sep="")
	write.table(res1.all,filename2,sep="\t",row.names=FALSE)


#simulate E
  
#create joint distributions where needed
sim.e.sp<-list()

win.size<-10#n
win.sim.e<-square(win.size)

#nb input parameter needs to be 2x point where closers pcf=1

#Thectardis
sim.e.sp[[1]]<-rSpeciesTC(round(win.size^2*0.2216398*0.2571073),3,2,11.45,sqrt(0.1523),7.23,sqrt(3.47))

#Bradgatia
sim.e.sp[[2]]<-rSpeciesTC(round(win.size^2*2.003623*0.1391414),8,2,8.85,sqrt(1.525866),5.76,sqrt(1.525866))

#Beothukis #nb uses bradgatia so incorporates the bivariates
hpBack<-density(sim.e.sp[[2]],win=win.sim.e,0.5)
sim.e.sp[[3]]<-rpoint(round(win.size^2*0.7447096),hpBack,win=win.sim.e)
marks(sim.e.sp[[3]])<-rnorm(sim.e.sp[[3]]$n,8.5,sqrt(12.61))

#ChaP

ppp.par1<-rSpeciesTC(round(win.size^2*3.37779*0.1723134),round(0.369313/0.1723134),1,13.5,sqrt(26.5),7.02,sqrt(3.95))
sim.e.sp[[4]]<-rSpeciesRC(ppp.par1,2.5,0.5,3.8,sqrt(0.52))

#ChaS


ppp.par2<-rSpeciesTC(round(win.size^2*1.516016*0.2275827),round(0.3099841/0.2275827),2,7.6,sqrt(8.94),3.76,sqrt(0.57))
sim.e.sp[[5]]<-rSpeciesRC(ppp.par2,1.5,0.34,2.65,sqrt(0.106))
plot(pcf(rescale(ic.species3$E_CharniodiscusS,100),stoyan=0.8),ylim=c(0,2))
plot(pcf(sim.e.sp[[5]],stoyan=0.8),add=TRUE,col=2)

#Fractossumm

ppp.par3<-rSpeciesTC(round(win.size^2*13.27179*0.159),ceiling(0.42/0.159),2.5,16.61,sqrt(29.3),10.78,sqrt(8.12))
sim.e.sp[[6]]<-rSpeciesRC(ppp.par3,1.5,0.34,7.21,sqrt(2.93))

plot(pcf(rescale(ic.species3$E_Fractofusus,100),stoyan=0.8),ylim=c(0,2))
plot(pcf(sim.e.sp[[6]],stoyan=0.8),add=TRUE,col=2)


#Primo

ppp.par4<-rSpeciesTC(round(win.size^2*1.852908*0.243),ceiling(0.426/0.243),4,7.1,sqrt(7.54),3.7,sqrt(1.03))
sim.e.sp[[7]]<-rSpeciesRC(ppp.par4,1.5,0.34,2.1,sqrt(0.2288))

plot(pcf(rescale(ic.species3$E_Primo1,100),stoyan=0.8),ylim=c(0,2))
plot(pcf(sim.e.sp[[7]],stoyan=0.8),add=TRUE,col=2)

#Plum

sim.e.sp[[8]]<-rSpeciesTC(round(win.size^2*0.6383*0.239),ceiling((1-0.239)/0.239),3,8.4,sqrt(11.25),3.19,sqrt(1.03))

plot(pcf(rescale(ic.species3$E_Primo2,100),stoyan=0.8),ylim=c(0,2))
plot(pcf(sim.e.sp[[8]],stoyan=0.8),add=TRUE,col=2)

#par(mfrow=c(5,5))
for(i in 1:25)
{
	maxWindSize<-10
	coordsMin<-sample(maxWindSize-1,2)
	sim.e.bi.s1<-ppp(coords(ppp1)[,1],coords(ppp1)[,2],marks=marks(ppp1),window = owin(c(coordsMin[1],coordsMin[1]+1),c(coordsMin[2],coordsMin[2]+1)))
		
	sim.e.bi.table<-cbind(coords(sim.e.bi.s1)[,1]-coordsMin[1],coords(sim.e.bi.s1)[,2]-coordsMin[2],marks(sim.e.bi.s1))
	colnames(sim.e.bi.table)<-c("x","y","Length","species")
	print(i)
	res1<-add.fossil.dim(sim.e.bi.table)
	filename1<-paste("C://Users/emily/Dropbox/Projects/d_susannaCFDcommunities/simEbivDim4_",paste(i),".txt",sep="")
	write.table(res1,filename1,sep="\t",row.names=FALSE)
}


sim.e.bi<-superimpose(Thectardis=sim.e.sp.bi[[1]],Bradgatia=sim.e.sp.bi[[2]],Beothukis=sim.e.sp.bi[[3]],CharnioP=sim.e.sp.bi[[4]],CharnioS=sim.e.sp.bi[[5]],Fractofusus=sim.e.sp.bi[[6]],Primocandelabrum=sim.e.sp.bi[[7]],Plumeropriscum=sim.e.sp.bi[[8]])
	
	sim.e.bi.table.all<-cbind(coords(sim.e.bi)[,1],coords(sim.e.bi)[,2],marks(sim.e.bi))
	colnames(sim.e.bi.table.all)<-c("x","y","Length","species")
	print(i)
	res1.all<-add.fossil.dim(sim.e.bi.table.all)
	filename2<-paste("C://Users/emily/Dropbox/Projects/d_susannaCFDcommunities/simEbivDim4_ALL_",paste(i),".txt",sep="")
	write.table(res1.all,filename2,sep="\t",row.names=FALSE)




#Subsampling to a meter squared
ppp1<-sim.e.bi



###univariate

	sim.e.1<-superimpose(Thectardis=sim.e.sp[[1]],Bradgatia=sim.e.sp[[2]],Beothukis=sim.e.sp[[3]],CharnioP=sim.e.sp[[4]],CharnioS=sim.e.sp[[5]],Fractofusus=sim.e.sp[[6]],Primocandelabrum=sim.e.sp[[7]],Plumeropriscum=sim.e.sp[[8]])
	sim.e.1.table<-cbind(coords(sim.e.1),marks(sim.e.1))
	colnames(sim.e.1.table)<-c("x","y","Length","species")

	filename1<-paste("D://Dropbox/Projects/d_susannaCFDcommunities/simE_",paste(i),".txt",sep="")
	write.table(sim.e.1.table,filename1,sep="\t",row.names=FALSE)


####not sure what doing
#### bivariate
sim.e.sp.bi<-sim.e.sp


#generate hp
hp.e.all<-density(superimpose(sim.e.sp[[4]],sim.e.sp[[5]],sim.e.sp[[6]],sim.e.sp[[7]]),0.5)
all.pts<-superimpose(sim.e.sp[[4]],sim.e.sp[[5]],sim.e.sp[[6]],sim.e.sp[[7]])
bi.e.all.pts<-round(win.size^2*summary(all.pts)$intensity*0.5)

bi.fea1<-rpoint(bi.e.all.pts*sim.e.sp[[7]]$n/all.pts$n,hp.e.all,win=win.sim.e)
marks(bi.fea1)<-rnorm(bi.fea1$n,2.1,sqrt(0.28))
sim.e.sp.bi[[7]]<-superimpose(sim.e.sp[[7]],bi.fea1)

bi.fea2<-rpoint(bi.e.all.pts*sim.e.sp[[8]]$n/all.pts$n,hp.e.all,win=win.sim.e)
marks(bi.fea2)<-rnorm(bi.fea2$n,3.1,sqrt(1.03))
sim.e.sp.bi[[8]]<-superimpose(sim.e.sp[[8]],bi.fea2)

bi.chaP<-rpoint(bi.e.all.pts*sim.e.sp[[4]]$n/all.pts$n,hp.e.all,win=win.sim.e)
marks(bi.chaP)<-rnorm(bi.chaP$n,3.8,sqrt(0.52))
sim.e.sp.bi[[4]]<-superimpose(sim.e.sp[[4]],bi.chaP)

bi.chaS<-rpoint(bi.e.all.pts*sim.e.sp[[5]]$n/all.pts$n,hp.e.all,win=win.sim.e)
marks(bi.chaS)<-rnorm(bi.chaS$n,2.65,sqrt(0.11))
sim.e.sp.bi[[5]]<-superimpose(sim.e.sp[[5]],bi.chaS)

sim.e.bi<-superimpose(Thectardis=sim.e.sp.bi[[1]],Bradgatia=sim.e.sp.bi[[2]],Beothukis=sim.e.sp.bi[[3]],CharnioP=sim.e.sp.bi[[4]],CharnioS=sim.e.sp.bi[[5]],Fractofusus=sim.e.sp.bi[[6]],Primocandelabrum=sim.e.sp.bi[[7]],Plumeropriscum=sim.e.sp.bi[[8]])
sim.e.bi.table<-cbind(coords(sim.e.bi),marks(sim.e.bi))
colnames(sim.e.bi.table)<-c("x","y","Length","species")

i<-1

filename1<-paste("D://Dropbox/Projects/d_susannaCFDcommunities/simE_",paste(i),".txt",sep="")
write.table(sim.e.bi.table,filename1,sep="\t",row.names=FALSE)



######

rSpeciesHP<-function(den1,r1,mean1,std1)
{
	hp<-density(rpoints(round(win.size^2*0.7447096)),0.5)


}

rSpeciesHP(58,0.5,

rSpeciesTC<-function(npar,nclust,r1,parMean,parStd,kidsMean,kidsStd)
{
	win.test<-win.sim.e#Window(ppp.par)
	ppp.par1<-(rpoint(npar,win=win.test))
	par.coords<-coords(ppp.par1)
	print(nrow(par.coords))
	marks(ppp.par1)<-rnorm(ppp.par1$n,parMean,parStd)
	res.hc.x<-list()
	for(i in 1:nrow(par.coords))
	{
		res.hc.x[[i]]<-runifdisc(nclust, r1, centre=c(par.coords[i,1],par.coords[i,2]))
	}
	res.hc.z2<-res.hc.x[[1]]
	for(i in 2:(length(res.hc.x)-1))
	{
	#
		res.hc.z2<-superimpose(res.hc.x[[i]],res.hc.z2,W=win.test)
		#print(table(duplicated(res.hc.z2)))

	}
	ppp.kids<-res.hc.z2
	print(ppp.kids$n)
	marks(ppp.kids)<-rnorm(ppp.kids$n,kidsMean,kidsStd)
	
	ppp.res<-superimpose(ppp.par1,ppp.kids)#res.hc.y[[i]]#res.hc.z2

	par(mfrow=c(1,3))
	plot(ppp.par1,legend=FALSE)
	plot(ppp.kids,legend=FALSE)
	plot(ppp.res,legend=FALSE)
	return(ppp.res)

}

test1<-rSpeciesTC(10,3,2,5,1,2,0.1) 

rSpeciesRC<-function(ppp.par1,nclust2,r2,newkidsMean,newkidsStd)
{
	win.test<-win.sim.e#Window(ppp.par1)
	par.coords<-coords(ppp.par1)
	res.hc.x<-list()
	for(i in 1:nrow(par.coords))
	{
		res.hc.x[[i]]<-runifdisc(nclust2, r2, centre=c(par.coords[i,1],par.coords[i,2]))
	}
	res.hc.z2<-res.hc.x[[1]]
	for(i in 2:(length(res.hc.x)-1))
	{
	#
		res.hc.z2<-superimpose(res.hc.x[[i]],res.hc.z2,W=win.test)
		#print(table(duplicated(res.hc.z2)))

	}
	ppp.kids<-res.hc.z2
	print(ppp.kids$n)
	marks(ppp.kids)<-rnorm(ppp.kids$n,newkidsMean,newkidsStd)
	
	ppp.res<-superimpose(ppp.par1,ppp.kids)#res.hc.y[[i]]#res.hc.z2

	par(mfrow=c(1,3))
	plot(ppp.par1,legend=FALSE)
	plot(ppp.kids,legend=FALSE)
	plot(ppp.res,legend=FALSE)
	return(ppp.res)

}

rSpeciesRC(test1,3,1,0.5,0.1)


rSpeciesClusters<-function(ppp.par,nclust,r1)
{
	par.coords<-coords(ppp.par)

	res.hc.x<-list()

	n<-nclust#15
	radius<-r1#0.1
	for(i in 1:nrow(par.coords))
	{
		res.hc.x[[i]]<-runifdisc(n, radius, centre=c(par.coords[i,1],par.coords[i,2]))
	}

	#for(i in 1:length(res.hc.x))
	#{
	#	plot(res.hc.x[[i]],cols=3,pch=16,add=TRUE)
	#}

	#join parents with babies
	##need to subtrace empty ones?

	#i<-5
	#res.hc.y<-superimpose(res.hc.x[[i]],ppp(par.coords[i,1],par.coords[i,2],window = win.test),W=win.test)

	res.hc.y<-res.hc.x##list()

	##for(i in 1:length(res.hc.x))
	##{
	##	res.hc.y[[i]]<-superimpose(res.hc.x[[i]],ppp(par.coords[i,1],par.coords[i,2],window = win.test),W=win.test)
	##	marks(res.hc.y[[i]])<-as.factor(paste("sp",i))
###
	##}
#
	#merge and impose as different species
	res.hc.z<-list()


	for(i in 1:(length(res.hc.y)-1))
	{
		res.hc.z[[i]]<-superimpose(res.hc.y[[i]],res.hc.y[[i+1]],W=win.test)
	}
#
	#clusters with each cluster their own species

	res.hc.z2<-res.hc.z[[1]]
#
	for(i in 2:(length(res.hc.y)-1))
	{
#
		res.hc.z2<-superimpose(res.hc.z[[i]],res.hc.z2,W=win.test)
	}

	ppp.res<-res.hc.z2#res.hc.y[[i]]#res.hc.z2
	
	
	par(mfrow=c(1,3))
	plot(ppp.par,legend=FALSE)
	plot(ppp.res,legend=FALSE)
	plot(envelope(ppp.res,pcf,nsim=99,nrank=4),ylim=c(0,2),legend=FALSE)
	return(ppp.res)

}


sim.e.s1<-ppp(coords(sim.e.1)[,1],coords(sim.e.1)[,2],marks=marks((sim.e.1)),window = rect(3,4,4,5))


subsamplePPP<-function(ppp1,samplesize,maxWindSize)
{
	coordsMin<-sample(maxWindSize-1,2)
	sim.e.s1<-ppp(coords(ppp1)[,1],coords(ppp1)[,2],marks=marks(ppp1),window = owin(c(coordsMin[1],coordsMin[1]+1),c(coordsMin[2],coordsMin[2]+1)))
	


}



