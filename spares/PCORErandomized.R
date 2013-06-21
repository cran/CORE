PCORErandomized <-
function(procno=1,COREobj,boundaries,nprocs=1,rngoffset=0){
need2know<-c("seedme","nshuffle","pow","shufflemethod","input","coreTable",  
	"minscore")
for(item in need2know)assign(item,COREobj[[item]])
minscore<-max(minscore,min(coreTable[,"score"]))
#rm(COREobj)
set.seed(seedme)
myshuffles<-nshuffle%/%nprocs
shuffleres<-nshuffle%%nprocs
shuffleskip<-myshuffles*(procno-1)+min(shuffleres,procno-1)
if(procno<=shuffleres)myshuffles<-myshuffles+1
weightList<-vector(mode="list",length=nrow(coreTable))
weight<-input[,"weight"]
for(i in 1:nrow(coreTable)){
	za<-input[input[,"chrom"]==coreTable[i,"chrom"],]
	cigt<-za[,"start"]<=coreTable[i,"start"]&za[,"end"]>=coreTable[i,"end"]
	weight[input[,"chrom"]==coreTable[i,"chrom"]][cigt]<-0      
	weightList[[i]]<-matrix(ncol=2,
		data=c(which(input[,"chrom"]==coreTable[i,"chrom"])[cigt],
		weight[input[,"chrom"]==coreTable[i,"chrom"]][cigt]))
}
simscores<-matrix(ncol=nshuffle,nrow=nrow(coreTable))
chrmax<-nrow(boundaries)
advanceRNG(randopt=shufflemethod,nrand=shuffleskip+rngoffset,
	nevents=nrow(input))
for(shuffle in 1:myshuffles){
if(shufflemethod=="SIMPLE")z<-
	cbind(randomEventMoves(input[,"end"]-input[,"start"]+1,boundaries),
	input[,"weight"])
if(shufflemethod=="RESCALE")z<-
	cbind(randomRescaledEventMoves(input[,c("start","end","chrom","weight")],
	boundaries),input[,"weight"])
dimnames(z)[[2]]<-c("start","end","chrom","weight")
chu<-unique(z[,"chrom"])
chc<-rep(0,length(chu))
for(i in 1:length(chu))chc[i]<-sum(z[z[, "chrom"]==chu[i],"weight"])
scoremat<-matrix(nrow=nrow(coreTable),ncol=length(chu),data=0)
for(ich in 1:length(chu)){
	wza<-which(z[,"chrom"]==chu[rev(order(chc))][ich])
	za<-z[wza,c("start","end","weight"),drop=F]
	zaf<-za[za[,"weight"]>0,,drop=F]
	y<-cbind(c(zaf[,"end"]+1,zaf[,"start"]),c(-zaf[,"weight"],zaf[,"weight"]))
	y<-y[order(y[,1]),,drop=F]
	cy2<-cumsum(y[,2])
	mscore<-max(cy2)
	for(i in 1:nrow(coreTable)){
		if(mscore<minscore){
			scoremat[i:nrow(coreTable),ich]<-0
			break
		}	
		scoremat[i,ich]<-mscore
		if(i==nrow(coreTable))break
		whereIwas<-which(wza%in%weightList[[i]][,1])
		if(length(whereIwas)>0){
			za[whereIwas,"weight"]<-0
			zaf<-za[za[,"weight"]>0,,drop=F]
			y<-cbind(c(zaf[,"end"]+1,zaf[,"start"]),c(-zaf[,"weight"],zaf[,"weight"]))
			y<-y[order(y[,1]),,drop=F]
			cy2<-cumsum(y[,2])
			mscore<-max(cy2)
		}
	}
}
simscores[,shuffle]<-apply(scoremat,1,max)
}
return(simscores)
}
