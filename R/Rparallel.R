Rparallel<-function(randfun,distrib,doshuffles,nshuffle,dataIn,returnme,
	boundaries,njobs){
        if(distrib=="Rparallel"&doshuffles!="NO"){
                #ncores<-min(njobs,switch(doshuffles=="FROMSCRATCH",
                #        nshuffle,nshuffle-dataIn$nshuffle),parallel::detectCores())
                #cl<-parallel::makeCluster(getOption("cl.cores",ncores))
                ncores<-min(njobs,ifelse(doshuffles=="FROMSCRATCH",
												 nshuffle,nshuffle-dataIn$nshuffle))
								cl<-parallel::makeCluster(ncores)
								parallel::clusterEvalQ(cl=cl,expr=library(CORE))
        }
        if(substring(distrib,1,3)=="mpi"&doshuffles!="NO"){
								WrRmpi()
								ncores<-min(njobs,ifelse(doshuffles=="FROMSCRATCH", 
                         nshuffle,nshuffle-dataIn$nshuffle))
								nshuffle<-returnme$nshuffle
								if(doshuffles=="FROMSCRATCH"){
												returnme$nshuffle<-returnme$nshuffle%/%ncores
												if(nshuffle%%ncores!=0)returnme$nshuffle<-returnme$nshuffle+1
								}
								if(doshuffles=="ADD"){
												returnme$nshuffle<-(returnme$nshuffle-dataIn$nshuffle)%/%ncores
												if((nshuffle-dataIn$nshuffle)%%ncores!=0)returnme$nshuffle<-returnme$nshuffle+1
								}
								save(ncores,doshuffles,dataIn,randfun,boundaries,
												returnme,file=paste(getwd(),"/tempMPI",sep=""))
								if(distrib=="mpi.ge"){
												sink(paste(getwd(),"/Rmpi.sh",sep=""))
												cat("mpirun -np $NSLOTS Rscript $PWD/Rmpi.R")
												sink()	
												system(paste("qsub -V -cwd -sync y -l m_mem_free=2G -pe mpi",as.character(ncores),"Rmpi.sh"))
								}
								load(paste(getwd(),"/mygather.temp",sep=""))
								gather.result<-get("gather.result")
								if(doshuffles=="FROMSCRATCH")mpidata<- gather.result[1:(nshuffle*nrow(returnme$coreTable))]
								if(doshuffles=="ADD")mpidata<- gather.result[1:((nshuffle-dataIn$nshuffle)*nrow(returnme$coreTable))]
        }
        if(doshuffles=="FROMSCRATCH"){
                returnme$simscores<-switch(distrib,
                        vanilla=randfun(COREobj=returnme,boundaries=boundaries),
                        Rparallel=matrix(nrow=nrow(returnme$coreTable),
                                data=unlist(parallel::parSapply(cl=cl,X=1:ncores,FUN=randfun,
                                COREobj=returnme,boundaries=boundaries,nprocs=ncores))),
												mpi.ge=matrix(nrow=nrow(returnme$coreTable),data=mpidata)
                        )
        }
        else if(doshuffles=="ADD"){
                returnme$simscores<-
                cbind(dataIn$simscores[1:nrow(returnme$coreTable),,drop=F],switch(distrib,
                        vanilla=randfun(COREobj=returnme,boundaries=boundaries,
                                rngoffset=dataIn$nshuffle),
                        Rparallel=matrix(nrow=nrow(returnme$coreTable),
                                data=unlist(parallel::parSapply(cl=cl,X=1:ncores,
                                FUN=randfun,COREobj=returnme,boundaries=boundaries,
                                rngoffset=dataIn$nshuffle,nprocs=ncores))),
												mpi.ge=matrix(nrow=nrow(returnme$coreTable),data=mpidata)
                ))
        }
        if("simscores"%in%names(returnme))returnme$p<-
                (rowSums(returnme$simscores>returnme$coreTable[,"score"])+1)/
                (ncol(returnme$simscores)+2)
        if(exists("cl"))stopCluster(cl)
				if(substring(distrib,1,3)=="mpi"){
								returnme$nshuffle<-nshuffle
								system("rm mygather.temp")
        				system("rm tempMPI")
								system("rm Rmpi.*")
				}
	return(returnme)
}
