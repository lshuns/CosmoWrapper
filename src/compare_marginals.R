#
#
# script to plot the marginal dz constraints compared to the prior
#
#
library(matrixStats) 

#Define the base variables 
chaindir.base<-'COSMOPIPE/GoldSet_GOLDSET/work_KV450/MCMC/output/sci_UNBLINDED/'
chainfile<-'KV450__HEADER.txt'
parameters<-c("D_z1","D_z2","D_z3","D_z4","D_z5","A_IA","Omega_m","sigma8","S8","n_s","Ac","h","ln10^{10}A_s")
parameters<-c("A_IA","n_s","h","Omega_m","sigma8","S8","ln10^{10}A_s")
param.expressions<-c("italic(A)[IA]","italic(n)[s]","h","Omega[m]","sigma[8]","italic(S)[8]","ln10^10*italic(A)[s]")
gauss.base<-"grep gaussian_prior_name CHAINDIR/log.param | awk -F= '{print $2}' | sed 's/\\[/c\\(/g' | sed 's/\\]/\\)/g' "
prior.base<-"grep \"parameters\\['PARAMETER'\\]\" CHAINDIR/log.param | awk -F= '{print $2}' | sed 's/\\[/c\\(/g' | sed 's/\\]/\\)/g' "
n.smooth=40
A_IA.prior<-function(X) dunif(X,min=-6,max=6)

#Setup the plot device
if (!interactive()) { 
  png(file=paste0("compare_marginals.png"),height=4*220,width=6*220,res=220)
} else if (dev.cur()==1) { 
  x11(height=5,width=7)
}
lay.mat<-matrix(c(1:c(length(parameters)+1),rep(0,3)),byrow=T,ncol=3)
lay.mat<-lay.mat[-nrow(lay.mat),]
lay.mat[which(lay.mat==0)]<-length(parameters)+1
layout(lay.mat)
par(mar=c(1.5,1.5,1.0,1.0),oma=c(1.0,1.5,0,0))

#Read in the command line arguments 
if (!interactive()) { 
  goldsetlist<-commandArgs(T)
} else { 
  goldsetlist<-helpRfuncs::vecsplit(readline(prompt="Enter the goldclass list: "),' ')
}
if (length(goldsetlist)==0) { 
  goldsetlist<-c("noDEEP2_shift_small_error",'noDEEP2_shift')
  goldsetlist<-c("NONE","NONE_shift","Fid_shift","noDEEP2_shift_small_error","nozCOSMOS","nozCOSMOS_shift",'noVVDS_shift','speczquality4','multispec3','noDEEP2_shift')
  goldsetlist<-c("NONE","NONE_shift","Fid_shift","noDEEP2_shift","nozCOSMOS_shift",'noVVDS_shift','speczquality4','multispec3')
}
fullsetlist<-c("NONE","NONE_shift","Fid_shift","noDEEP2_shift_small_error","nozCOSMOS","nozCOSMOS_shift",'noVVDS_shift','speczquality4','multispec3',"noDEEP2_shift",'noDEEP2','Fid')
fullsetlabs<-c("KV450-DIR","KV450-DIR dz","SOM Fiducial dz","SOM noDEEP2 dz","SOM nozCOSMOS","SOM nozCOSMOS dz",'SOM noVVDS dz','SOM speczquality4','SOM multispec3',"SOM noDEEP2 dz",'SOM Fiducial')
fullsetordr<-c(1,2,10,9,3,4,5,6,7,8,8)
collist<-c("darkblue","darkred","orange","darkgreen",RColorBrewer::brewer.pal(8,"Set2"))
fullcolordr<-c(1,2,3,9,9,5,6,8,9,7,9)

gauss.priors<-colours<-prior<-chain<-list()

goldsetlabs<-goldsetordr<-NULL
for (set in goldsetlist) { 
  set.ind<-which(fullsetlist==set)
  goldsetordr<-c(goldsetordr,fullsetordr[set.ind])
  goldsetlabs<-c(goldsetlabs,fullsetlabs[set.ind])
  colours[[set]]<-collist[fullcolordr[set.ind]]
}
goldsetordr<-order(goldsetordr)

cat('\n param ')
for (set in goldsetlist[goldsetordr]) { 
  cat(set," ")
  #Read in the chains
  chaindir<-gsub("GOLDSET",set,chaindir.base)
  if (!file.exists(paste0(chaindir,chainfile))) { 
    next
  } 
  chain[[set]]<-data.table::fread(paste0(chaindir,chainfile),fill=TRUE)
  #Fix the chain header if needed
  if (colnames(chain[[set]])[1]=="#") { 
    colnames(chain[[set]])<-c(colnames(chain[[set]])[-1],"#")
  }
  #Check for wayward commas...
  colnames(chain[[set]])<-gsub(",","",colnames(chain[[set]]))
  
  gauss.priors[[set]]<-eval(parse(text=system(helpRfuncs::vgsub(c("GOLDSET","CHAINDIR"),
                                                                c(set,chaindir),
                                                                gauss.base),intern=T)))
  #Select the colour 
  for (i in 1:length(parameters)) { 
    param<-parameters[i]
    #Read in the prior
    prior.vec<-system(helpRfuncs::vgsub(c("GOLDSET","CHAINDIR","PARAMETER"),
                                        c(set,chaindir,param),prior.base),intern=T)
    #Remove any "None" values 
    prior.vec<-gsub("None","NA",prior.vec)
    #Assign the prior vector
    prior[[set]][[param]]$base<-eval(parse(text=prior.vec))

    #Get the min and max values 
    if (length(prior[[set]][[param]]$base)!=0) { 
      prior[[set]][[param]]$xrange<-as.numeric(prior[[set]][[param]]$base[2:3])
    } else { 
      prior[[set]][[param]]$xrange<-c(NA,NA)
    } 
    #Does this parameter have a prior? 
    if (all(is.finite(prior[[set]][[param]]$xrange))) {
      #Define the prior x-values 
      step<-abs(diff(prior[[set]][[param]]$xrange))/n.smooth/2
      prior[[set]][[param]]$eval<-data.frame(x=seq(min(prior[[set]][[param]]$xrange)-5*step,max(prior[[set]][[param]]$xrange)+5*step,by=step))
      #Is the prior gaussian?
      prior[[set]][[param]]$gauss<-param%in%gauss.priors[[set]]
      if (prior[[set]][[param]]$gauss) { 
        prior[[set]][[param]]$eval$y<-dnorm(prior[[set]][[param]]$eval$x,mean=as.numeric(prior[[set]][[param]]$base[1]),
                                            sd=as.numeric(prior[[set]][[param]]$base[4]))
      } else { 
        prior[[set]][[param]]$eval$y<-dunif(prior[[set]][[param]]$eval$x,min=min(prior[[set]][[param]]$xrange),
                                            max=max(prior[[set]][[param]]$xrange))
      }
    } else { 
      prior[[set]][[param]]$eval<-data.frame(x=seq(0,1),y=c(NA,NA))
    }
  }
}

#Loop over the parameters 
for (i in 1:length(parameters)) { 
  param<-parameters[i]
  cat("\n",param," ")
  param.expr<-param.expressions[i]
  xmin<-Inf
  xmax<--Inf
  ymax<-0
  for (set in goldsetlist) { 
    if (length(chain[[set]])==0) { 
      next
    }
    if (length(chain[[set]][[param]])==0) { 
      next
    }
    if (all(is.finite(prior[[set]][[param]]$xrange))) { 
      xmax<-max(xmax,max(prior[[set]][[param]]$xrange))
      xmin<-min(xmin,min(prior[[set]][[param]]$xrange))
      ymax<-max(ymax,prior[[set]][[param]]$eval$y)
    } 
    xmax<-max(xmax,Hmisc::wtd.quantile(chain[[set]][[param]],weight=chain[[set]]$weight,normwt=TRUE,probs=0.999))
    xmin<-min(xmin,Hmisc::wtd.quantile(chain[[set]][[param]],weight=chain[[set]]$weight,normwt=TRUE,probs=0.001))
    chaindens<-density(chain[[set]][[param]],weight=chain[[set]]$weight,n=100)
    labloc<-ifelse(which.max(chaindens$y)>length(chaindens$y)/2,'topleft','topright')
    ymax<-max(ymax,max(chaindens$y))
  }
  ymax<-ymax*1.4
  step<-(xmax-xmin)/n.smooth

  #setup the plot 
  magicaxis::magplot(NA,side=2:4,labels=c(T,F,F),majorn=3,
                     xlab='',ylab=ifelse(which(lay.mat==i,arr.ind=T)[,2]==1,'PDF',''),
                     mtline=1.5,ylim=c(0,ymax),xlim=c(xmin-5*step,xmax+5*step))
  magicaxis::magaxis(side=1,labels=c(T),majorn=5)
  #Add the zero line
  abline(v=0,col='grey',lty=3)

  astro::label(labloc,lab=eval(parse(text=paste0('expression(',param.expr,')'))),cex=1.5)
  #overlay the prior
  lines(prior[["Fid_shift"]][[param]]$eval,lwd=1.2,lty=3)
  for (set in goldsetlist[goldsetordr]) { 
    if (length(chain[[set]])==0) { 
      cat("NA±NA ")
      next
    }
    if (any(!is.finite(chain[[set]]$weight))) { 
      chain[[set]]$weight[which(!is.finite(chain[[set]]$weight))]<-0
    }
    cat(paste0(round(digits=3,weighted.mean(chain[[set]][[param]],chain[[set]]$weight)),'±',
               round(digits=3,weightedSd(chain[[set]][[param]],chain[[set]]$weight/min(chain[[set]]$weight))),' '))
    lines(density(chain[[set]][[param]],weight=chain[[set]]$weight,n=1e3,
                    from=xmin,to=xmax,bw=step/sqrt(12),kernel='rect'),
          lwd=1.5,col=seqinr::col2alpha(colours[[set]],0.8))
  }
}
cat("\n")
  
#Draw the legend 
legend.col<-NULL
for (set in goldsetlist) 
  legend.col<-c(legend.col,colours[[set]])
plot(1,type='n',axes=F)
legend('left',legend=c("Prior","MarginalPosterior",goldsetlabs),lwd=2,
       col=c("black","black",legend.col),lty=c(3,1,rep(1,length(goldsetlist))),
       bty='n',pch=NA,ncol=2,cex=1.2)

#close the plot device 
if (!interactive()) {
  dev.off()
}

