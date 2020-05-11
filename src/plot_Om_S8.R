#
#
# Function plots Om vs S8 from MCMC Chains 
#
#

library(magicaxis)

#Load in the bespoke magcon function 
source("magcon.R") 

png(file='plot_Om_S8_summary.png',height=3*220,width=4*220,res=220)
par(mar=c(3.5,3.0,0.5,0.5))

RdBu<-colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))
xmin<-0.09
xmax<-0.52
ymin<-0.59
ymax<-0.96
buff<-c(0.05,0.05)
alpha<-c(0.8,0.5)
text.cex<-0.8

chains<-paste0("ChainDir/",c("GoldSet_NONE.txt","GoldSet_Fiducial.txt","GoldSet_NONE_shift.txt","GoldSet_noDEEP2.txt","base_plikHM_TTTEEE_lowl_lowE.txt"))
labels<-c("KV450-DIR",expression(paste("SOM-Gold-Fiducial ",delta*italic(z))),expression(paste("KV450-DIR ",delta*italic(z))),expression(paste("SOM-Gold-noDEEP2 ",delta*italic(z))),"Planck-Legacy")
plot.order<-c(1,4,2,3,5)
legend.order<-c(1,4,3,2,5)
doim<-c(F,F,F,F,F)
plty<-rbind(c(3,2),c(1,1),c(1,1),c(1,1),c(1,1))
fill<-c(F,T,T,F,T)
joint.smoothing<-c(T,T,T,T,F)

imset<-list('#000000', #KVD
            '#ff7f00', #FID
            '#984ea3', #KVDS
            '#4daf4a', #NODEEP
            '#e41a1c') #PL

if (!exists('h')) { 
  htext<-'Optimal Smoothing'
} else { 
  htext<-'Input Smoothing'
}

conlist<-list()
hlist<-matrix(NA,nrow=length(plot.order),ncol=2)
for (i in order(plot.order)) { 
  #Read the catalogue
  cat<-try(data.table::fread(chains[i])) 
  if (class(cat)[1]=='try-error') next
  #Fix the headers 
  if (colnames(cat)[1]=="V1") { 
    names<-system(paste("awk '{print $1}' ",gsub('.txt','.paramnames',chains[i])),intern=TRUE)
    colnames(cat)<-c("weights","chi2",names)
  }
  colnames(cat)<-gsub(',','',colnames(cat),fixed=T)
  colnames(cat)<-gsub('*','',colnames(cat),fixed=T)
  if (length(cat$Omega_m)==0) { 
    if (length(cat$omegam)!=0) {
      cat$Omega_m<-cat$omegam
    } else if (length(cat$omegamh2)!=0) { 
      cat$Omega_m<-cat$omegamh2
    } else { 
      stop("No omegamh2 or omegam")
    }
  } else { 
    colnames(cat)<-c(colnames(cat)[-1],colnames(cat)[1])
  }
  if (length(cat$S8)==0) { 
    cat$S8<-cat$sigma8*sqrt(cat$Omega_m/0.3)
  }
  if (length(cat$weights)==0) { 
    cat$weights<-rep(1,nrow(cat))
  }
  #Define the image colours 
  imcols<-seqinr::col2alpha(imset[[i]],alpha=seq(0,alpha[1]-min(0.2*i,0.5),length=1e3))
  #Plot the contours
  if (i==plot.order[1]) { 
    #Start the plot from scratch
    plot(NA,type='n',xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab="",ylab="",axes=FALSE)
    S8vec<-try(invisible(helpRfuncs::vecsplit(system(paste0("awk '/^S8/{i++}i==2{print $0; exit}' GoldSet_Fid_shift/work_KV450/MCMC/output/sci_UNBLINDED/parameter_table.txt"),intern=T),' ')))
    S8vec<-S8vec[which(S8vec!='')]
    S8vec<-gsub(',','',S8vec)
    S8<-as.numeric(S8vec[2])
    S8lo<-as.numeric(S8vec[3])
    S8hi<-as.numeric(S8vec[4])
    rect(yb=S8+S8lo,xl=0,yt=S8+S8hi,xr=1.,col=seqinr::col2alpha('grey',0.8),border=F)
    if (exists('h')) { 
      res=magcon(cat$Omega_m,cat$S8,weights=cat$weights,conlevels=c(diff(pnorm(c(-1,1))),diff(pnorm(c(-2,2)))),xlim=c(xmin,xmax),ylim=c(ymin,ymax),
            nbins=0,ngrid=1000,fill.col=seqinr::col2alpha(imset[i],c(alpha[1],alpha[2])),imcol=imcols,doim=doim[i],col=NA,fill=fill[i],lwd=2,
            add=TRUE,barposition='bottomright',barorient='h',h=h,dobar=F)
    } else { 
      res=magcon(cat$Omega_m,cat$S8,weights=cat$weights,conlevels=c(diff(pnorm(c(-1,1))),diff(pnorm(c(-2,2)))),xlim=c(xmin,xmax),ylim=c(ymin,ymax),
            nbins=0,ngrid=1000,fill.col=seqinr::col2alpha(imset[i],c(alpha[1],alpha[2])),imcol=imcols,doim=doim[i],col=NA,fill=fill[i],lwd=2,
            add=TRUE,barposition='bottomright',barorient='h',dobar=F)
    }
    hlist[i,]=res$h
  } else { 
    if (joint.smoothing[i]) { 
      h=hlist[plot.order[1],]
    }
    #Add to existing plot 
    if (exists('h')) { 
      res=magcon(cat$Omega_m,cat$S8,weights=cat$weights,conlevels=c(diff(pnorm(c(-1,1))),diff(pnorm(c(-2,2)))),
             nbins=0,ngrid=1000,fill.col=seqinr::col2alpha(imset[i],c(alpha[1],alpha[2])),imcol=imcols,doim=doim[i],add=T,col=NA,fill=fill[i],h=h,lwd=2,dobar=F)
    } else { 
      res=magcon(cat$Omega_m,cat$S8,weights=cat$weights,conlevels=c(diff(pnorm(c(-1,1))),diff(pnorm(c(-2,2)))),
             nbins=0,ngrid=1000,fill.col=seqinr::col2alpha(imset[i],c(alpha[1],alpha[2])),imcol=imcols,doim=doim[i],add=T,col=NA,fill=fill[i],lwd=2,dobar=F)
    } 
    if (joint.smoothing[i]) { 
      rm('h')
    }
    hlist[i,]=res$h
  }
  conlist[[i]]<-res$contours
}
#Draw the contour boundaries 
for (i in order(plot.order,decreasing=T)) { 
  if (length(conlist)<i || length(conlist[[i]])==0) next
  for (j in 1:length(conlist[[i]])) { 
    if (!grepl("#",imset[i],fixed=T)) { 
      lines(conlist[[i]][[j]],col=seqinr::col2alpha(paste0('dark',imset[i]),0.4),lwd=1,lty=plty[i,j])
    } else { 
      lines(conlist[[i]][[j]],col=seqinr::col2alpha(imset[[i]],0.4),lwd=1,lty=plty[i,j])
    }
  }
}
#Draw the colourbar
if (any(doim)) { 
  res=magcon(cat$Omega_m,cat$S8,weights=cat$weights,conlevels=c(diff(pnorm(c(-1,1))),diff(pnorm(c(-2,2)))),xlim=c(xmin,xmax),ylim=c(ymin,ymax),add=TRUE,
             nbins=0,ngrid=1000,imcol=seqinr::col2alpha("black",alpha=seq(1,0,len=1e3)),doim=F,docon=F,dobar=T,lwd=2,barposition='bottomright',barorient='h')
}

h<-matrixStats::colMaxs(hlist)
#Define the x and y buffers 
xbuff<-abs(diff(c(xmin,xmax)))*buff[1]
ybuff<-abs(diff(c(ymin,ymax)))*buff[2]
#Draw the smoothing Kernel
text(xmax-xbuff-h[1]*3,ymax-ybuff,lab=htext,cex=text.cex-0.2)
text(xmax-xbuff-h[1]*3,ymax-ybuff-max(h[2]*7,0.02),lab='Kernel',cex=text.cex-0.2)
#Draw the smoothing kernel 
for (i in which(!joint.smoothing)) { 
  htmp<-hlist[i,]
  if (!is.null(htmp) && !all(is.na(htmp))) {
    magcon(rep(xmax-xbuff-h[1]*3,1e4),rep(ymax-ybuff-max(h[2]*3.5,0.01),1e4),conlevels=c(diff(pnorm(c(-1,1))),diff(pnorm(c(-2,2)))),xlim=c(0.09,0.52),ylim=c(0.59,0.94),
          nbins=0,ngrid=1000,imcol=c(NA,rev(RdBu(1e3))),doim=F,add=T,h=htmp,lwd=1.5,dobar=F,col=imset[[i]])
  }
}
htmp<-hlist[which(joint.smoothing)[1],]
if (!is.null(htmp) && !all(is.na(htmp))) {
  magcon(rep(xmax-xbuff-h[1]*3,1e4),rep(ymax-ybuff-max(h[2]*3.5,0.01),1e4),conlevels=c(diff(pnorm(c(-1,1))),diff(pnorm(c(-2,2)))),xlim=c(0.09,0.52),ylim=c(0.59,0.94),
        nbins=0,ngrid=1000,imcol=c(NA,rev(RdBu(1e3))),doim=F,add=T,h=htmp,lwd=1.5,dobar=F,col='black')
}
text(rep(xmin+xbuff/2,length(legend.order)),seq(ymax-ybuff,0.85,len=length(legend.order)),
     labels=labels[order(legend.order)],col=unlist(imset)[order(legend.order)],pos=4,cex=text.cex)

magaxis(xlab=expression(Omega[m]),side=1,labels=T,mtline=1.8,lab.cex=1.5)
magaxis(ylab=expression(italic(S)[8]*"="*sigma[8]*sqrt(Omega[m]/0.3)),side=2:4,labels=c(T,F,F),mtline=1.5)

dev.off()
