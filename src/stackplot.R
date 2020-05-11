#
#
# Creates the S8 Stack plot 
#
#

library(magicaxis)
root<-"./GoldSet_"
root<-"./COSMOPIPE/GoldSet_"
file<-"/work_KV450/MCMC/output/sci_UNBLINDED/parameter_table.txt"
command<-"'/^S8/{i++}i==2{print $0; exit}' "   

png(file='stackplot.png',height=3*220,width=3*220,res=220)
col=c('#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628', '#f781bf', '#FFD700')
goldsets<-c("NONE","NONE_shift","Fid_shift","noDEEP2_shift","nozCOSMOS_shift",'noVVDS_shift','speczquality4','multispec3','noDEEP2','Fid')
goldlabs<-c("KV450-DIR","KV450-DIR","SOM Fiducial","SOM noDEEP2","SOM nozCOSMOS",'SOM noVVDS','SOM speczquality4','SOM multispec3',"SOM noDEEP2","SOM Fiducial")
goldcols<-c("black",col[4],col[5],col[3],'black','black','black','black',col[3],'black')
goldyvals<-c(9,4,0,1,2,3,7,8,21,6)
goldpch=c(1,1,23,1,1,1,1,1,1,1,1)
breaklines<-c(5,10)
ymin<-14
ymax<--1

par(mar=c(3.2,0.5,0.5,0.5))
layout(1) 
plot(axes=F,xlim=c(0.63,1.075),ylim=c(ymin,ymax),NA,xlab='',ylab='',type='n')
text(0.65,breaklines[1]/2,label=expression(delta*italic(z)!=0),cex=0.8,srt=90)
if (breaklines[1]!=breaklines[2]) { 
  text(0.65,(breaklines[2]+breaklines[1])/2,label=expression(delta*italic(z)==0),cex=0.8,srt=90)
}
text(0.65,(breaklines[2]+ymin)/2,label="Ext",cex=0.8,srt=90)

#Make the grey shading and separation lines
S8vec<-try(invisible(helpRfuncs::vecsplit(system(paste0('awk ',command,root,"Fid_shift",file),intern=T),' ')))
S8vec<-S8vec[which(S8vec!='')]
S8vec<-gsub(',','',S8vec)
S8<-as.numeric(S8vec[2])
S8lo<-as.numeric(S8vec[3])
S8hi<-as.numeric(S8vec[4])
rect(xl=S8+S8lo,yb=20,xr=S8+S8hi,yt=-5,col=seqinr::col2alpha('grey',0.7),border=NA)
abline(v=S8,col=seqinr::col2alpha('darkgrey',0.8),lwd=2)
#Draw the axes
magaxis(side=c(1,3),labels=c(T,F),xlab=expression(italic(S[8])*"="*sigma[8]*sqrt(Omega[m]/0.3)),majorn=5,minorn=5,frame.plot=T)

#Planck
planck.S8<-0.812
planck.S8lo<-0.802-planck.S8
planck.S8hi<-0.822-planck.S8
planck.y<-ymin-1
points(pch=1,col=col[1],planck.S8,planck.y,lwd=2,cex=1.2)
magerr(pch=1,col=col[1],planck.S8,planck.y,xlo=planck.S8lo,xhi=planck.S8hi,lwd=2)
text(0.95,planck.y,label="Planck-Legacy",cex=0.8)
#Hildebrandt2019
h19.S8<-0.737
h19.S8lo<-0.036
h19.S8hi<-0.040
h19.y<-ymin-3
#Draw the rectangle around the H2020 and KV450DIR points
rect(xl=h19.S8-h19.S8lo-0.0075,xr=h19.S8+h19.S8hi+0.0075,yb=h19.y-2.5,yt=h19.y+0.5,col='darkblue',lwd=1,lty=2,density=0,border=T)
points(pch=1,col='black',h19.S8,h19.y,lwd=2,cex=1.2)
magerr(pch=1,col='black',h19.S8,h19.y,xlo=h19.S8lo,xhi=h19.S8hi,lwd=2)
text(0.95,h19.y,label="Hildebrandt+ (2020)",cex=0.8)

#Hildebrandt2019 noDEEP
h19.S8<-0.761
h19.S8lo<-0.037
h19.S8hi<-0.041
h19.y<-ymin-2
points(pch=1,col='black',h19.S8,h19.y,lwd=2,cex=1.2)
magerr(pch=1,col='black',h19.S8,h19.y,xlo=h19.S8lo,xhi=h19.S8hi,lwd=2)
text(0.95,h19.y,label="H+20 noDEEP2",cex=0.8)

#Draw the lines 
abline(h=breaklines,lty=2,lwd=2,col='black')

#Gold Samples
for (i in 1:length(goldsets)) { 
  samp<-goldsets[i] 
  yval<-goldyvals[i] 
  #Read the chain
  S8vec<-try(invisible(helpRfuncs::vecsplit(system(paste0('awk ',command,root,samp,file),intern=T),' ')))
  if (length(S8vec)!=0) { 
    S8vec<-S8vec[which(S8vec!='')]
    S8vec<-gsub(',','',S8vec)
    S8<-as.numeric(S8vec[2])
    S8lo<-as.numeric(S8vec[3])
    S8hi<-as.numeric(S8vec[4])
    if (goldpch[i]==23) { 
      points(pch=goldpch[i],col=goldcols[i],x=S8,y=yval,lwd=1,bg=goldcols[i],cex=1.5)
    } else { 
      points(pch=goldpch[i],col=goldcols[i],x=S8,y=yval,lwd=2,bg=goldcols[i],cex=1.2)
    }
    magerr(pch=goldpch[i],col=goldcols[i],x=S8,y=yval,xlo=S8lo,xhi=S8hi,lwd=2)
    cat(paste0(goldlabs[i],': ',S8,'+',S8hi,S8lo,'\n'))
  }
  text(0.95,yval,label=goldlabs[i],cex=0.8)
}
dev.off()
