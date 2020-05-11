# Script to construct specz samples from list driven simulations
#
#
binlim<-c(0.1,0.3,0.5,0.7,0.9,1.2)

inputs<-commandArgs(T)
#Read in catalogue
cat<-helpRfuncs::read.file(inputs[1],data.table=F,extname='OBJECTS')
print(colnames(cat))
for (i in 1:5) {
  #Write the tomo catalogue:
  tomo<-cat[which(cat$Z_B > binlim[i] & cat$Z_B <= binlim[i+1]),]
  outend<-helpRfuncs::vecsplit(inputs[1],'.',-1)
  outfile<-sub(paste0(".",outend),paste0("_TOMO",i,'.',outend),inputs[1])
  if (outfile==inputs[1]) { stop("unable to create a distinct outfile because infile has no extension?!") }
  if (!is.na(inputs[2])) {
    outfile<-helpRfuncs::vecsplit(outfile,'/',-1)
    outfile=file.path(inputs[2],outfile)
  }
  print(outfile)
  helpRfuncs::write.file(file=outfile,tomo)
}
