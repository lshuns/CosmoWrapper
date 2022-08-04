#
# Script to construct survey goldclass. 
# Adapted from script construct_survey_goldclass.R & others 
#

inputs<-commandArgs(T)

library(kohonen)
library(plotrix)
library(helpRfuncs)
library(data.table)

#Function for constructing new output paths /*fold*/{{{
outputpath<-function(newstring,catalogue,output.folder,format) { 
  if (!dir.exists(output.folder)) { 
    dir.create(output.folder,recursive=TRUE)
  }
  if (grepl('\\.csv',catalogue,ignore.case=TRUE)) { 
    extension<-'.csv'
  } else if (grepl('\\.cat',catalogue,ignore.case=TRUE)){
    warning("Cannot write to LDAC; preparing to write as FITS instead!") 
    extension<-'.fits'
    catalogue<-sub(".cat",".fits",catalogue)
  } else if (grepl('\\.cat',catalogue,ignore.case=TRUE)) { 
    extension<-'.cat'
  } else if (grepl('\\.asc',catalogue,ignore.case=TRUE)) { 
    extension<-'.asc'
  } else if (grepl('\\.fits',catalogue,ignore.case=TRUE)) { 
    extension<-'.fits'
  } else if (grepl('\\.Rdata',catalogue,ignore.case=TRUE)) { 
    extension<-'.Rdata'
  } else { 
    stop(paste0("Unknown File Extension:\n",catalogue))
  }
  if (!missing(format)) { 
    outname<-gsub(extension,paste0(newstring,format),catalogue)
  } else {
    outname<-gsub(extension,paste0(newstring,extension),catalogue)
  }
  outname<-vecsplit(outname,'/',-1)
  if (outname==catalogue) { 
    stop(paste0("Cannot create unique filename for output catalogue!!\n",catalogue))
  }
  outname<-file.path(output.folder,outname)
  return=outname
}
#/*fend*/}}}

#Set the default parameters /*fold*/ {{{
force<-FALSE
do.surveyClass<-FALSE
do.redshiftClass<-TRUE
blind.list<-c("A","B","C")
count.variable<-'recal_weight'
output.folder<-"DR4_SOMweight_Results/"
spec.variable<-'redshift'
QC.expr<-"abs(mean(spec$redshift,na.rm=T)-weighted.mean(phot$Z_B,phot$recal_weight,na.rm=T))<=0.61" 
tomo.bin<-c(4000,2200,2800,4200,2000)
tomo.lim<-c(0.1,0.3,0.5,0.7,0.9,1.2)
nz.format<-'.fits'
#/*fend*/}}}

##Define the various SOM GoldClass sets /*fold*/{{{
#class.sets<-rbind(c("Fid",          "phot$SurveyGoldFlag>=0",         "spec$z_Flag>0"      ),
#                  c("noVVDS",       "!phot$SurveyGoldFlag%in%c(0,16)","spec$SurveyFlag!=16"),
#                  c("nozCOSMOS",    "!phot$SurveyGoldFlag%in%c(0,2)", "spec$SurveyFlag!=2" ),
#                  c("noDEEP2",      "!phot$SurveyGoldFlag%in%c(0,4)", "spec$SurveyFlag!=4" ),
#                  c("speczquality4","phot$RedshiftGoldClass>=2",      "spec$z_Flag>=4"     ),
#                  c("multispec3",   "phot$SurveyGoldClass>=3",        "spec$z_Flag>0"      ))
##/*fend*/}}}
#Define the JLvdB SOM GoldClass sets /*fold*/{{{
class.sets<-rbind(c("nQ4",            "phot$RedshiftGoldClass>=0",  "spec$z_Flag>=4",                     "Zbest"),
                  c("Fid",            "phot$RedshiftGoldClass>=0",  "spec$z_Flag>=3",                     "Zbest"),
                  c("plusPAUS",       "phot$RedshiftGoldClass>=0",  "spec$z_Flag>=2",                     "Zbest"),
                  c("plusPAUSCOS15",  "phot$RedshiftGoldClass>=0",  "spec$z_Flag>=1",                     "Zbest"),
                  c("onlyPAUS",       "phot$RedshiftGoldClass>=0",  "bitwAnd(spec$Zsource,2L)==2L",       "Zpaus"),
                  c("onlyCOS15",      "phot$RedshiftGoldClass>=0",  "bitwAnd(spec$Zsource,4L)==4L",       "Zcos15"),
                  c("onlyPAUSCOS15",  "phot$RedshiftGoldClass>=0",  "spec$Zsource>1",                     "ZbestPhot"),
                  c("Fid_noDEVILS",   "phot$RedshiftGoldClass>=0",  "(spec$z_Flag>=3)&(spec$source!=10)", "Zbest"))
#/*fend*/}}}

#Define the datasets /*fold*/{{{
labels<-c("CDFS","zCOSMOS","DEEP2","G15DEEP","VVDS")
dataset<-c(0,100,200,400,500,1000)
#/*fend*/}}}

#Read the options /*fold*/{{{
while (length(inputs)!=0) {
  #Check the options syntax /*fold*/{{{
  while (length(inputs)!=0 && inputs[1]=='') { inputs<-inputs[-1] }  
  if (!grepl('^-',inputs[1])) {
    print(inputs)
    stop(paste("Incorrect options provided!",
               "Check the lengths for each option!\n"))
  }
  #Check for test variable /*fold*/{{{
  if (any(inputs=='--test')) { 
    inputs<-c(inputs[which(inputs=='--test')],inputs[which(inputs!='--test')])
  }
  #/*fend*/}}}
  #/*fend*/}}}
  if (inputs[1]=='-p') { 
    #/*fold*/{{{
    inputs<-inputs[-1]
    phot.catalogue<-inputs[1]
    inputs<-inputs[-1]
    #/*fend*/}}}
  } else if (inputs[1]=='-cv') { 
    #/*fold*/{{{
    warning("Modifying count variable from Fiducial?!") 
    cat("WARNING: Modifying count variable from Fiducial?!\n") 
    inputs<-inputs[-1]
    count.variable<-inputs[1]
    inputs<-inputs[-1]
    #/*fend*/}}}
  } else if (inputs[1]=='-qc') { 
    #/*fold*/{{{
    warning("Modifying QC Expression from Fiducial?!") 
    cat("WARNING: Modifying QC Expression from Fiducial?!\n") 
    inputs<-inputs[-1]
    QC.expr<-inputs[1]
    inputs<-inputs[-1]
    QC.string<-inputs[1]
    if (QC.string=='') { 
      stop("Error in QC variable specification. Syntax is -qc 'expression' 'QC_ID_string'")
    }
    inputs<-inputs[-1]
    #/*fend*/}}}
  } else if (inputs[1]=='-s') { 
    #/*fold*/{{{
    inputs<-inputs[-1]
    spec.catalogue<-inputs[1]
    inputs<-inputs[-1]
    #/*fend*/}}}
  } else if (inputs[1]=='--outputpath') { 
    #/*fold*/{{{
    inputs<-inputs[-1]
    output.folder<-inputs[1]
    inputs<-inputs[-1]
    #/*fend*/}}}
  } else if (inputs[1]=='--som') { 
    #/*fold*/{{{
    inputs<-inputs[-1]
    trained.som<-inputs[1]
    inputs<-inputs[-1]
    #/*fend*/}}}
  } else if (inputs[1]=='--tomobins') { 
    #/*fold*/{{{
    warning("Modifying tomographic cluster numbers from Fiducial?!") 
    cat("WARNING: Modifying tomographic cluster nubmers from Fiducial?!\n") 
    print(inputs)
    inputs<-inputs[-1]
    if (any(grepl('^-',inputs))) { 
      tomo.bin<-(inputs[1:(which(grepl('^-',inputs))[1]-1)])
      inputs<-inputs[-(1:(which(grepl('^-',inputs))[1]-1))]
    } else { 
      tomo.bin<-(inputs)
      inputs<-NULL
    } 
    if (length(tomo.bin)==1 && grepl(" ",tomo.bin)) { 
      tomo.bin<-vecsplit(tomo.bin," ")
    }
    tomo.bin<-as.numeric(tomo.bin)
    if (any(is.na(tomo.bin))) { 
      stop("--tomobins specified incorrectly: resulted in NA tomobin.")
    }
    #/*fend*/}}}
  } else if (inputs[1]=='--blinds') { 
    #/*fold*/{{{
    warning("Modifying Blinds away from Fiducial?!") 
    cat("WARNING: Modifying Blinds away from Fiducial?!\n") 
    inputs<-inputs[-1]
    if (any(grepl('^-',inputs))) { 
      blind.list<-inputs[1:(which(grepl('^-',inputs))[1]-1)]
      inputs<-inputs[-(1:(which(grepl('^-',inputs))[1]-1))]
    } else { 
      blind.list<-inputs
      inputs<-NULL
    } 
    #/*fend*/}}}
  } else if (inputs[1]=='--tomolims') { 
    #/*fold*/{{{
    warning("Modifying tomographic bins from Fiducial?!") 
    cat("WARNING: Modifying tomographic bins from Fiducial?!\n") 
    inputs<-inputs[-1]
    if (any(grepl('^-',inputs))) { 
      tomo.lim<-as.numeric(inputs[1:(which(grepl('^-',inputs))[1]-1)])
      inputs<-inputs[-(1:(which(grepl('^-',inputs))[1]-1))]
    } else { 
      tomo.lim<-as.numeric(inputs)
      inputs<-NULL
    } 
    #/*fend*/}}}
  } else if (inputs[1]=='--force') { 
    #/*fold*/{{{
    inputs<-inputs[-1]
    force<-TRUE
    #/*fend*/}}}
  } else if (inputs[1]=='--nzformat') { 
    #/*fold*/{{{
    inputs<-inputs[-1]
    nz.format<-inputs[1]
    inputs<-inputs[-1]
    #/*fend*/}}}
  } else {
    stop(paste("Unknown option",inputs[1]))
  }
}
#/*fend*/}}}

#Check the input tomolims /*fold*/{{{
maxBin<-length(tomo.lim)-1
if (length(tomo.bin)!=length(tomo.lim)-1) { 
  stop("mismatched tomo.lim and tomo.bin lengths!")
}
#/*fend*/}}}

#Load in the Fiducial SOM /*fold*/{{{
cat("Load SOM Data\n")
nam<-load(trained.som)
if (nam=='nam') { 
  #Whoops, we overwrote the SOM with it's name!
  load(trained.som)
  fid.som<-get('nam')
} else if (nam!='fid.som') { 
  fid.som<-get(nam)
  rm(nam)
}
#/*fend*/}}}
if (length(fid.som$whiten.param)==0) { 
  cat("WARNING: This SOM contains no whitening parameters! Assuming that it is an old format!")
  fid.som$whiten.param<-fid.som$rescale.param
}

#Read in the full spectroscopic catalogue /*fold*/{{{
if (!exists('spec')) { 
  spec<-read.file(spec.catalogue,extname='OBJECTS')
} else { 
  cat("Using 'spec' from Global Environment!\n")
}
#/*fend*/}}}

#Read in the photometric catalogue /*fold*/{{{
if (!exists("phot")) {
  cat("Read Phot Data\n")
  phot<-read.file(phot.catalogue,extname='OBJECTS')
} else { 
  cat("Using 'phot' from Global Environment!\n")
}
#/*fend*/}}}

#Flag any bad GAAP magnitudes /*fold*/{{{
print(colnames(phot))
phot.good.rows<-matrixStats::rowAlls(phot[,paste0('MAG_GAAP_',c('u','g','r','i','Z','Y','J','H','Ks'))]>0)
#/*fend*/}}}

#Parse the full spectro cat to the Fiducial SOM /*fold*/{{{
if (!exists('spec.fid.som')) { 
  cat("Parse Spec Data\n")
  spec.fid.som<-kohparse(som=fid.som,data=spec,data.missing=-99,data.threshold=c(0,40),n.cores=64)
} else { 
  cat("Using 'spec.fid.som' from Global Environment!\n")
}
spec.fid.som$hclust<-hclust(dist(spec.fid.som$codes[[1]]))
#/*fend*/}}}
#Parse the full photometric cat to the SOM /*fold*/{{{
if (nrow(phot)!=length(fid.som$unit.classif)) { 
  if (!exists("phot.fid.som")) { 
    #Parse the full Photom cat to the Fiducial SOM
    cat("Parse phot Data\n")
    phot.fid.som<-kohparse(som=fid.som,data=phot,data.missing=-99,data.threshold=c(0,40),n.cores=64)
  } else { 
  cat("Using 'phot.fid.som' from Global Environment!\n")
  }
} else { 
  cat("No need to parse phot Data\n")
  phot.fid.som<-fid.som
}
phot.fid.som$hclust<-hclust(dist(phot.fid.som$codes[[1]]))
#/*fend*/}}}

#Define the output catalogue name /*fold*/{{{
outname.phot<-outputpath("_allgoldclass",phot.catalogue,output.folder)
outname.spec<-outputpath("_allgoldclass",spec.catalogue,output.folder)
#/*fend*/}}}
if (!file.exists(outname.phot)|!file.exists(outname.spec)) {   
  #Construct the GOLD Flags /*fold*/{{{
  cat("Output File does not exist!\nCreating GoldFlag and GoldClass\n")
  for (blind in blind.list) { 
    if (blind != "NONE") { 
      blind.count.variable<-paste0(count.variable,'_',blind)
    } else { 
      blind.count.variable<-count.variable
    } 
    if (do.surveyClass) { 
      #Survey GoldClass /*fold*/{{{
      #Initialise the gold flag /*fold*/{{{
      phot[[paste0("SurveyGoldFlag_",blind)]]<-0
      phot[[paste0("SurveyGoldClass_",blind)]]<-0
      spec[[paste0("SurveyFlag_",blind)]]<-0
      #/*fend*/}}}
      #Loop through each dataset /*fold*/ {{{ 
      for (i in 1:(length(dataset)-1)) { 
        cat(paste0(labels[i],' only & blind ',blind,' & '))
        #Loop through the tomographic bins /*fold*/{{{
        for (tomo in 1:maxBin) { 
          #Define the spec sample /*fold*/{{{
          spec.index<-which(spec$Source >= dataset[i] & spec$Source < dataset[i+1] & spec$Z_B>tomo.lim[tomo] & spec$Z_B<=tomo.lim[tomo+1])
          #/*fend*/}}}
          #Define the phot sample /*fold*/{{{
          phot.index<-which(phot$Z_B>tomo.lim[tomo] & phot$Z_B<=tomo.lim[tomo+1] & phot[[blind.count.variable]]>0 & phot.good.rows)
          #/*fend*/}}}
          #Generate the SOM groupings /*fold*/{{{
          spec.tmp.som<-generate.kohgroups(spec.fid.som,n.cluster.bins=tomo.bin[tomo],n.cores=64,subset=spec.index,quiet=T)
          phot.tmp.som<-generate.kohgroups(phot.fid.som,n.cluster.bins=tomo.bin[tomo],n.cores=64,subset=phot.index,quiet=T)
          #/*fend*/}}}
          #Generate the weights /*fold*/{{{
          spec.tmp.weights<-generate.kohgroup.property(som=spec.tmp.som,
                                                       data=spec,
                                                       n.cluster.bins=tomo.bin[tomo],n.cores=64,
                                                       expression="nrow(data)",quiet=TRUE)$property$value.1
          phot.tmp.weights<-generate.kohgroup.property(som=phot.tmp.som,
                                                       data=phot,
                                                       n.cluster.bins=tomo.bin[tomo],n.cores=64,
                                                       expression="sum(data[[blind.count.variable]])",quiet=TRUE)$property$value.1
          tmp.weights<-phot.tmp.weights/spec.tmp.weights
          tmp.weights<-list(refr.weight=1/tmp.weights[phot.tmp.som$clust.classif],train.weight=tmp.weights[spec.tmp.som$clust.classif])
          tmp.weights$refr.weight[which(!is.finite(tmp.weights$refr.weight))]<-0
          tmp.weights$train.weight[which(!is.finite(tmp.weights$train.weight))]<-0
          #/*fend*/}}}
          #Delta n_eff /*fold*/{{{
          cat(" $ ")
          cat(tround(digits=1,100*
              (sum(phot[[blind.count.variable]][which(tmp.weights$refr.weight!=0)])^2/
               sum(phot[[blind.count.variable]][which(tmp.weights$refr.weight!=0)]^2))/
              (sum(phot[[blind.count.variable]][phot.index])^2/
               sum(phot[[blind.count.variable]][phot.index]^2))
              ))
          cat(" $ & ")
          #/*fend*/}}}
          #Represented by {"CDFS","zCOSMOS","DEEP2","G15DEEP","VVDS"} == {1,2,4,8,16} /*fold*/{{{
          phot[[paste0("SurveyGoldFlag_",blind)]][phot.index]<-phot[[paste0("SurveyGoldFlag_",blind)]][phot.index]+ifelse(tmp.weights$refr.weight[phot.index]!=0,1,0)*2^(i-1)
          phot[[paste0("SurveyGoldClass_",blind)]][phot.index]<-phot[[paste0("SurveyGoldClass_",blind)]][phot.index]+ifelse(tmp.weights$refr.weight[phot.index]!=0,1,0)
          spec[[paste0("SurveyFlag_",blind)]][spec.index]<-spec[[paste0("SurveyFlag_",blind)]][spec.index]+2^(i-1)
          #/*fend*/}}}
        }
        #/*fend*/}}}
        cat("\n")
      }
      #/*fend*/}}}
      #/*fend*/}}}
    }
    if (do.redshiftClass) {
      #Redshift GoldClass /*fold*/ {{{ 
      #Initialise the gold flag /*fold*/{{{
      phot[[paste0("RedshiftGoldClass_",blind)]]<-0
      spec[[paste0("RedshiftGoldClass_",blind)]]<-0
      #/*fend*/}}}
      #Loop through nQ flags /*fold*/ {{{
      for (nQ in c(1,2,3,4)) { 
        cat(paste0("nQ≥",nQ,' only & blind ',blind,' & '))
        #Loop through the tomographic bins /*fold*/{{{
        for (tomo in 1:maxBin) { 
          #Define the spec sample /*fold*/{{{
          spec.index<-which(spec$z_Flag >= nQ & spec$Z_B>tomo.lim[tomo] & spec$Z_B<=tomo.lim[tomo+1])
          #/*fend*/}}}
          #Define the phot sample /*fold*/{{{
          phot.index<-which(phot$Z_B>tomo.lim[tomo] & phot$Z_B<=tomo.lim[tomo+1] & phot[[blind.count.variable]]>0 & phot.good.rows)
          #/*fend*/}}}
          #Generate the SOM groupings /*fold*/{{{
          spec.tmp.som<-generate.kohgroups(spec.fid.som,n.cluster.bins=tomo.bin[tomo],n.cores=64,subset=spec.index,quiet=TRUE)
          phot.tmp.som<-generate.kohgroups(phot.fid.som,n.cluster.bins=tomo.bin[tomo],n.cores=64,subset=phot.index,quiet=TRUE)
          #/*fend*/}}}
          #Generate the weights /*fold*/{{{
          spec.tmp.weights<-generate.kohgroup.property(som=spec.tmp.som,
                                                       data=spec,
                                                       n.cluster.bins=tomo.bin[tomo],n.cores=64,
                                                       expression="nrow(data)",quiet=TRUE)$property$value.1
          phot.tmp.weights<-generate.kohgroup.property(som=phot.tmp.som,
                                                       data=phot,
                                                       n.cluster.bins=tomo.bin[tomo],n.cores=64,
                                                       expression="sum(data[[blind.count.variable]])",quiet=TRUE)$property$value.1
          tmp.weights<-phot.tmp.weights/spec.tmp.weights
          tmp.weights<-list(refr.weight=1/tmp.weights[phot.tmp.som$clust.classif],train.weight=tmp.weights[spec.tmp.som$clust.classif])
          tmp.weights$refr.weight[which(!is.finite(tmp.weights$refr.weight))]<-0
          tmp.weights$train.weight[which(!is.finite(tmp.weights$train.weight))]<-0
          #/*fend*/}}}
          #Delta n_eff /*fold*/{{{
          cat(" $ ")
          cat(tround(digits=1,100*
              (sum(phot[[blind.count.variable]][which(tmp.weights$refr.weight!=0)])^2/
               sum(phot[[blind.count.variable]][which(tmp.weights$refr.weight!=0)]^2))/
              (sum(phot[[blind.count.variable]][phot.index])^2/
               sum(phot[[blind.count.variable]][phot.index]^2))
              ))
          cat(" $ & ")
          #/*fend*/}}}
          #Represented by nQ={3,4} == {1,2} /*fold*/{{{
          phot[[paste0("RedshiftGoldClass_",blind)]][phot.index]<-phot[[paste0("RedshiftGoldClass_",blind)]][phot.index]+ifelse(tmp.weights$refr.weight[phot.index]!=0,1,0)
          spec[[paste0("RedshiftGoldClass_",blind)]][spec.index]<-spec[[paste0("RedshiftGoldClass_",blind)]][spec.index]+ifelse(tmp.weights$train.weight[spec.index]!=0,1,0)
          #/*fend*/}}}
        }
        cat("\n")
        #/*fend*/}}}
      }
      #/*fend*/}}}
      #/*fend*/}}}
    }
  }
  #/*fend*/}}}
  #Output the catalogues /*fold*/ {{{
  outname<-outputpath("_allgoldclass",phot.catalogue,output.folder)
  write.file(file=outname,phot)
  outname<-outputpath("_allgoldclass",spec.catalogue,output.folder)
  write.file(file=outname,spec)
  #/*fend*/}}}
} else { 
  #Load the GOLD Flag catalogue /*fold*/{{{
  cat("Output File does exist!\n")
  outname<-outputpath("_allgoldclass",phot.catalogue,output.folder)
  cat("Reading goldclasses photometric sample from:\n")
  cat(paste(outname,"\n"))
  phot.out<-read.file(file=outname,extname='OBJECTS')

  outname<-outputpath("_allgoldclass",spec.catalogue,output.folder)
  cat("Reading goldclasses spectroscopic sample from:\n")
  cat(paste(outname,"\n"))
  spec.out<-read.file(file=outname,extname='OBJECTS')
  if (nrow(phot)!=nrow(phot.out)) { stop("The input and 'output' phot catalogues have different length!") }
  if (nrow(spec)!=nrow(spec.out)) { stop("The input and 'output' phot catalogues have different length!") }
  phot<-phot.out
  spec<-spec.out
  #Flag any bad GAAP magnitudes /*fold*/{{{
  phot.good.rows<-matrixStats::rowAlls(phot[,paste0('MAG_GAAP_',c('u','g','r','i','Z','Y','J','H','Ks'))]>0)
  #/*fend*/}}}
  rm('phot.out')
  rm('spec.out')
  #/*fend*/}}}
}

#Construct the Nz for each tomographic bin and GoldClass /*fold*/{{{
cat("\n\n GoldClass & Blind & Tomographic Bin & Delta n_{\\rm eff} \\\n")
for (blind in blind.list) { 
  for (set in 1:nrow(class.sets)) { 
    #Initialise the weights and flags /*fold*/ {{{
    spec[[paste0("Weight_SOM_",class.sets[set,1],"_",blind)]]<-0
    phot[[paste0("Flag_SOM_",class.sets[set,1],"_",blind)]]<-0
    spec[[paste0("GroupFactor_SOM_",class.sets[set,1],"_",blind)]]<-0
    phot[[paste0("GroupFactor_SOM_",class.sets[set,1],"_",blind)]]<-0
    #/*fend*/}}}
    #Initialise the expressions /*fold*/ {{{ 
    #Count Variable {{{
    if (blind != "NONE") { 
      blind.count.variable<-paste0(count.variable,'_',blind)
    } else { 
      blind.count.variable<-count.variable
    }
    #}}}
    spec.redshift<-class.sets[set,4]
    #QC Expression {{{
    goldset.QC.expr<-gsub(spec.variable,spec.redshift,QC.expr)
    blind.QC.expr<-gsub(count.variable,blind.count.variable,goldset.QC.expr)
    write(paste0(class.sets[set,1],": ",blind.QC.expr),stderr())
    split.expr<-split.expr(blind.QC.expr,ignore='abs')
    split.names<-names(split.expr$components)
    keep<-rep(TRUE,length(split.names))
    for (ind in 1:length(split.names)) {
      keep[ind]<-grepl(split.names[ind],split.expr$replace.expr)
    }
    split.expr$components<-split.expr$components[keep]

    if (any(grepl("spec",split.expr$components)&grepl("phot",split.expr$components))) { 
      stop("Spec and Phot expressions are mixed together!")
    }
    spec.expressions<-split.expr$components[which(grepl("spec",split.expr$components)&!grepl("full.spec",split.expr$components))]
    phot.expressions<-split.expr$components[which(grepl("phot",split.expr$components)&!grepl("full.phot",split.expr$components))]
    full.spec.expressions<-split.expr$components[which(grepl("full.spec",split.expr$components))]
    full.phot.expressions<-split.expr$components[which(grepl("full.phot",split.expr$components))]
    #}}}
    #Update the expressions for kohgroup.property {{{
    spec.expressions<-gsub("spec$","data$",spec.expressions,fixed=T)
    phot.expressions<-gsub("phot$","data$",phot.expressions,fixed=T)
    spec.expressions<-gsub("spec[[","data[[",spec.expressions,fixed=T)
    phot.expressions<-gsub("phot[[","data[[",phot.expressions,fixed=T)
    full.spec.expressions<-gsub("full.spec$","spec$",full.spec.expressions,fixed=T)
    full.phot.expressions<-gsub("full.phot$","phot$",full.phot.expressions,fixed=T)
    full.spec.expressions<-gsub("full.spec[[","spec[[",full.spec.expressions,fixed=T)
    full.phot.expressions<-gsub("full.phot[[","phot[[",full.phot.expressions,fixed=T)
    #}}}
    #Spectroscopic QC Expression {{{
    #blind.spec.expression<-gsub(count.variable,blind.count.variable,spec.expressions)
    blind.spec.expression<-spec.expressions
    blind.spec.expression<-gsub("SurveyFlag",paste0("SurveyFlag_",blind),blind.spec.expression)
    blind.spec.expression<-gsub("Class",paste0("Class_",blind),blind.spec.expression)
    #}}}
    #Photometric QC Expression {{{
    #blind.phot.expression<-gsub(count.variable,blind.count.variable,phot.expressions)
    blind.phot.expression<-phot.expressions
    blind.phot.expression<-gsub("Flag",paste0("Flag_",blind),blind.phot.expression)
    blind.phot.expression<-gsub("Class",paste0("Class_",blind),blind.phot.expression)
    #}}}
    #Spectroscopic Selection Expression {{{
    blind.spec.selection<-gsub(count.variable,blind.count.variable,class.sets[set,3])
    blind.spec.selection<-gsub("SurveyFlag",paste0("SurveyFlag_",blind),blind.spec.selection)
    blind.spec.selection<-gsub("Class",paste0("Class_",blind),blind.spec.selection)
    #}}}
    #Photometric Selection Expression {{{
    blind.phot.selection<-gsub(count.variable,blind.count.variable,class.sets[set,2])
    blind.phot.selection<-gsub("Flag",paste0("Flag_",blind),blind.phot.selection)
    blind.phot.selection<-gsub("Class",paste0("Class_",blind),blind.phot.selection)
    #}}}
    #Check the spec and phot selection expressions {{{
    if (length(eval(parse(text=blind.spec.selection)))==0) { 
      stop(paste("\nExpression: ",blind.spec.selection, "results in a length 0 vector!\n"))
    }
    if (length(eval(parse(text=blind.phot.selection)))==0) { 
      stop(paste("\nExpression: ",blind.phot.selection, "results in a length 0 vector!\n"))
    }
    #}}}}
    #/*fend*/}}}
    cat(paste(class.sets[set,1],'&',blind,'&'))
    for (tomo in 1:maxBin) { 
      #Define the output Nz name /*fold*/{{{
      outname.nz<-outputpath(paste0('_',class.sets[set,1],"_blind",blind,"_TOMO",tomo,"_Nz"),
                             phot.catalogue,output.folder=output.folder,format=nz.format)
      #/*fend*/}}}
      #Select the relevant photometric and spectroscopic data /*fold*/{{{
      tomo.index<-which(phot$Z_B>tomo.lim[tomo] & phot$Z_B<=tomo.lim[tomo+1] & 
                        phot[[blind.count.variable]]>0 & phot.good.rows)
      phot.index<-which(phot$Z_B>tomo.lim[tomo] & phot$Z_B<=tomo.lim[tomo+1] & 
                        phot[[blind.count.variable]]>0 & eval(parse(text=blind.phot.selection)) & phot.good.rows)
      spec.index<-which(spec$Z_B>tomo.lim[tomo] & spec$Z_B<=tomo.lim[tomo+1] & eval(parse(text=blind.spec.selection)))
      #/*fend*/}}}
      #Group the SOM pixels and calculate the properties /*fold*/{{{
      #Run the spec expressions {{{
      spec.tmp.vals<-generate.kohgroup.property(som=spec.fid.som,data=spec,n.cluster.bins=tomo.bin[tomo],n.cores=64,subset=spec.index,
                                                  expression=c(blind.spec.expression,"nrow(data)"),
                                                  expr.label=c(names(blind.spec.expression),"Nspec"),quiet=TRUE)
      spec.tmp.groups<-spec.tmp.vals$som$clust.classif
      qc.frame<-as.data.table(spec.tmp.vals$property)
      #}}}
      #Run the phot expressions {{{
      phot.tmp.vals<-generate.kohgroup.property(phot.fid.som,data=phot,n.cluster.bins=tomo.bin[tomo],n.cores=64,subset=phot.index,
                                                  expression=c(blind.phot.expression,'sum(data[[blind.count.variable]])'),
                                                  expr.label=c(names(blind.phot.expression),"sumLFweight"),quiet=TRUE)
      phot.tmp.groups<-phot.tmp.vals$som$clust.classif
      qc.frame<-as.data.table(cbind(qc.frame,phot.tmp.vals$property[,-which(colnames(phot.tmp.vals$property)=='group.id')]))
      #}}}
      if (any(grepl("full.spec",split.expr$components))) {
        #Run the training cat QC components {{{
        full.spec.tmp.vals<-NULL
        for (expr in full.spec.expressions) {
          full.spec.tmp.vals<-cbind(full.spec.tmp.vals,eval(parse(text=expr)))
        }
        colnames(full.spec.tmp.vals)<-names(split.expr$components)[which(grepl("full.spec",split.expr$components))]
        #}}}
        #Add to the QC frame {{{
        qc.frame<-as.data.table(cbind(qc.frame,full.spec.tmp.vals))
        #}}}
      }
      if (any(grepl("full.phot",split.expr$components))) {
        #Run the reference cat QC components {{{
        full.phot.tmp.vals<-NULL
        for (expr in full.phot.expression) {
          full.phot.tmp.vals<-cbind(full.phot.tmp.vals,eval(parse(text=expr)))
        }
        colnames(full.phot.tmp.vals)<-names(split.expr$components)[which(grepl("full.phot",split.expr$components))]
        #}}}
        #Add to the QC frame {{{
        qc.frame<-as.data.table(cbind(qc.frame,full.phot.vals))
        #}}}
      }
      #/*fend*/}}}
      #QC the SOM groups /*fold*/{{{
      qc.res.expr<-paste0("data.table(group.id=group.id,QCeval=",split.expr$replace.expr,")")
      qc.result<-qc.frame[,eval(parse(text=qc.res.expr))]
      #}}}
      #Check that the QC returned something {{{
      if (length(qc.result$QCeval)==0) {
        stop("QC result is length 0! There are no good associations remaining!")
      }
      #}}}
      #Check that the QC value is valid {{{
      if (all(is.na(qc.result$QCeval))) {
        stop("All QC results are NA! There are no good associations remaining!")
      }
      #}}}
      #Compute the weights with QC /*fold*/{{{
      phot.tmp.vals$property$sumLFweight[which(!qc.result$QCeval)]<-NA
      spec.tmp.vals$property$Nspec[which(!qc.result$QCeval)]<-NA
      tmp.clust.weights<-phot.tmp.vals$property$sumLFweight/spec.tmp.vals$property$Nspec
      tmp.weights<-list(train.weight=tmp.clust.weights[spec.tmp.vals$som$clust.classif],
                        refr.weight=1/tmp.clust.weights[phot.tmp.vals$som$clust.classif])
      #Check for non finite weights /*fold*/{{{
      if (any(!is.finite(tmp.weights$refr.weight))) { 
        warning("Setting non-finite refr.weight to 0!\n") 
        tmp.weights$refr.weight[which(!is.finite(tmp.weights$refr.weight))]<-0
      }
      if (any(!is.finite(tmp.weights$train.weight))) { 
        warning("Setting non-finite train.weight to 0!\n") 
        tmp.weights$train.weight[which(!is.finite(tmp.weights$train.weight))]<-0
      }
      #/*fend*/}}}
      #Print some dimension statistics {{{
      #print(str(tmp.weights))
      #print(length(spec.index))
      #print(length(phot.index))
      #print(length(phot[[paste0("Flag_SOM_",class.sets[set,1],"_",blind)]][phot.index]))
      #print(length(ifelse(tmp.weights$refr.weight>0,1,0)))
      #}}}
      #/*fend*/}}}
      #Print the neffective /*fold*/{{{
      cat(" $ ")
      cat(tround(digits=1,100*
          (sum(phot[[blind.count.variable]][which(tmp.weights$refr.weight!=0)])^2/
           sum(phot[[blind.count.variable]][which(tmp.weights$refr.weight!=0)]^2))/
          (sum(phot[[blind.count.variable]][tomo.index])^2/
           sum(phot[[blind.count.variable]][tomo.index]^2))
          ))
      cat(" $ & ")
      #/*fend*/}}}
      #Construct the Nz /*fold*/{{{
      nzdist<-weighted.hist(spec[[spec.redshift]],w=tmp.weights$train.weight/sum(tmp.weights$train.weight),
                    breaks=seq(0,6.001,by=0.05),plot=F)
      #/*fend*/}}}
      #Output the Nz #/*fold*/{{{
      write.file(file=outname.nz,data.frame(binstart=nzdist$breaks[1:length(nzdist$mids)],
                                         density=nzdist$density),quote=F,row.names=F,col.names=F)
      #/*fend*/}}}
      #Assign the spectroscopic weights /*fold*/{{{
      spec[[paste0("Weight_SOM_",class.sets[set,1],"_",blind)]][spec.index]<-tmp.weights$train.weight[spec.index]
      spec[[paste0("GroupFactor_SOM_",class.sets[set,1],"_",blind)]][spec.index]<-spec.tmp.groups[spec.index]
      #/*fend*/}}}
      #Assign the photometric Flags /*fold*/{{{
      phot[[paste0("Flag_SOM_",class.sets[set,1],"_",blind)]][phot.index]<-ifelse(tmp.weights$refr.weight[phot.index]>0,1,0)
      phot[[paste0("GroupFactor_SOM_",class.sets[set,1],"_",blind)]][phot.index]<-phot.tmp.groups[phot.index]
      #/*fend*/}}}
    }
    cat("\n")
  }
}
#/*fend*/}}}

#Output the catalogues /*fold*/ {{{
outname<-outputpath("_allgoldclass",phot.catalogue,output.folder)
write.file(file=outname,phot)
outname<-outputpath("_allgoldclass",spec.catalogue,output.folder)
write.file(file=outname,spec)
#/*fend*/}}}

#Finish 
