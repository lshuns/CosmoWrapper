#
# Script to install the required R packages for CosmoWrapper
#

#Check for an R installation 
P_R=`which R`
if [ "${P_R}" == ""] 
then 
  echo ERROR: There is no R executable in the PATH
  exit 1 
fi 

#Install the required packages 
R --no-restore --no-save <<EOF
qrequire<-function(...) suppressWarnings(suppressPackageStartupMessages(require(...)))
require.and.load<-function(name,githubrep,force=FALSE) { 
  if (!qrequire(name,character.only=TRUE) || force) { 
    if (!missing(githubrep)) { 
      devtools::install_github(paste(githubrep,name,sep='/'),upgrade='always')
    } else { 
      install.packages(name,repos='https://cloud.r-project.org/')
    }
    if (grepl('/',name) & !missing(githubrep)) { 
      name<-rev(strsplit(name,'/')[[1]])[1]
    }
    if (!qrequire(name,character.only=TRUE)) { 
      stop(paste("Failed to install package",name))
    }
  }
}
require.and.load('devtools') 
require.and.load('RColorBrewer') 
require.and.load('matrixStats')
require.and.load('data.table') 
require.and.load('helpRfuncs','AngusWright') 
require.and.load('FITSio') 
require.and.load('foreach')
require.and.load('doParallel')
require.and.load('itertools')
require.and.load('astro')
require.and.load('kohonen/kohonen','AngusWright') 

EOF

#Finish 

