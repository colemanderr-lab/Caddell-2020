#This is code used to make figures 4E and 4F
#author: "Elizabeth Purdom"

# Read in the readInMetabolite
metabolite<-read.csv("metabolite_peakheight_032420.csv",header=FALSE,skip=2,row.names=1)
metaboliteMeta<-read.csv("metabolite_peakheight_032420.csv",header=FALSE,nrow=2,row.names=1)
metaboliteMeta<-data.frame(t(metaboliteMeta))

metaboliteMeta$SampleType<-sapply(strsplit(as.character(metaboliteMeta$Sample),"-"),.subset2,1)
metaboliteMeta$SampleType<-factor(metaboliteMeta$SampleType,levels=c("Rhizo",  "Root",  "Soil"),labels=c("Rhizosphere",  "Root",  "Soil"))

metaboliteMeta$Timepoint<-sapply(strsplit(as.character(metaboliteMeta$Sample),"-"),.subset2,2)
metaboliteMeta$Timepoint<-factor(metaboliteMeta$Timepoint,levels=c("0hr","24hr"),labels=c("Week8","Post24hr"))
metaboliteMeta$Treatment<-sapply(strsplit(as.character(metaboliteMeta$Sample),"-"),.subset2,3)
metaboliteMeta$Treatment<-factor(metaboliteMeta$Treatment,levels=c("D","W"),labels=c("Drought","Watered"))
metaboliteMeta$Rep<-sapply(strsplit(as.character(metaboliteMeta$Sample),"-"),.subset2,4)
metaboliteMeta$SampleId<-with(metaboliteMeta,paste(SampleType,Timepoint, Treatment, Rep,sep="-"))

colnames(metabolite)<-metaboliteMeta$SampleId
rownames(metaboliteMeta)<-metaboliteMeta$SampleId
library(SummarizedExperiment)
metaboliteSE<-SummarizedExperiment(metabolite,colData=metaboliteMeta)

# > names(table(metaboliteMeta$Label))
#  [1] "Rhizo-0hr-D"  "Rhizo-0hr-W"  "Rhizo-24hr-D" "Rhizo-24hr-W" "Root-0hr-D"   "Root-0hr-W"
#  [7] "Root-24hr-D"  "Root-24hr-W"  "Soil-0hr-D"   "Soil-0hr-W"   "Soil-24hr-D"  "Soil-24hr-W"
# > names(table(microbe$Group))
#  [1] "TP24h_Drought_Rhizosphere" "TP24h_Drought_Root"        "TP24h_Drought_Soil"
#  [4] "TP24h_Watered_Rhizosphere" "TP24h_Watered_Root"        "TP24h_Watered_Soil"
#  [7] "TP8_Drought_Rhizosphere"   "TP8_Drought_Root"          "TP8_Drought_Soil"
# [10] "TP8_Watered_Rhizosphere"   "TP8_Watered_Root"          "TP8_Watered_Soil"


### Read in microbial
microbe<-read.csv("microbiome_ASV_abundance.csv")
microbe$Timepoint<-factor(microbe$Timepoint,levels=c(1,2),labels=c("Week8","Post24hr"))
microbe$SampleId<-with(microbe,paste(SampleType,Timepoint, Treatment, Rep,sep="-"))

allOTU<-as.character(sort(unique(microbe$OTU)))
OTUData<-do.call("cbind",by(microbe,as.factor(microbe$SampleId),function(x){
    x$Abundance[match(allOTU,as.character(x$OTU))]
},simplify=FALSE))
rownames(OTUData)<-allOTU
OTUMeta<-unique(microbe[,c("OTU","Domain","Phylum","Class","Order","Family","Genus","Species")])
OTUMeta<-OTUMeta[match(allOTU,OTUMeta$OTU),]
row.names(OTUMeta)<-OTUMeta[,1]
sampleMeta<-do.call("rbind",by(microbe,microbe$SampleId,function(x){
    unique(x[,c("Sample", "Date", "Variety", "Plot", "Treatment", "SampleType", "Rep", "Timepoint" , "Group","SampleTreatment","SampleId")])
},simplify=FALSE))

library(SummarizedExperiment)
OTUSE<-SummarizedExperiment(OTUData,colData=sampleMeta,rowData=OTUData)

#Combine them across phyllum etc. (could do this better perhaps with `phyloseq`)

#### Combine across phyllum etc:
combinePhylo<-function(columnName){
    out<-do.call("rbind",by(OTUData,OTUMeta[,columnName],colSums))
    missing<-grep("Unknown",rownames(out))
    if(length(missing)>0) return(out[-missing,])
    else return(out)
}
domainOTU<-combinePhylo("Domain")
phylumOTU<-combinePhylo("Phylum")
classOTU<-combinePhylo("Class")
orderOTU<-combinePhylo("Order")
familyOTU<-combinePhylo("Family")
genusOTU<-combinePhylo("Genus")
speciesOTU<-combinePhylo("Species")
# domain>kingdom>pylum>class>order>family>genus>species
library(MultiAssayExperiment)
allAssays<-MultiAssayExperiment(experiments=SimpleList(
    domain=domainOTU,
    phylum=phylumOTU,
    class=classOTU,
    order=orderOTU,
    family=familyOTU,
    genus=genusOTU,
    species=speciesOTU,
    microbe=OTUSE,
    metabolite=metaboliteSE),colData=sampleMeta)
experiments(allAssays)

#This is way slower than it should be!
takeDiff<-function(data,meta,type=c("Treatment","Timepoint")){
    type<-match.arg(type)
    
    #check all right names
    #Treatment  SampleType Rep Timepoint
    if(!all(c("SampleId","Treatment","SampleType","Rep","Timepoint") %in% colnames(meta))) stop("error missing names")
    if(!all(colnames(data) %in% meta[,"SampleId"])) stop("samples not in meta data")
    if(!all(meta[,"SampleId"] %in% colnames(data))) stop("samples not in data")
    
    # make in same order
    meta<-meta[match(colnames(data),meta[,"SampleId"]),]
    if(nrow(meta)!=ncol(data)) stop("error in aligning data")
    
    #Reorder, so can just subtract two matrices
    ord<-switch(type,"Treatment"=order(meta[,"SampleType"],meta[,"Timepoint"],meta[,"Rep"],meta[,"Treatment"]),
                "Timepoint"=order(meta[,"SampleType"],meta[,"Treatment"],meta[,"Rep"],meta[,"Timepoint"])
    )
    data<-data[,ord]
    meta<-meta[ord,]
    #group by everything BUT relevant variable
    fac<-switch(type,
                "Treatment"=paste(meta[,"SampleType"],meta[,"Timepoint"],meta[,"Rep"],sep="-"),
                "Timepoint"=paste(meta[,"SampleType"],meta[,"Treatment"],meta[,"Rep"],sep="-"))
    
    if(type=="Treatment"){
        if(!all(fac[meta[,"Treatment"]=="Drought"]==fac[meta[,"Treatment"]=="Watered"])) stop("error in subsetting data")
        out<-(data[,meta[,"Treatment"]=="Drought"]-data[,meta[,"Treatment"]=="Watered"])
        colnames(out)<-fac[meta[,"Treatment"]=="Drought"]
        
    }
    if(type=="Timepoint"){
        if(!all(fac[meta[,"Timepoint"]=="Post24hr"]==fac[meta[,"Timepoint"]=="Week8"])) stop("error in subsetting data")
        out<-(data[,meta[,"Timepoint"]=="Post24hr"]-data[,meta[,"Timepoint"]=="Week8"])
        colnames(out)<-fac[meta[,"Timepoint"]=="Post24hr"]
        
    }
    return(out)
}


## Combine Class and metabolite data
combinedClass<-rbind(log(classOTU+1),log(metabolite))
combinedClassDiffTreat<-takeDiff(combinedClass,sampleMeta,"Treatment")
combinedClassDiffTP<-takeDiff(combinedClass,sampleMeta,"Timepoint")


#Color palette for top 10 phyla
phyllumColor <- c(
  "Other"="#4d5b6a", # Other
  "TM7"="#599ec4" # TM7
  ,Nitrospirae="#a1c5c8" # Nitrospirae
  ,Verrucomicrobia="#7997a1" # Verrucomicrobia
  ,Firmicutes="#427e83" # Firmicutes
  ,Gemmatimonadetes="#84b59c" # Gemmatimonadetes
  ,Bacteroidetes="#919f70" # Bacteroidetes
  ,Chloroflexi="#686a47" # Chloroflexi
  ,Acidobacteria="#cb9e59" # Acidobacteria
  ,Proteobacteria="#a04344" # Proteobacteria
  ,Actinobacteria="#50414c") # Actinobacteria
phyllumFac<-OTUMeta$Phylum[match(row.names(classOTU),OTUMeta$Class)]
phyllumFac[!phyllumFac %in% names(phyllumColor)]<-"Other"
#Put NA for metabolite
phyllumFac<-c(phyllumFac,rep(NA,nrow(metabolite)))

#Colors for bacteria vs metabolite
bacMetab<-c(Bacteria = "#8d8089" ,Metabolites = "#e5e0db")
dataSource<-rep(c("Bacteria","Metabolites"),times=c(nrow(classOTU),nrow(metabolite)))


annRowClass<-data.frame(dataSource=as.factor(dataSource), Phyllum=phyllumFac)
rownames(annRowClass)<-rownames(combinedClass)
annColors<-list(dataSource=bacMetab, Phyllum=phyllumColor)


library(pheatmap)

# As far as heatmap dimensions, I believe the image for root expression you sent was 6x12 (uncropped). 
# I think we will use something closer to 8x12 for the figure 4. 
#' However, I am hoping to expand the width of 'of interest' and 'datasource'.
# The final relative dimensions of the plot will be closer to 
# 3 (tree):1 (of interest):1 (datasource) :5 (heatmap), as in the figure 4 outline. 
# If it is easier, I can just stretch the image to the desired dimensions after I get it from you. 

### All Samples
### TO DO: average of replicates
pdf("paperFigures/allTimespoints.pdf",width=8,height=12)
dataSet<-t(scale(t(combinedClass)))
bks<-clusterExperiment::setBreaks(dataSet,breaks=.999,makeSymmetric=TRUE)
fullTree<-pheatmap(dataSet,cluster_cols=TRUE,annotation_names_col=TRUE,main="",
                   scale="none",color=colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(50),
                   annotation_row=annRowClass,annotation_colors=annColors,annotation_legend=TRUE,
                   treeheight_row=120,breaks=bks,
                   fontsize_row=3, fontsize_col=6)
dev.off()

### Root only
pdf("paperFigures/RootOnly.pdf",width=8,height=12)
whSample<-"Root"
dataSet<-scale(t(combinedClass))
dataSet<-t(dataSet[grep(whSample,rownames(dataSet)),])
#bks<-clusterExperiment::setBreaks(dataSet,breaks=.999,makeSymmetric=TRUE)
rootTree<-pheatmap(dataSet,cluster_cols=TRUE,annotation_names_col=TRUE,main="",
                   scale="none",color=colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(50),
                   annotation_row=annRowClass,annotation_colors=annColors,annotation_legend=TRUE,
                   treeheight_row=120,breaks=bks,
                   fontsize_row=3, fontsize_col=6)
dev.off()

### Root Zoom
pdf("paperFigures/RootZoom.pdf",width=8,height=8)
rootClusterHeatmap<-cutree(rootTree$tree_row,k=65)
whActino<-rootClusterHeatmap[which(names(rootClusterHeatmap)=="Actinobacteria")]
names(rootClusterHeatmap[rootClusterHeatmap==whActino])
dataSet<-scale(t(combinedClass))
whSample<-"Root"
dataSet<-t(dataSet[grep(whSample,rownames(dataSet)),])
dataSet<-dataSet[rootClusterHeatmap==whActino,]
bks<-clusterExperiment::setBreaks(dataSet,breaks=.999,makeSymmetric=TRUE)
rootZoomTree<-pheatmap(dataSet,cluster_cols=rootTree$tree_col,annotation_names_col=TRUE,main="",
                   scale="none",color=colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(50),
                   treeheight_row=120,breaks=bks,treeheight_col=0,
                   fontsize_row=6, fontsize_col=6)
dev.off()
