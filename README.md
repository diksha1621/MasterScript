# MasterScript



####caper PGLS model

#set libraries

library(phylosignal)

library(adephylo)

library(ape)

library(phylobase)

library(caper)



#read in file

CNV_df<- read.csv("NEXUS_for_phylosignal_small.csv")



#normallize genome size by the size of rDNA cassette 

cassette_size_in_bp<- CNV_df$rDNA_CN*9100 #assume 9.1kb for average single cassette size

genome_size_in_bp<- CNV_df$size*1000000 # Mbp convert to bp

genome_length_w_cassette<- genome_size_in_bp + cassette_size_in_bp #add the casette size to the genome size

df_with_rDNA_size1<-cbind(CNV_df, genome_size_in_bp)

df_with_rDNA_size<-cbind(df_with_rDNA_size1, genome_length_w_cassette)



#read in the tree

phy<- as.character(CNV_df[1,3])

phy<- read.tree(text = phy)



#check the data

length(phy$tip.label)



#read in the data

dat_WO_CN_length <- data.frame(taxa=df_with_rDNA_size$X, genome_size=df_with_rDNA_size$genome_size_in_bp, CNV=as.numeric(df_with_rDNA_size$rDNA_CN))

dat_W_CN_length <- data.frame(taxa=df_with_rDNA_size$X, genome_size=df_with_rDNA_size$genome_length_w_cassette, CNV=as.numeric(df_with_rDNA_size$rDNA_CN))



#match data to tips

cdat_WO_CN<- comparative.data(data=dat_WO_CN_length, phy=phy, names.col="taxa", vcv=TRUE)

cdat_W_CN<- comparative.data(data=dat_W_CN_length, phy=phy, names.col="taxa", vcv=TRUE)



#pgls

mod_WO <- pgls(formula = log(CNV) ~ log(genome_size), data = cdat_WO_CN)

summary(mod_WO)

mod_W <- pgls(formula= log(CNV) ~ log(genome_size), data = cdat_W_CN)

summary(mod_W)
