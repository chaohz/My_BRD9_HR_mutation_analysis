 /opt/R/3.5.1/bin/R

setwd("/home/chao/A3B/MAF/data/OV/mutation")
#===================================================================
#===================================================================
mutData1 = read.csv(dir()[5], header=F, sep="\t", stringsAsFactors=F)
mutData2=mutData1[-(1:5),1]
rm(mutData1)
mutData3=array(mutData2,dim=c(120,length(mutData2)/120))
rm(mutData2)
mutData=as.data.frame(t(mutData3),stringsAsFactors=F)
rm(mutData3)
mutM = mutData[-1,]
names(mutM) = as.vector(as.matrix(mutData[1,]))
rm(mutData)

MCWmotif = mutM

patients_with_BRD9_M_OV = substr(mutM[which(mutM[ ,"SYMBOL"] == "BRD9"), "Tumor_Sample_Barcode"] ,1 ,15)
patients_with_BRD9_M_OV_index = patients_with_BRD9_M_OV[!duplicated(patients_with_BRD9_M_OV)]
patients_with_BRCA_M_OV = substr(mutM[which((mutM[ ,"SYMBOL"] == "BRCA1") | (mutM[ ,"SYMBOL"] == "BRCA2")), "Tumor_Sample_Barcode"] ,1 ,15)
patients_with_BRCA_M_OV_index = patients_with_BRCA_M_OV[!duplicated(patients_with_BRCA_M_OV)]

patients_with_M_OV = substr(mutM[, "Tumor_Sample_Barcode"] ,1 ,15)
patients_with_M_OV_index = patients_with_M_OV[!duplicated(patients_with_M_OV)]


setwd("/data/BRCA/BRCA_SNV_MUTECT2")
#===================================================================
#===================================================================
mutData1 = read.csv(dir()[5], header=F, sep="\t", stringsAsFactors=F)
mutData2=mutData1[-(1:5),1]
rm(mutData1)
mutData3=array(mutData2,dim=c(120,length(mutData2)/120))
rm(mutData2)
mutData=as.data.frame(t(mutData3),stringsAsFactors=F)
rm(mutData3)
mutM = mutData[-1,]
names(mutM) = as.vector(as.matrix(mutData[1,]))
rm(mutData)

MCWmotif = mutM

patients_with_BRD9_M_BR = substr(mutM[which(mutM[ ,"SYMBOL"] == "BRD9"), "Tumor_Sample_Barcode"] ,1 ,15)
patients_with_BRD9_M_BR_index = patients_with_BRD9_M_BR[!duplicated(patients_with_BRD9_M_BR)]
patients_with_BRCA_M_BR = substr(mutM[which((mutM[ ,"SYMBOL"] == "BRCA1") | (mutM[ ,"SYMBOL"] == "BRCA2")), "Tumor_Sample_Barcode"] ,1 ,15)
patients_with_BRCA_M_BR_index = patients_with_BRCA_M_BR[!duplicated(patients_with_BRCA_M_BR)]

patients_with_M_BR = substr(mutM[, "Tumor_Sample_Barcode"] ,1 ,15)
patients_with_M_BR_index = patients_with_M_BR[!duplicated(patients_with_M_BR)]


#setwd("E:/TCGA/data - 副本/OV/mutation")
#setwd("/home/chao/A3B/MAF/data/BRCA/BRCA_SNV_MUTECT2")
setwd("/home/chao/A3B/MAF/data/PAAD/mutation")
#===================================================================
#===================================================================
mutData1 = read.csv(dir()[5], header=F, sep="\t", stringsAsFactors=F)
mutData2=mutData1[-(1:5),1]
rm(mutData1)
mutData3=array(mutData2,dim=c(120,length(mutData2)/120))
rm(mutData2)
mutData=as.data.frame(t(mutData3),stringsAsFactors=F)
rm(mutData3)
mutM = mutData[-1,]
names(mutM) = as.vector(as.matrix(mutData[1,]))
rm(mutData)

MCWmotif = mutM

patients_with_BRD9_M_PA = substr(mutM[which(mutM[ ,"SYMBOL"] == "BRD9"), "Tumor_Sample_Barcode"] ,1 ,15)
patients_with_BRD9_M_PA_index = patients_with_BRD9_M_PA[!duplicated(patients_with_BRD9_M_PA)]
patients_with_BRCA_M_PA = substr(mutM[which((mutM[ ,"SYMBOL"] == "BRCA1") | (mutM[ ,"SYMBOL"] == "BRCA2")), "Tumor_Sample_Barcode"] ,1 ,15)
patients_with_BRCA_M_PA_index = patients_with_BRCA_M_PA[!duplicated(patients_with_BRCA_M_PA)]

patients_with_M_PA = substr(mutM[, "Tumor_Sample_Barcode"] ,1 ,15)
patients_with_M_PA_index = patients_with_M_PA[!duplicated(patients_with_M_PA)]

#setwd("E:/TCGA/data - 副本/OV/mutation")
#setwd("/home/chao/A3B/MAF/data/BRCA/BRCA_SNV_MUTECT2")
setwd("/home/chao/A3B/MAF/data/PRAD/mutation")
#===================================================================
#===================================================================
mutData1 = read.csv(dir()[5], header=F, sep="\t", stringsAsFactors=F)
mutData2=mutData1[-(1:5),1]
rm(mutData1)
mutData3=array(mutData2,dim=c(120,length(mutData2)/120))
rm(mutData2)
mutData=as.data.frame(t(mutData3),stringsAsFactors=F)
rm(mutData3)
mutM = mutData[-1,]
names(mutM) = as.vector(as.matrix(mutData[1,]))
rm(mutData)

MCWmotif = mutM


patients_with_BRD9_M_PR = substr(mutM[which(mutM[ ,"SYMBOL"] == "BRD9"), "Tumor_Sample_Barcode"] ,1 ,15)
patients_with_BRD9_M_PR_index = patients_with_BRD9_M_PR[!duplicated(patients_with_BRD9_M_PR)]
patients_with_BRCA_M_PR = substr(mutM[which((mutM[ ,"SYMBOL"] == "BRCA1") | (mutM[ ,"SYMBOL"] == "BRCA2")), "Tumor_Sample_Barcode"] ,1 ,15)
patients_with_BRCA_M_PR_index = patients_with_BRCA_M_PR[!duplicated(patients_with_BRCA_M_PR)]

patients_with_M_PR = substr(mutM[, "Tumor_Sample_Barcode"] ,1 ,15)
patients_with_M_PR_index = patients_with_M_PR[!duplicated(patients_with_M_PR)]

#===================================================================
#combine four cancer types
#===================================================================
#patients_with_BRD9_M_index = union(union(patients_with_BRD9_M_OV_index, patients_with_BRD9_M_BR_index), union(patients_with_BRD9_M_PA_index, patients_with_BRD9_M_PR_index))
#patients_with_BRCA_M_index = union(union(patients_with_BRCA_M_OV_index, patients_with_BRCA_M_BR_index), union(patients_with_BRCA_M_PA_index, patients_with_BRCA_M_PR_index))
#patients_with_M_index = union(union(patients_with_M_OV_index, patients_with_M_BR_index), union(patients_with_M_PA_index, patients_with_M_PR_index))

patients_with_BRD9_M_index = patients_with_BRD9_M_PR_index
patients_with_BRCA_M_index = patients_with_BRCA_M_PR_index
patients_with_M_index = patients_with_M_PR_index

#HR_sig_num = fit_res$contribution["Signature.3", ]
#A3B_sig_num = fit_res$contribution["Signature.13", ]

#patients_with_BRD9_M_HRsig = intersect(patients_with_BRD9_M_index, names(HR_sig_num))
#patients_with_BRCA_M_HRsig = intersect(patients_with_BRCA_M_index, names(HR_sig_num))
patients_with_BRD9_M_HRsig = patients_with_BRD9_M_index
patients_with_BRCA_M_HRsig = patients_with_BRCA_M_index

patients_with_M_HRsig = patients_with_M_index
#HR_sig_num_with_M_HRsig = HR_sig_num[patients_with_M_HRsig]

patients_with_M_HRsig_high = patients_with_M_index


patients_with_M_HRsig_high_BRD9_BRCA = intersect(patients_with_M_HRsig_high, intersect(patients_with_BRD9_M_HRsig, patients_with_BRCA_M_HRsig))
patients_with_M_HRsig_high_BRD9 = intersect(patients_with_M_HRsig_high, patients_with_BRD9_M_HRsig)[- which(intersect(patients_with_M_HRsig_high, patients_with_BRD9_M_HRsig) %in% patients_with_M_HRsig_high_BRD9_BRCA)]
patients_with_M_HRsig_high_BRCA = intersect(patients_with_M_HRsig_high, patients_with_BRCA_M_HRsig)[- which(intersect(patients_with_M_HRsig_high, patients_with_BRCA_M_HRsig) %in% patients_with_M_HRsig_high_BRD9_BRCA)]
patients_with_M_HRsig_high_null = patients_with_M_HRsig_high[- which(patients_with_M_HRsig_high %in% union(union(patients_with_M_HRsig_high_BRD9_BRCA, patients_with_M_HRsig_high_BRD9), patients_with_M_HRsig_high_BRCA))]



BRD9_mutation_array = data.frame(rep(0,length(patients_with_M_HRsig_high)))
rownames(BRD9_mutation_array ) = patients_with_M_HRsig_high
names(BRD9_mutation_array ) = "BRD9"
BRD9_mutation_array[union(patients_with_M_HRsig_high_BRD9, patients_with_M_HRsig_high_BRD9_BRCA),] = 1

BRCA_mutation_array = data.frame(rep(0,length(patients_with_M_HRsig_high)))
rownames(BRCA_mutation_array) = patients_with_M_HRsig_high
names(BRCA_mutation_array) = "BRCA"
BRCA_mutation_array[union(patients_with_M_HRsig_high_BRCA, patients_with_M_HRsig_high_BRD9_BRCA),] = 1

mydata = data.frame(BRD9_mutation_array, BRCA_mutation_array)
mydata = table(mydata)

#chisq.test(mydata)
fisher.test(mydata,alternative ="less")
