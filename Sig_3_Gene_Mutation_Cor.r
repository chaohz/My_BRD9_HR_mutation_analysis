 /opt/R/3.5.1/bin/R
load("/home/chao/A3B/MAF/data/UVM/mutation/96_trinucleodide_mutation.Rdata")

library(BSgenome)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)
library(MutationalPatterns)



sp_url <- paste("https://cancer.sanger.ac.uk/cancergenome/assets/",
 "signatures_probabilities.txt", sep = "")
 
cancer_signatures = read.table(sp_url, sep = "\t", header = TRUE)
new_order = match(row.names(mut_mat), cancer_signatures$Somatic.Mutation.Type)
cancer_signatures = cancer_signatures[as.vector(new_order),]
row.names(cancer_signatures) = cancer_signatures$Somatic.Mutation.Type
cancer_signatures = as.matrix(cancer_signatures[,4:33])

fit_res <- fit_to_signatures(mut_mat, cancer_signatures)
#fit_res$contribution
#fit_res$reconstructed


#================================================================================================
#================================================================================================
HR_sig_num = fit_res$contribution["Signature.3", ]

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

#patients_with_EGFR_M = substr(mutM[which((mutM[ ,"SYMBOL"] == "EGFR") & (mutM[ ,"Variant_Classification"] != "RNA") & (mutM[ ,"Variant_Classification"] != "Silent")), "Tumor_Sample_Barcode"] ,1 ,15)
#patients_with_EGFR_M = substr(mutM[which((mutM[ ,"SYMBOL"] == "EGFR") & (mutM[ ,"Variant_Classification"] != "Silent")), "Tumor_Sample_Barcode"] ,1 ,15)
patients_with_EGFR_M = substr(mutM[which(mutM[ ,"SYMBOL"] == "BRD9"), "Tumor_Sample_Barcode"] ,1 ,15)
patients_with_EGFR_M = patients_with_EGFR_M[!duplicated(patients_with_EGFR_M)]

patients_with_M = substr(mutM[, "Tumor_Sample_Barcode"] ,1 ,15)
patients_with_M = patients_with_M[!duplicated(patients_with_M)]

rm(mutM)

commonpatients = intersect(names(HR_sig_num), patients_with_M)
patients_EGFR_M = intersect(patients_with_EGFR_M, commonpatients)
HR_sig_num = HR_sig_num[commonpatients]

mutation_array = data.frame(rep(0,length(commonpatients)))
rownames(mutation_array) = commonpatients 
names(mutation_array) = "mutation"
mutation_array[patients_EGFR_M,] = 1

m_threshold = fivenum(HR_sig_num)[4]
#m_threshold = mean(HR_sig_num)
#m_threshold = sort(HR_sig_num)[round(length(HR_sig_num)*0.70)]

express_array = data.frame(rep(0,length(commonpatients)))
rownames(express_array) = commonpatients
names(express_array) = "expression"
express_array[which(HR_sig_num > m_threshold), ] = 1

#=====================================================================================
mydata = data.frame(mutation_array, express_array)
#mydata = table(mutation_array[,1], express_array[,1])
mydata = table(mydata)

#chisq.test(mydata)
fisher.test(mydata,alternative ="greater")
