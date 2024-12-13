###### THIS SCRIPT INVOLVES A TEST RUN OF SKAT-O
#### TEST WITH THE PLINK OUTPUT FROM annotation_test.sh
### WITH HELP FROM https://cran.r-project.org/web/packages/SKAT/SKAT.pdf & https://groups.google.com/g/skat_slee/c/5nGua8r76QQ

library(SKAT)
library(data.table)

##### READ IN PLINK FILES
#out_test.bim, .bed, fmt_out_test.fam
plink_path="" #insert your plink_path here
bim_file<-file.path(paste0(plink_path,"out_test.bim"))
bed_file<-file.path(paste0(plink_path,"out_test.bed"))
fam_file<-file.path(paste0(plink_path,"pheno_test.fam"))#version with updated phenotypes and sex
set_id<-file.path(paste0(plink_path,"rsid_list_set.txt"))#setID file from annotation.sh
ssd_path<-file.path(paste0(plink_path,"test.ssd"))
file_info_path<-file.path(paste0(plink_path,"file_info.txt"))
cov_path<- paste0(plink_path,"test_pca.eigenvec")

###################### STEP 1: Generate_SSD_SetID to create SNP sets
####### WORKS WELL
Generate_SSD_SetID(bed_file, bim_file,fam_file,set_id,
  paste0(plink_path,"test.ssd"), #SSD FILE
  paste0(plink_path,"file_info.txt")) #SSD INFO FILE

##################### STEP 2: Get PLINK info IN FORMAT FOR THE SKAT TEST
# ssd file from previous step
ssd_file<- Open_SSD(ssd_path, file_info_path)
#plink fam file
fam_read<-Read_Plink_FAM(fam_file, Is.binary=TRUE)
phen_dt<-fam_read$Phenotype
# plink covariate file
cov_input<-Read_Plink_FAM_Cov(fam_file, cov_path, Is.binary=TRUE, flag1=0, cov_header=FALSE) #recode phenotype

###################### STEP 3: RUN SKAT ON NULL MODEL
y1<-as.factor(cov_input$Phenotype)
x1<-cov_input$COV1
x2<-cov_input$COV2
x3<-cov_input$COV3
 
#case 2 control 1
obj<-SKAT_Null_Model(cov_input[,6] ~ cov_input[,1]+cov_input[,2]+cov_input[,3],  #trait ~ covariates matrix #couuld try with cov1:3 too
                 data=cov_input,
                 out_type="D", #dichotomous file
                Adjustment=TRUE) #covariates
#cov_input[,6] ~ cov_input[1]
SKAT.SSD.All(ssd_file, obj, kernel = "linear.weighted") #the number of markers to test for an association after excluding non-polymorphic or high missing rates markers 

###################### STEP 4: GET MAP
skatbin_res<-SKATBinary.SSD.All(ssd_file, obj)

Get_EffectiveNumberTest(skatbin_res$results[,8], alpha=0.05, Is.MidP=TRUE)

###################### STEP 5: MAP QQPLOT
QQPlot_Adj(skatbin_res$results[,2],skatbin_res$results[,8],main="QQ plot")




