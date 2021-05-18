########################################安装包##################################

install.packages("usethis")
install.packages("devtools")
library(usethis)
library(devtools)
devtools::install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)
library("meta")
library(MRPRESSO)
##################################数据提取和处理################################

#diabetes     Nonalcoholic fatty liver disease 
ao <- available_outcomes()
a_exp_dat <- extract_instruments(outcomes='finn-a-NAFLD',clump = F,r2=0.001,kb=10000)
b_out_dat <- extract_outcome_data(snps = a_exp_dat$SNP,
                                      outcomes ='ebi-a-GCST006867' )
a_exp_dat <- extract_instruments(outcomes='bbj-a-101',clump = F,r2=0.01,kb=5000)
b_out_dat <- extract_outcome_data(snps = a_exp_dat$SNP,
                                  outcomes ='ebi-a-GCST006867')



a_exp_dat <- extract_instruments(outcomes='finn-a-H7_GLAUCOMA')
b_out_dat <- extract_outcome_data(snps = a_exp_dat$SNP,
                                  outcomes ='finn-a-H7_MYOPIA' )


a_exp_dat <- extract_instruments(outcomes='finn-a-H7_MYOPIA',clump = F,r2=0.01,kb=5000)
b_out_dat <- extract_outcome_data(snps = a_exp_dat$SNP,
                                  outcomes ='finn-a-H7_GLAUCOMA' )



dat<- harmonise_data(
  exposure_dat = a_exp_dat, 
  outcome_dat = b_out_dat
)
mr(
  dat,parameters = default_parameters(),method_list = subset(mr_method_list(), 
                                                             use_by_default)$obj
)
res<- mr(dat)
write.table(dat,file="F:/dat.csv",sep = ",",row.names=FALSE) 




read_exposure_data(
  filename,
  clump = FALSE,
  sep = " ",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  units_col = "units",
  ncase_col = "ncase",
  ncontrol_col = "ncontrol",
  samplesize_col = "samplesize",
  gene_col = "gene",
  id_col = "id",
  min_pval = 1e-200,
  log_pval = FALSE
)

#OR值
generate_odds_ratios(res)
res_single <- mr_singlesnp(dat)

a=generate_odds_ratios(res_single)
lnor<-log(a[,"or"])

lnuci<-log(a[,"or_uci95"])

lnlci<- log(a[,"or_lci95"])

selnor<- (lnuci-lnlci)/(2*1.96)

metal <- metagen(lnor,selnor, studlab=res_single[,"SNP"], sm="OR")


setwd("F:/project/");

adata=read.table("a.txt")


##################################孟德尔随机化分析##############################

setwd("F:/project")
dat1<-read.csv('a.csv')
mr(
  dat,parameters = default_parameters(),method_list = subset(mr_method_list(), 
  use_by_default)$obj
)
res<- mr(dat)
mrinput<-mr_input(
  bx = dat$beta.exposure,
  bxse = dat$se.exposure,
  by = dat$beta.outcome,
  byse = dat$se.outcome,
  
)

mr_allmethods(mrinput, method = "all")
####################################敏感性分析##################################

#OR值
generate_odds_ratios(res)
#异质性检验
mr_heterogeneity(dat)
#水平多效性
mr_pleiotropy_test(dat)
#敏感性分析

##################################数据可视化####################################

#散点图
p1 <- mr_scatter_plot(res, dat)
p1[[1]]

#森林图
b=match(unique(dat$SNP), dat$SNP)
c=dat[b,]
res_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(res_single)
forest(res_single)
p2[[1]]

#留一法

res_loo <- mr_leaveoneout(c)
p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]]

#倒漏斗图

res_single <- mr_singlesnp(dat)
p4<- mr_funnel_plot(res_single)
p4[[1]]

#
devtools::install_github("rondolab/MR-PRESSO",force = TRUE)
mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
          SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE,
          DISTORTIONtest = TRUE, data = dat, NbDistribution = 1000,
          SignifThreshold = 0.05)
