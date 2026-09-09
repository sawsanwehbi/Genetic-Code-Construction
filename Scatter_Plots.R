library(common)
library(ggplot2)
library(broom)

# Data from Wehbi et al. 2024 available in sawsanwehbi/Pfam-age-classification
Pfam_ConAAC <- read.csv('Pfam_data_ancestralAAC.csv', header = T)
Clan_ConAAfreq_Df <- read.csv('Clan_data_ancestralAAC.csv')
# Some Post-LUCA pfams reclassified as pre-LACA or pre-LBCA based on the presence of 
# 2 LBCA or LACA nodes branching near the root.
AncientPostLUCA <- read.csv('AncientPostLUCA.csv', header = T)


#### Pre-LACA and Pre-LBCA clan classification ######
postLUCAclans <- Clan_ConAAfreq_Df$Clans[which(Clan_ConAAfreq_Df$Clan_ancestor == 'postLUCA')]
LBCAclans <- vector()
LACAclans <- vector()
for ( i in 1:length(postLUCAclans)){
  if (any(Pfam_ConAAC$ancestor[which(Pfam_ConAAC$clans == postLUCAclans[i])] =='LBCA'))
  {LBCAclans[i] <- postLUCAclans[i]} 
  else if (any(Pfam_ConAAC$ancestor[which(Pfam_ConAAC$clans == postLUCAclans[i])] =='LACA'))
  {LACAclans[i] <- postLUCAclans[i]} else {next} }
LBCAclans <- LBCAclans[-which(is.na(LBCAclans))]
LACAclans <- LACAclans[-which(is.na(LACAclans))]
BroaddiversifiedLACAclans <- LACAclans[which(grepl('CL....', LACAclans))] #20
BroaddiversifiedLBCAclans <- LBCAclans[which(grepl('CL....', LBCAclans))] #156
#Broaddiversifiedclans <- c(BroaddiversifiedLBCAclans, BroaddiversifiedLACAclans) 
## broader classification of ancient postLUCA

## PostLUCA clans with >1 LACA/LBCA or 1 preLACA/preLBCA pfam; stricter classification
MultipleLACAclans <- sapply(1:length(BroaddiversifiedLACAclans), function(i){
  length(which(Pfam_ConAAC$ancestor[which(Pfam_ConAAC$clans %in% BroaddiversifiedLACAclans[i])] == 'LACA'))
})
MultipleLBCAclans <- sapply(1:length(BroaddiversifiedLBCAclans), function(i){
  length(which(Pfam_ConAAC$ancestor[which(Pfam_ConAAC$clans %in% BroaddiversifiedLBCAclans[i])] == 'LBCA'))
})
diversifiedLACAclans <- BroaddiversifiedLACAclans[which(MultipleLACAclans  > 1)] #4
diversifiedLACAclans <-  unique(c(diversifiedLACAclans ,unique(Pfam_ConAAC$clans[which(Pfam_ConAAC$pfamIDs %in% AncientPostLUCA$PFAM_IDs[which(
  AncientPostLUCA$classified_ancestor == 'preLACA')])]) )) #53
diversifiedLACAclans <- Clan_ConAAfreq_Df$Clans[which( Clan_ConAAfreq_Df$Clans %in% diversifiedLACAclans & 
                                                         Clan_ConAAfreq_Df$Clan_ancestor == 'postLUCA')] #31

diversifiedLBCAclans <- BroaddiversifiedLBCAclans[which(MultipleLBCAclans  > 1)] #77
diversifiedLBCAclans <-  unique(c(diversifiedLBCAclans , unique(Pfam_ConAAC$clans[which(Pfam_ConAAC$pfamIDs %in% AncientPostLUCA$PFAM_IDs[which(
  AncientPostLUCA$classified_ancestor == 'preLBCA')])])))#364
diversifiedLBCAclans <- Clan_ConAAfreq_Df$Clans[which( Clan_ConAAfreq_Df$Clans %in% diversifiedLBCAclans & 
                                                         Clan_ConAAfreq_Df$Clan_ancestor == 'postLUCA')] #251
diversifiedclans <- c(diversifiedLBCAclans, diversifiedLACAclans) 
# diversifiedLBCAclans and diversifiedLACAclans are preLBCA and preLACA clans

BroaddiversifiedLACAclans <- unique(c(BroaddiversifiedLACAclans, diversifiedLACAclans )) #39
BroaddiversifiedLBCAclans <- unique(c(BroaddiversifiedLBCAclans, diversifiedLBCAclans ))  #315
Broaddiversifiedclans <- c(BroaddiversifiedLBCAclans, BroaddiversifiedLACAclans) 
# Broaddiversifiedclans include diverged post-LUCA clans that do not exactly fit or pre-LACA/pre-LBCA criteria
# but we are unsure enough about their age that we can exclude them when using postLUCA AAC for usage ratios

## Revised post-LUCA AAC and SE
# In Wehbi et al. 2024 we had 1232 post-LUCA clans in the denominator
# After filtering 39 Diverged LACA and 315 diverged LBCA clans (which include pre-LACA & preLBCA),
# we now have 878 post-LUCA clans in denominator
NoDivpostLUCA_Clans_ConAAC <- colWeightedMeans(as.matrix(Clan_ConAAfreq_Df[which(Clan_ConAAfreq_Df$Clan_ancestor == 'postLUCA' &
                               !Clan_ConAAfreq_Df$Clans %in% Broaddiversifiedclans ),2:21]),
                               Clan_ConAAfreq_Df$Clan_Conlength[which(Clan_ConAAfreq_Df$Clan_ancestor == 'postLUCA' &
                                !Clan_ConAAfreq_Df$Clans %in% Broaddiversifiedclans)] )
LUCA_Clans_ConAAC <- colWeightedMeans(as.matrix(Clan_ConAAfreq_Df[which(Clan_ConAAfreq_Df$Clan_ancestor == 'LUCA'),2:21]),
                                      Clan_ConAAfreq_Df$Clan_Conlength[which(Clan_ConAAfreq_Df$Clan_ancestor == 'LUCA')] )
##SE
LUCAweightedse <- vector()
for (colnb in 2:21) {
  LUCAweightedse[colnb] <- weighted_se(Clan_ConAAfreq_Df[,colnb][which(Clan_ConAAfreq_Df$Clan_ancestor == 'LUCA')],
                                       Clan_ConAAfreq_Df$Clan_Conlength[which(Clan_ConAAfreq_Df$Clan_ancestor == 'LUCA')])}
NoDivpostLUCAweightedse <- vector()
for (colnb in 2:21){
  NoDivpostLUCAweightedse[colnb] <- weighted_se(Clan_ConAAfreq_Df[,colnb][which(Clan_ConAAfreq_Df$Clan_ancestor == 'postLUCA' &
                                    !Clan_ConAAfreq_Df$Clans %in% Broaddiversifiedclans )],
                                                Clan_ConAAfreq_Df$Clan_Conlength[which(Clan_ConAAfreq_Df$Clan_ancestor == 'postLUCA' &
                                                                                         !Clan_ConAAfreq_Df$Clans %in% Broaddiversifiedclans)])}

## Revised LUCA usage and SE
NoDiv_LUCAusage <- LUCA_Clans_ConAAC/NoDivpostLUCA_Clans_ConAAC
NodivLUCAclanratio_var <-  (LUCAweightedse[-1]^2)/((NoDivpostLUCA_Clans_ConAAC)^2) +
  (NoDivpostLUCAweightedse[-1]^2)*((LUCA_Clans_ConAAC)^2)/((NoDivpostLUCA_Clans_ConAAC)^4)
NodivLUCAclanratio_se <- sqrt(NodivLUCAclanratio_var)

NodivLUCAclanratio_se <- NodivLUCAclanratio_se[match(AA_properties$Letter,names(NodivLUCAclanratio_se))]
NoDiv_LUCAusage <- NoDiv_LUCAusage[match(  AA_properties$Letter, names(NoDiv_LUCAusage))]
# NoDiv_LUCAusage & NodivLUCAclanratio_se are renamed as LUCA_usage & LUCA_usageSE in AA_properties_Fig1_SupFig1.csv


#### Manuscript Figures ####
AA_properties_plots <- read.csv( '../Tryptophan paper/AA_properties_Fig1_SupFig1.csv')

## Figure 1A plot
pd <- position_dodge(width = 0.5)
NoDivNSA_model <- lm( AA_properties_plots$LUCA_usage ~ AA_properties_plots$Nb_nonH_sidechain,
                      weights = 1/(AA_properties_plots$LUCA_usageSE^2))
NoDivNSA_confinterval <- broom::augment(NoDivNSA_model, interval="confidence")
NoDivNSA_modelclass1 <- lm( AA_properties_plots$LUCA_usage[which(AA_properties_plots$aaRS_Class == 1 | AA_properties_plots$aaRS_Class == "1_2")] ~ 
                              AA_properties_plots$Nb_nonH_sidechain[which(AA_properties_plots$aaRS_Class == 1 | AA_properties_plots$aaRS_Class == "1_2")],
                            weights = 1/(AA_properties_plots$LUCA_usageSE[which(AA_properties_plots$aaRS_Class == 1 | AA_properties_plots$aaRS_Class == "1_2")]^2))
NoDivNSA_modelclass2 <- lm( AA_properties_plots$LUCA_usage[which(AA_properties_plots$aaRS_Class == 2 | AA_properties_plots$aaRS_Class == "1_2")] ~ 
                              AA_properties_plots$Nb_nonH_sidechain[which(AA_properties_plots$aaRS_Class == 2 | AA_properties_plots$aaRS_Class == "1_2")],
                            weights = 1/(AA_properties_plots$LUCA_usageSE[which(AA_properties_plots$aaRS_Class == 2 | AA_properties_plots$aaRS_Class == "1_2")]^2))
NoDivNSA_model_confintervalclass1 <- broom::augment(NoDivNSA_modelclass1, interval="confidence")
NoDivNSA_model_confintervalclass2 <- broom::augment(NoDivNSA_modelclass2, interval="confidence")
summary(NoDivNSA_model)

NoDivLUCA_NSA_plot <- ggplot(AA_properties_plots, aes(y = as.numeric(LUCA_usage), 
   x = as.numeric(Nb_nonH_sidechain), label = AA)) + 
  xlab('Number of non-H in side chain') + ylab('LUCA clan usage')  + 
  theme(legend.position="none") + scale_x_continuous(breaks = c(0,2,4,6,8,10)) + 
  geom_text(color='blue', size = 12,  position = pd) +
  #annotate(geom="text", y=0.83, x=1, label="Class II", color="#619CFF", size = 12) +
  #annotate(geom="text", y=0.78, x=1, label="Class I", color="#F8766D", size = 12) +
  #annotate(geom="text", y=1.1, x=180, label="Class I & II", color="#00BA38", size = 12) +
  geom_errorbar(aes(ymin=as.numeric(LUCA_usage)-LUCA_usageSE, ymax=as.numeric(LUCA_usage)+LUCA_usageSE), position = pd) +
  theme(axis.text=element_text(size=24),axis.title=element_text(size=32,face="bold"),
        legend.text=element_text(size=14), legend.title = element_text(size=16,face="bold")) + 
  geom_line(aes(y = predict(NoDivNSA_model)), linewidth = 1,color = 'black') +
  annotate(geom="text", x=4, y=0.7, label=paste0(paste0('Weighted R',supsc('2')), "= 0.55"), color= "black", size = 14) +
  annotate(geom="text", x=4, y=0.64, label="p = 0.0002", color= "black", size = 14) +
  #annotate(geom="text", x=70, y=1.2, label="b)", color="black", size = 16) +
  geom_ribbon(aes(ymin=NoDivNSA_confinterval$.lower, ymax=NoDivNSA_confinterval$.upper), colour=NA, alpha=0.3)


### Figure 1B plot 
BiosynReqmodified <-  AA_properties_plots[- which(AA_properties_plots$AA == 'N' | AA_properties_plots$AA == 'Q'),]
# LysRS-I is presumed to precede LysRS-II, so K is treated as Class I
BiosynReqmodified$aaRS_Class[9] <- 1

ModifiedNoDivNSA_model <- lm(BiosynReqmodified$LUCA_usage ~ BiosynReqmodified$Nb_nonH_sidechain,
                             weights = 1/(BiosynReqmodified$LUCA_usageSE^2))
ModifiedNoDivNSA_modelclass2 <- lm(BiosynReqmodified$LUCA_usage[which(BiosynReqmodified$aaRS_Class == 2)] ~ 
                                     BiosynReqmodified$Nb_nonH_sidechain[which(BiosynReqmodified$aaRS_Class == 2)],
                                   weights = 1/(BiosynReqmodified$LUCA_usageSE[which(BiosynReqmodified$aaRS_Class == 2  )]^2))
ModifiedNoDivNSA_modelclass1 <- lm(BiosynReqmodified$LUCA_usage[which(BiosynReqmodified$aaRS_Class == 1 )] ~ 
                                     BiosynReqmodified$Nb_nonH_sidechain[which(BiosynReqmodified$aaRS_Class == 1  )],
                                   weights = 1/(BiosynReqmodified$LUCA_usageSE[which(BiosynReqmodified$aaRS_Class == 1 )]^2))
ModifiedNoDivNSA_model_confintervalclass1 <- broom::augment(ModifiedNoDivNSA_modelclass1, interval="confidence")
ModifiedNoDivNSA_model_confintervalclass2 <- broom::augment(ModifiedNoDivNSA_modelclass2, interval="confidence")
summary(ModifiedNoDivNSA_modelclass1 )

ModifiedNoDivLUCA_NSA_plot <- ggplot(BiosynReqmodified, aes(y = as.numeric(LUCA_usage), 
  x = as.numeric(Nb_nonH_sidechain), label = AA)) + 
  xlab('Number of non-H in side chain') + ylab('LUCA clan usage')  + 
  theme(legend.position="none") + geom_text(aes(colour = factor(aaRS_Class)), size = 16, position = pd) + 
  labs(color='aaRS Class') + scale_x_continuous(breaks = c(0,2,4,6,8,10)) +
  annotate(geom="text", y=0.8, x=1, label="Class II", color="#00BFC4", size = 12) +
  annotate(geom="text", y=0.75, x=1, label="Class I", color="#F8766D", size = 12) +
  geom_errorbar(aes(ymin=as.numeric(LUCA_usage)-LUCA_usageSE, ymax=as.numeric(LUCA_usage)+LUCA_usageSE), position = pd) +
  theme(axis.text=element_text(size=24),axis.title=element_text(size=32,face="bold"),
        legend.text=element_text(size=14), legend.title = element_text(size=16,face="bold")) + 
  #geom_line(aes(y = predict(NoDivNSA_model)), linewidth = 1,color = 'black') +
  #annotate(geom="text", x=4, y=0.73, label=paste0(paste0('Weighted R',supsc('2')), "= 0.69"), color= "black", size = 14) +
  #annotate(geom="text", x=4, y=0.67, label="p = 5e-6", color= "black", size = 14) +
  #annotate(geom="text", x=70, y=1.2, label="b)", color="black", size = 16) +
  #geom_ribbon(aes(ymin=NoDivNSA_confinterval$.lower, ymax=NoDivNSA_confinterval$.upper), colour=NA, alpha=0.3)
  geom_line(data= BiosynReqmodified[which(BiosynReqmodified$aaRS_Class==1  ),] ,aes(y = predict(ModifiedNoDivNSA_modelclass1)), linewidth = 1,color = "#F8766D") +
  geom_line(data= BiosynReqmodified[which(BiosynReqmodified$aaRS_Class==2),], aes(y = predict(ModifiedNoDivNSA_modelclass2)), linewidth = 1,color = "#00BFC4") +
  annotate(geom="text", x=3, y=0.65, label=paste0(paste0('Weighted R',supsc('2')), "= 0.52"), color= "#00BFC4", size = 14) +
  annotate(geom="text", x=3, y=0.6, label="p = 0.04", color= "#00BFC4", size = 14) +
  annotate(geom="text", x=7, y=1.2, label=paste0(paste0('Weighted R',supsc('2')), "= 0.87"), color= "#F8766D", size = 14) +
  annotate(geom="text", x=7, y=1.15, label="p = 0.0002", color= "#F8766D", size = 14) +
  geom_ribbon(data=BiosynReqmodified[which(BiosynReqmodified$aaRS_Class==1  ),] , aes(ymin=ModifiedNoDivNSA_model_confintervalclass1$.lower, ymax=ModifiedNoDivNSA_model_confintervalclass1$.upper), colour="#F8766D", alpha=0.3) +
  geom_ribbon(data=BiosynReqmodified[which(BiosynReqmodified$aaRS_Class==2 ),] , aes(ymin=ModifiedNoDivNSA_model_confintervalclass2$.lower, ymax=ModifiedNoDivNSA_model_confintervalclass2$.upper), colour="#00BFC4", alpha=0.3)

## Supp Figure 1
NoDivConclan_wls_model <- lm(AA_properties_plots$Trifonov_order ~ AA_properties_plots$LUCA_usage)
NoDivConclan_wls_confinterval <- broom::augment(NoDivConclan_wls_model , interval="confidence")
summary(NoDivConclan_wls_model )
NoDivLUCAvsTrifonov_plot <- ggplot(AA_properties_plots, aes(y = as.numeric(Trifonov_order ), 
  x = as.numeric(LUCA_usage), label = AA)) + 
  ylab('Trifonov (2000) order') +
  xlab('LUCA clan usage') + 
  geom_text(color='blue', size = 12) +
  theme(axis.text=element_text(size=24),axis.title=element_text(size=32,face="bold"),
        legend.text=element_text(size=27), legend.title = element_text(size=16,face="bold")) + 
  scale_fill_brewer(palette="Dark2") +  theme(legend.position = 'none') +
  guides(color = guide_legend(override.aes = list(size = 16))) + 
  geom_line(aes(y = predict(NoDivConclan_wls_model)), linewidth = 1,color = 'black') +
  annotate(geom="text", y=26, x=0.9, label=paste0(paste0('R',supsc('2')), "= 0.27"), color="black", size = 14) +
  annotate(geom="text", y=23, x=0.9, label="p = 0.02", color="black", size = 14) +
  #annotate(geom="text", y=26, x=0.63, label="c)", color="black", size = 16) +
  geom_ribbon(aes(ymin=NoDivConclan_wls_confinterval$.lower, ymax=NoDivConclan_wls_confinterval$.upper), colour=NA, alpha=0.3)


### Figure 6
MinMeanDistances <- AA_properties_plots[-which( AA_properties_plots$AA == 'Q' | AA_properties_plots$AA == 'N'),]

cor.test(MinMeanDistances$mean_protozymedistance, MinMeanDistances$LUCA_usage, method='pearson')
cor.test(MinMeanDistances$min_protozymedistance, MinMeanDistances$LUCA_usage, method='pearson')
cor.test(MinMeanDistances$min_protozymedistance, MinMeanDistances$Recruitment_order, method='spearman')


min_model <- lm(MinMeanDistances$min_protozymedistance ~ MinMeanDistances$LUCA_usage)
mean_model <- lm(MinMeanDistances$mean_protozymedistance ~ MinMeanDistances$LUCA_usage)
# Pearson rho reported in figure
ggplot(MinMeanDistances, aes(x = -as.numeric(LUCA_usage),y = as.numeric(min_protozymedistance),
  label = AA)) +  geom_text( size = 10) +  ylab('Distance of aaRS from protozyme') +
  xlab('LUCA clan usage') +   xlim(-1.15, -0.6) + ylim(0,8) +
  scale_x_continuous(breaks = c(-1,-0.8,-0.6), labels = -c(-1.0,-0.8,-0.6)) +
  annotate(geom="text", x=-1, y=1, label=paste0("Mean"), color="purple", size = 11) +
  annotate(geom="text", x=-1, y=0.5, label=paste0("rho = 0.50"), color="purple", size = 11) +
  annotate(geom="text", x=-1, y=0, label="p = 0.03", color="purple", size = 11) + 
  annotate(geom="text", x=-0.8, y=1, label=paste0("Minimum"), color="#0071CE", size = 11) +
  annotate(geom="text", x=-0.8, y=0.5, label=paste0("rho = 0.38"), color="#0071CE", size = 11) +
  annotate(geom="text", x=-0.8, y=0, label="p = 0.1", color="#0071CE", size = 11) +
  geom_line(data= MinMeanDistances, aes(y = predict(min_model)), linewidth = 1,color = "#0071CE") +
  geom_line(data= MinMeanDistances, aes(y = predict(mean_model)), linewidth = 1,color = 'purple') +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=30,face="bold") ) 
# cosmetic edits to this plot where made in illustrator
