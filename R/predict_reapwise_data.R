



library(lfe)

load("M:/NewRRED_M_drive/CEB/AndrewCD/OMB_variability/cropdat_reapwise_for_REAP_NNET.Rda")

#data_clean
cropdat_w$total_precip <- rowSums(cropdat_w[,grep('prc', colnames(cropdat_w))])

#temp_data
w <- cropdat_w









