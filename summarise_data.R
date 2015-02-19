
# The models available are: strict_constant.log, strict_bs.log, ucld_constant.log, ucld_bs.log

#STORE OUTPUT IN SQLITE 3. MIGHT NEED TO RUN MLE ESTIMATES LONGER...


# For each folder: read log files and check diagnostics, save for each model the mean rate parameter (clockrate or ucldMean), and MLE

source('~/Desktop/bact_rates_scripts/proc_data.R')
beast_path <- '~/Desktop/phylo_programs/beast181/bin/beast'


########
## PROCESS CHOLERA ANALYSES
#######

##############
#Set up files and output matrix
dir_runs <- 'bacteria_rates/cholera/runs/'
all_log_files <- grep('log', dir(dir_runs), value = T)
mle_files <- grep('mle', all_log_files)
prior_files <- grep('prior', all_log_files)

cholera_results <- matrix(NA, 4, 5)
rownames(cholera_results) <- c('strict_constant', 'strict_bs', 'ucld_constant', 'ucld_bs')
colnames(cholera_results) <- c('Rate', 'Lower95', 'Upper95', 'MLE_path', 'MLE_ss')
cholera_results <- as.data.frame(cholera_results)


#strict_constant
s_c_log <- grep('strict_constant', all_log_files[-c(mle_files, prior_files)], value = TRUE)
s_c_mle <- grep('strict_constant', all_log_files[mle_files], value = TRUE)

s_c_diagnose <- diagnose_mcmc(read.table(paste0(dir_runs, s_c_log), head = T))
s_c_mle <- get_mle(paste0(dir_runs, s_c_mle), beast_path)

cholera_results[1, ] <- c(s_c_diagnose[grep('clock[.]rate', rownames(s_c_diagnose)), 1:3], s_c_mle)

#strict_bs
s_bs_log <- grep('strict_bs', all_log_files[-c(mle_files, prior_files)], value = TRUE)
s_bs_mle <- grep('strict_bs', all_log_files[mle_files], value = TRUE)

s_bs_diagnose <- diagnose_mcmc(read.table(paste0(dir_runs, s_bs_log), head = T))
s_bs_mle <- get_mle(paste0(dir_runs, s_bs_mle), beast_path)

cholera_results[2, ] <- c(s_bs_diagnose[grep('clock[.]rate', rownames(s_bs_diagnose)), 1:3], s_bs_mle)

#ucld_constant
uc_c_log <- grep('ucld_constant', all_log_files[-c(mle_files, prior_files)], value = TRUE)
uc_c_mle <- grep('ucld_constant', all_log_files[mle_files], value = TRUE)

uc_c_diagnose <- diagnose_mcmc(read.table(paste0(dir_runs, uc_c_log), head = T))
uc_c_mle <- get_mle(paste0(dir_runs, uc_c_mle), beast_path)

cholera_results[3, ] <- c(uc_c_diagnose[grep('ucld[.]mean', rownames(uc_c_diagnose)), 1:3], uc_c_mle)

#ucld_bs
uc_bs_log <- grep('ucld_bs', all_log_files[-c(mle_files, prior_files)], value = TRUE)
uc_bs_mle <- grep('ucld_bs', all_log_files[mle_files], value = TRUE)

uc_bs_diagnose <- diagnose_mcmc(read.table(paste0(dir_runs, uc_bs_log), head = T))
uc_bs_mle <- get_mle(paste0(dir_runs, uc_bs_mle), beast_path)

cholera_results[4, ] <- c(uc_bs_diagnose[grep('ucld[.]mean', rownames(uc_bs_diagnose)), 1:3], uc_bs_mle)


######
## PROCESS YERSINIA PESTIS ANALYSES
######


##############
#Set up files and output matrix
dir_runs <- 'bacteria_rates/y_pestis/runs/'
all_log_files <- grep('log', dir(dir_runs), value = T)
mle_files <- grep('mle', all_log_files)
prior_files <- grep('prior', all_log_files)

ypestis_results <- matrix(NA, 4, 5)
rownames(ypestis_results) <- c('strict_constant', 'strict_bs', 'ucld_constant', 'ucld_bs')
colnames(ypestis_results) <- c('Rate', 'Lower95', 'Upper95', 'MLE_path', 'MLE_ss')
ypestis_results <- as.data.frame(ypestis_results)

#strict_constant
s_c_log <- grep('strict_constant', all_log_files[-c(mle_files, prior_files)], value = TRUE)
s_c_mle <- grep('strict_constant', all_log_files[mle_files], value = TRUE)

s_c_diagnose <- diagnose_mcmc(read.table(paste0(dir_runs, s_c_log), head = T))
s_c_mle <- get_mle(paste0(dir_runs, s_c_mle), beast_path)

ypestis_results[1, ] <- c(s_c_diagnose[grep('clock[.]rate', rownames(s_c_diagnose)), 1:3], s_c_mle)

#strict_bs
s_bs_log <- grep('strict_bs', all_log_files[-c(mle_files, prior_files)], value = TRUE)
s_bs_mle <- grep('strict_bs', all_log_files[mle_files], value = TRUE)

s_bs_diagnose <- diagnose_mcmc(read.table(paste0(dir_runs, s_bs_log), head = T))
s_bs_mle <- get_mle(paste0(dir_runs, s_bs_mle), beast_path)

ypestis_results[2, ] <- c(s_bs_diagnose[grep('clock[.]rate', rownames(s_bs_diagnose)), 1:3], s_bs_mle)

#ucld_constant
uc_c_log <- grep('ucln_constant', all_log_files[-c(mle_files, prior_files)], value = TRUE)
uc_c_mle <- grep('ucln_constant', all_log_files[mle_files], value = TRUE)

uc_c_diagnose <- diagnose_mcmc(read.table(paste0(dir_runs, uc_c_log), head = T))
uc_c_mle <- get_mle(paste0(dir_runs, uc_c_mle), beast_path)

ypestis_results[3, ] <- c(uc_c_diagnose[grep('ucld[.]mean', rownames(uc_c_diagnose)), 1:3], uc_c_mle)

#ucld_bs
uc_bs_log <- grep('ucln_bs', all_log_files[-c(mle_files, prior_files)], value = TRUE)
uc_bs_mle <- grep('ucln_bs', all_log_files[mle_files], value = TRUE)

uc_bs_diagnose <- diagnose_mcmc(read.table(paste0(dir_runs, uc_bs_log), head = T))
uc_bs_mle <- get_mle(paste0(dir_runs, uc_bs_mle), beast_path)

ypestis_results[4, ] <- c(uc_bs_diagnose[grep('ucld[.]mean', rownames(uc_bs_diagnose)), 1:3], uc_bs_mle)


######
## PROCESS S.SONNEI ANALYSES
######

##############
#Set up files and output matrix
dir_runs <- 'bacteria_rates/shigella/runs/'
all_log_files <- grep('log', dir(dir_runs), value = T)
mle_files <- grep('mle', all_log_files)
prior_files <- grep('prior', all_log_files)

shigella_results <- matrix(NA, 4, 5)
rownames(shigella_results) <- c('strict_constant', 'strict_bs', 'ucld_constant', 'ucld_bs')
colnames(shigella_results) <- c('Rate', 'Lower95', 'Upper95', 'MLE_path', 'MLE_ss')
shigella_results <- as.data.frame(shigella_results)

#strict_constant
s_c_log <- grep('strict_constant', all_log_files[-c(mle_files, prior_files)], value = TRUE)
s_c_mle <- grep('strict_constant', all_log_files[mle_files], value = TRUE)

s_c_diagnose <- diagnose_mcmc(read.table(paste0(dir_runs, s_c_log), head = T))
s_c_mle <- get_mle(paste0(dir_runs, s_c_mle), beast_path)

shigella_results[1, ] <- c(s_c_diagnose[grep('clock[.]rate', rownames(s_c_diagnose)), 1:3], s_c_mle)

#strict_bs
s_bs_log <- grep('strict_bs', all_log_files[-c(mle_files, prior_files)], value = TRUE)
s_bs_mle <- grep('strict_bs', all_log_files[mle_files], value = TRUE)

s_bs_diagnose <- diagnose_mcmc(read.table(paste0(dir_runs, s_bs_log), head = T))
s_bs_mle <- get_mle(paste0(dir_runs, s_bs_mle), beast_path)

shigella_results[2, ] <- c(s_bs_diagnose[grep('clock[.]rate', rownames(s_bs_diagnose)), 1:3], s_bs_mle)

#ucld_constant
uc_c_log <- grep('ucld_constant', all_log_files[-c(mle_files, prior_files)], value = TRUE)
uc_c_mle <- grep('ucld_constant', all_log_files[mle_files], value = TRUE)

uc_c_diagnose <- diagnose_mcmc(read.table(paste0(dir_runs, uc_c_log), head = T))
uc_c_mle <- get_mle(paste0(dir_runs, uc_c_mle), beast_path)

shigella_results[3, ] <- c(uc_c_diagnose[grep('ucld[.]mean', rownames(uc_c_diagnose)), 1:3], uc_c_mle)

#ucld_bs
uc_bs_log <- grep('ucld_bs', all_log_files[-c(mle_files, prior_files)], value = TRUE)
uc_bs_mle <- grep('ucld_bs', all_log_files[mle_files], value = TRUE)

uc_bs_diagnose <- diagnose_mcmc(read.table(paste0(dir_runs, uc_bs_log), head = T))
uc_bs_mle <- get_mle(paste0(dir_runs, uc_bs_mle), beast_path)

shigella_results[4, ] <- c(uc_bs_diagnose[grep('ucld[.]mean', rownames(uc_bs_diagnose)), 1:3], uc_bs_mle)


#######
## PROCESS M. LEPRAE
######

######
# Set up files and output matrix
dir_runs <- 'bacteria_rates/M_leprae/runs/'
all_log_files <- grep('log', dir(dir_runs), value = T)
mle_files <- grep('[.]mle[.]', all_log_files)
prior_files <- grep('prior', all_log_files)

mleprae_results <- matrix(NA, 4, 5)
rownames(mleprae_results) <- c('strict_constant', 'strict_bs', 'ucld_constant', 'ucld_bs')
colnames(mleprae_results) <- c('Rate', 'Lower95', 'Upper95', 'MLE_path', 'MLE_ss')
mleprae_results <- as.data.frame(mleprae_results)

#strict_constant
s_c_log <- grep('strict_constant', all_log_files[-c(mle_files, prior_files)], value = TRUE)
s_c_mle <- grep('strict_constant', all_log_files[mle_files], value = TRUE)

s_c_diagnose <- diagnose_mcmc(read.table(paste0(dir_runs, s_c_log), head = T))
s_c_mle <- get_mle(paste0(dir_runs, s_c_mle), beast_path)

mleprae_results[1, ] <- c(s_c_diagnose[grep('clock[.]rate', rownames(s_c_diagnose)), 1:3], s_c_mle)

#strict_bs
s_bs_log <- grep('strict_bs', all_log_files[-c(mle_files, prior_files)], value = TRUE)
s_bs_mle <- grep('strict_bs', all_log_files[mle_files], value = TRUE)

s_bs_diagnose <- diagnose_mcmc(read.table(paste0(dir_runs, s_bs_log), head = T))
s_bs_mle <- get_mle(paste0(dir_runs, s_bs_mle), beast_path)

mleprae_results[2, ] <- c(s_bs_diagnose[grep('clock[.]rate', rownames(s_bs_diagnose)), 1:3], s_bs_mle)

#ucld_constant
uc_c_log <- grep('ucld_constant', all_log_files[-c(mle_files, prior_files)], value = TRUE)
uc_c_mle <- grep('ucld_constant', all_log_files[mle_files], value = TRUE)

uc_c_diagnose <- diagnose_mcmc(read.table(paste0(dir_runs, uc_c_log), head = T))
uc_c_mle <- get_mle(paste0(dir_runs, uc_c_mle), beast_path)

mleprae_results[3, ] <- c(uc_c_diagnose[grep('ucld[.]mean', rownames(uc_c_diagnose)), 1:3], uc_c_mle)

#ucld_bs
uc_bs_log <- grep('ucld_bs', all_log_files[-c(mle_files, prior_files)], value = TRUE)
uc_bs_mle <- grep('ucld_bs', all_log_files[mle_files], value = TRUE)

uc_bs_diagnose <- diagnose_mcmc(read.table(paste0(dir_runs, uc_bs_log), head = T))
uc_bs_mle <- get_mle(paste0(dir_runs, uc_bs_mle), beast_path)

mleprae_results[4, ] <- c(uc_bs_diagnose[grep('ucld[.]mean', rownames(uc_bs_diagnose)), 1:3], uc_bs_mle)


