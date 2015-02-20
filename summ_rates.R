library(coda)
library(RSQLite)
source('bact_rates_scripts/proc_data.R')

con <- dbConnect(RSQLite::SQLite(), 'log_summary_bacteria.db')

## This code should be repeated for each data set
all_comp_file_names <- grep('mleprae.+log', dir('bacteria_rates_allruns'), value = T)
mle_loc <- grep('[.]mle[.]', all_comp_file_names)

mle_files <- paste0('bacteria_rates_allruns/', all_comp_file_names[mle_loc])
log_files <- paste0('bacteria_rates_allruns/', all_comp_file_names[-mle_loc])

list_logs <- list()
for(i in 1:length(log_files)){
      log_temp <- read.table(log_files[i], head = T)
      list_logs[[i]] <- diagnose_mcmc(log_temp)
      rownames(list_logs[[i]]) <- gsub('[.]', '_', rownames(list_logs[[i]]))
      colnames(list_logs[[i]]) <- gsub('%', '_', colnames(list_logs[[i]]))
}
names(list_logs) <- gsub('^.+/|[.].+$', '', log_files)

for(i in 1:length(list_logs)){
      dbWriteTable(conn = con, name = names(list_logs)[i], value = data.frame(list_logs[[i]]), overwrite = T)
}




list_mles <- list()
for(i in 1:length(mle_files)){
      list_mles[[i]] <- get_mle(mle_files[i], 'beast')
}

matrix_mles <- concat_list(list_mles)
rownames(matrix_mles) <- gsub('^.+/', '', mle_files)

dbWriteTable(conn = con, name = 'mlep_mles', value = data.frame(matrix_mles), overwrite = T)
dbDisconnect(con)
