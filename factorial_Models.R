# ======================
# ---- main program ----  
source('bayesian_updates.R')

# prepare possible models 
models <- c("L", "D")
resp_update_pars <- c("R0", "R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9", "RA", "RB", "RC", "RD", "RE")
color_update_pars <- c("C0")
position_update_pars <- c("P0")

modelnames <- apply(expand.grid(apply(expand.grid(apply(expand.grid(models, resp_update_pars), 1, paste, collapse=""),
                                                  color_update_pars), 1, paste, collapse=""), position_update_pars), 1, paste, collapse="")

subfiles=dir('data', '*.dat')
data <- data.table(do.call('rbind',lapply(subfiles, readData)))
data <- unite(data, "sub_session", sub, session)

lsubs <- split(data, data$sub_session) # split data for each subject and session

# Now parallel optimization
# detection cores and distribution computing among cores
g <- expand.grid(dat = lsubs, mod = modelnames, stringsAsFactors = FALSE)
lpars <- parEstimate(g$dat, g$mod)

# now convert lists to data.frame
pars.df <- as.data.frame(do.call(rbind, lpars), stringsAsFactors = FALSE)
names(pars.df)[35:36] <- c("AIC", "BIC")
row.names(pars.df) <- apply(expand.grid(unique(pars.df$sub), modelnames), 1, paste, collapse="_")
# convert character columsn to numeric
pars.df[,1:36] <- sapply(pars.df[,1:36], as.numeric)
pars.df <- separate(pars.df, sub, c("sub", "session"))

# save it.
saveRDS(pars.df, file = './data/fitpars_resp.rds')


resp_update_pars <- c("R0")
color_update_pars <- c("C0", "C1", "C2", "C3", "C4", "C5", "C6", "C7")
position_update_pars <- c("P0")

modelnames <- apply(expand.grid(apply(expand.grid(apply(expand.grid(models, resp_update_pars), 1, paste, collapse=""),
                                                  color_update_pars), 1, paste, collapse=""), position_update_pars), 1, paste, collapse="")

# Now parallel optimization
# detection cores and distribution computing among cores
g <- expand.grid(dat = lsubs, mod = modelnames, stringsAsFactors = FALSE)
lpars <- parEstimate(g$dat, g$mod)

# now convert lists to data.frame
pars.df <- as.data.frame(do.call(rbind, lpars), stringsAsFactors = FALSE)
names(pars.df)[35:36] <- c("AIC", "BIC")
row.names(pars.df) <- apply(expand.grid(unique(pars.df$sub), modelnames), 1, paste, collapse="_")
# convert character columsn to numeric
pars.df[,1:36] <- sapply(pars.df[,1:36], as.numeric)
pars.df <- separate(pars.df, sub, c("sub", "session"))

# save it.
saveRDS(pars.df, file = './data/fitpars_col.rds')


resp_update_pars <- c("R0")
color_update_pars <- c("C0")
position_update_pars <- c("P0", "P1", "P2", "P3", "P4", "P5", "P6", 'P7', 'P8', 'P9')

modelnames <- apply(expand.grid(apply(expand.grid(apply(expand.grid(models, resp_update_pars), 1, paste, collapse=""),
                                                  color_update_pars), 1, paste, collapse=""), position_update_pars), 1, paste, collapse="")

# Now parallel optimization
# detection cores and distribution computing among cores
g <- expand.grid(dat = lsubs, mod = modelnames, stringsAsFactors = FALSE)
lpars <- parEstimate(g$dat, g$mod)

# now convert lists to data.frame
pars.df <- as.data.frame(do.call(rbind, lpars), stringsAsFactors = FALSE)
names(pars.df)[35:36] <- c("AIC", "BIC")
row.names(pars.df) <- apply(expand.grid(unique(pars.df$sub), modelnames), 1, paste, collapse="_")
# convert character columsn to numeric
pars.df[,1:36] <- sapply(pars.df[,1:36], as.numeric)
pars.df <- separate(pars.df, sub, c("sub", "session"))

# save it.
saveRDS(pars.df, file = './data/fitpars_pos.rds')


# Now run each comparison again with the other factor levels set to the winner for that factor
resp_update_pars <- c("R0", "R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9", "RA", "RB", "RC", "RD", "RE")
color_update_pars <- c("C3")
position_update_pars <- c("P7")

modelnames <- apply(expand.grid(apply(expand.grid(apply(expand.grid(models, resp_update_pars), 1, paste, collapse=""),
                                                  color_update_pars), 1, paste, collapse=""), position_update_pars), 1, paste, collapse="")

g <- expand.grid(dat = lsubs, mod = modelnames, stringsAsFactors = FALSE)
lpars <- parEstimate(g$dat, g$mod)

# now convert lists to data.frame
pars.df <- as.data.frame(do.call(rbind, lpars), stringsAsFactors = FALSE)
names(pars.df)[35:36] <- c("AIC", "BIC")
row.names(pars.df) <- apply(expand.grid(unique(pars.df$sub), modelnames), 1, paste, collapse="_")
# convert character columsn to numeric
pars.df[,1:36] <- sapply(pars.df[,1:36], as.numeric)
pars.df <- separate(pars.df, sub, c("sub", "session"))

# save it.
saveRDS(pars.df, file = './data/fitpars_opt_resp.rds')

# Now run each comparison again with the other factor levels set to the winner for that factor
resp_update_pars <- c("R4")
color_update_pars <- c("C0", "C1", "C2", "C3", "C4", "C5", "C6", "C7")
position_update_pars <- c("P7")

modelnames <- apply(expand.grid(apply(expand.grid(apply(expand.grid(models, resp_update_pars), 1, paste, collapse=""),
                                                  color_update_pars), 1, paste, collapse=""), position_update_pars), 1, paste, collapse="")

g <- expand.grid(dat = lsubs, mod = modelnames, stringsAsFactors = FALSE)
lpars <- parEstimate(g$dat, g$mod)

# now convert lists to data.frame
pars.df <- as.data.frame(do.call(rbind, lpars), stringsAsFactors = FALSE)
names(pars.df)[35:36] <- c("AIC", "BIC")
row.names(pars.df) <- apply(expand.grid(unique(pars.df$sub), modelnames), 1, paste, collapse="_")
# convert character columsn to numeric
pars.df[,1:36] <- sapply(pars.df[,1:36], as.numeric)
pars.df <- separate(pars.df, sub, c("sub", "session"))

# save it.
saveRDS(pars.df, file = './data/fitpars_opt_col.rds')

# Now run each comparison again with the other factor levels set to the winner for that factor
resp_update_pars <- c("R4")
color_update_pars <- c("C3")
position_update_pars <- c("P0", "P1", "P2", "P3", "P4", "P5", "P6", 'P7', 'P8', 'P9')

modelnames <- apply(expand.grid(apply(expand.grid(apply(expand.grid(models, resp_update_pars), 1, paste, collapse=""),
                                                  color_update_pars), 1, paste, collapse=""), position_update_pars), 1, paste, collapse="")

g <- expand.grid(dat = lsubs, mod = modelnames, stringsAsFactors = FALSE)
lpars <- parEstimate(g$dat, g$mod)

# now convert lists to data.frame
pars.df <- as.data.frame(do.call(rbind, lpars), stringsAsFactors = FALSE)
names(pars.df)[35:36] <- c("AIC", "BIC")
row.names(pars.df) <- apply(expand.grid(unique(pars.df$sub), modelnames), 1, paste, collapse="_")
# convert character columsn to numeric
pars.df[,1:36] <- sapply(pars.df[,1:36], as.numeric)
pars.df <- separate(pars.df, sub, c("sub", "session"))

# save it.
saveRDS(pars.df, file = './data/fitpars_opt_pos.rds')

# The best model for the position based updating changed when running with the optimal parameters for 
# response and color, so now we run the optimization for response and color again with the new best model 
# for position

# Now run each comparison again with the other factor levels set to the winner for that factor
resp_update_pars <- c("R0", "R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9", "RA", "RB", "RC", "RD", "RE")
color_update_pars <- c("C3")
position_update_pars <- c("P4")

modelnames <- apply(expand.grid(apply(expand.grid(apply(expand.grid(models, resp_update_pars), 1, paste, collapse=""),
                                                  color_update_pars), 1, paste, collapse=""), position_update_pars), 1, paste, collapse="")

g <- expand.grid(dat = lsubs, mod = modelnames, stringsAsFactors = FALSE)
lpars <- parEstimate(g$dat, g$mod)

# now convert lists to data.frame
pars.df <- as.data.frame(do.call(rbind, lpars), stringsAsFactors = FALSE)
names(pars.df)[35:36] <- c("AIC", "BIC")
row.names(pars.df) <- apply(expand.grid(unique(pars.df$sub), modelnames), 1, paste, collapse="_")
# convert character columsn to numeric
pars.df[,1:36] <- sapply(pars.df[,1:36], as.numeric)
pars.df <- separate(pars.df, sub, c("sub", "session"))

# save it.
saveRDS(pars.df, file = './data/fitpars_opt_r2_resp.rds')

# Now run each comparison again with the other factor levels set to the winner for that factor
resp_update_pars <- c("R4")
color_update_pars <- c("C0", "C1", "C2", "C3", "C4", "C5", "C6", "C7")
position_update_pars <- c("P4")

modelnames <- apply(expand.grid(apply(expand.grid(apply(expand.grid(models, resp_update_pars), 1, paste, collapse=""),
                                                  color_update_pars), 1, paste, collapse=""), position_update_pars), 1, paste, collapse="")

g <- expand.grid(dat = lsubs, mod = modelnames, stringsAsFactors = FALSE)
lpars <- parEstimate(g$dat, g$mod)

# now convert lists to data.frame
pars.df <- as.data.frame(do.call(rbind, lpars), stringsAsFactors = FALSE)
names(pars.df)[35:36] <- c("AIC", "BIC")
row.names(pars.df) <- apply(expand.grid(unique(pars.df$sub), modelnames), 1, paste, collapse="_")
# convert character columsn to numeric
pars.df[,1:36] <- sapply(pars.df[,1:36], as.numeric)
pars.df <- separate(pars.df, sub, c("sub", "session"))

# save it.
saveRDS(pars.df, file = './data/fitpars_opt_r2_col.rds')
