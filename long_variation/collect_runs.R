##### Collect per-task results into the same objects the local run produced
# Each array task wrote sim_results/run_XXXX.rds containing the list returned
# by run_one_iter(): elements acc, cm, params, empiricals.

files <- list.files("sim_results", pattern = "^run_[0-9]+\\.rds$", full.names = TRUE)

if (length(files) == 0) stop("No sim_results/run_*.rds files found.")

# read in numeric task order (run_0001, run_0002, ...) so rows line up with iterations
ids   <- as.integer(gsub("^run_0*([0-9]+)\\.rds$", "\\1", basename(files)))
files <- files[order(ids)]
ids   <- sort(ids)

cat("Found", length(files), "result files (iterations",
    min(ids), "to", max(ids), ").\n")

# warn about any gaps in the expected 1..max range (a task that died / didn't run)
missing <- setdiff(seq_len(max(ids)), ids)
if (length(missing) > 0) {
  warning("Missing iterations: ", paste(missing, collapse = ", "))
}

results_list <- lapply(files, readRDS)

### reshape results  (identical to the local script)
output            <- as.data.frame(do.call(rbind, lapply(results_list, `[[`, "acc")))
cm_results        <- lapply(results_list, `[[`, "cm")
param_results     <- lapply(results_list, `[[`, "params")
empirical_results <- lapply(results_list, `[[`, "empiricals")

acc <- colMeans(output)
print(acc)

write.csv(output, "output_fix_long.csv")
saveRDS(cm_results, "cm_results_fix_long.rds")
saveRDS(param_results, "param_results_fix_long.rds")

cat("Wrote output_fix_long.csv, cm_results_fix_long.rds, param_results_fix_long.rds\n")
