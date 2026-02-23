#call this script from the root of the project with:
library(sleuth)

# Define metadata: 0030/0044 are 2dpi, 0033/0045 are 6dpi
#this will need to be updated if you change the sample names or add more samples
#i did this because sleuth needs a specific format for the metadata, and this is a simple way to create it
sample_id <- c("SRR5660030", "SRR5660033", "SRR5660044", "SRR5660045")
condition <- c("2dpi", "6dpi", "2dpi", "6dpi")
path <- paste0("results/", sample_id, "_kallisto")

#this was done manually, but you could also automate it by listing the files in the results directory and extracting the sample names
s2c <- data.frame(sample=sample_id, condition=condition, path=path, stringsAsFactors=F)

# Run Sleuth
#i was able to do this by following the sleuth tutorial, but you may need to adjust the model formula if you have more complex experimental designs
so <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE)
so <- sleuth_fit(so, ~condition, 'full')
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')

# Extract results
#thought process: i want to get the results of the likelihood ratio test, which compares the full model (with condition) to the reduced model 
# (without condition). this will tell me if there are any genes that are differentially expressed between the conditions. 
# i will then filter for significant hits and write them to a file in a format that can be easily read by downstream analysis tools.
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
significant <- sleuth_table[sleuth_table$qval < 0.05, ]

# Write significant hits to report format
write.table(significant[, c("target_id", "test_stat", "pval", "qval")], 
            file="results/sleuth_results.txt", sep="\t", quote=F, row.names=F)