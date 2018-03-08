setwd("~/Dropbox/glasso/concord_results/")

out_file = "chosen_incov.csv"
for (subdir in list.dirs(recursive = FALSE)) {
	for (cov_file in list.files(subdir, pattern = "lambda0.0002_nonzero*")){
		load(paste(subdir, cov_file, sep = "/"))
		write.csv(inv_cov, file = paste(subdir, out_file, sep = "/"))
	}	
}

