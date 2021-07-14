ADDITIONAL_FLAGS  =  ' --exclude data/uk10k/exclude.txt'

ld_clump_vcf  =  function (dat = NULL, clump_kb = 10000, clump_r2 = 0.001, clump_p = 0.99, 
          pop = "EUR", access_token = NULL,
          vcffile_pattern = NULL, plink_bin = NULL, ...) 
{
    stopifnot("rsid" %in% names(dat))
    stopifnot(is.data.frame(dat))
    if (!"pval" %in% names(dat)) {
        if ("p" %in% names(dat)) {
            warning("No 'pval' column found in dat object. Using 'p' column.")
            dat[["pval"]] <- dat[["p"]]
        }
        else {
            warning("No 'pval' column found in dat object. Setting p-values for all SNPs to clump_p parameter.")
            dat[["pval"]] <- clump_p
        }
    }
    if (!"id" %in% names(dat)) {
        dat$id <- random_string(1)
    }
    if (is.null(vcffile_pattern)) {
        access_token = check_access_token()
    }
    ids <- unique(dat[["id"]])
    res <- list()
    for (i in 1:length(ids)) {
        x <- subset(dat, dat[["id"]] == ids[i])
        if (nrow(x) == 1) {
            message("Only one SNP for ", ids[i])
            res[[i]] <- x
        }
        else {
            message("Clumping ", ids[i], ", ", nrow(x), " variants, using ", 
                    pop, " population reference")
            vcffile  =  sprintf( vcffile_pattern, ids[i] )
            res[[i]] <- ld_clump_local_vcf(x, clump_kb = clump_kb, 
                                           clump_r2 = clump_r2,
                                           clump_p = clump_p,
                                           vcffile = sprintf( vcffile_pattern,
                                                              as.numeric(ids[i]) ),
                                           plink_bin = plink_bin,
                                           ... )
        }
    }
    res <- dplyr::bind_rows(res)
    return(res)
}



ld_clump_local_vcf  =  function (dat, clump_kb, clump_r2, clump_p, vcffile, plink_bin,
                                 additional_flags = ADDITIONAL_FLAGS ) 
{
    shell <- ifelse(Sys.info()["sysname"] == "Windows", "cmd", 
                    "sh")
    fn <- tempfile()
    write.table(data.frame(SNP = dat[["rsid"]], P = dat[["pval"]]), 
                file = fn, row.names = F, col.names = T, quote = F)
    fun2 <- paste0(shQuote(plink_bin, type = shell),
                   " --vcf ", shQuote(vcffile, type = shell),
                   " --clump ", shQuote(fn, type = shell),
                   " --clump-p1 ", clump_p,
                   " --clump-r2 ", clump_r2,
                   " --clump-kb ", clump_kb,
                   additional_flags,
                   " --out ", shQuote(fn, type = shell))
    print(fun2)
    system(fun2)
    res <- read.table(paste(fn, ".clumped", sep = ""), header = T)
    unlink(paste(fn, "*", sep = ""))
    y <- subset(dat, !dat[["rsid"]] %in% res[["SNP"]])
    if (nrow(y) > 0) {
        message("Removing ", length(y[["rsid"]]), " of ", nrow(dat), 
                " variants due to LD with other variants or absence from LD reference panel")
    }
    return(subset(dat, dat[["rsid"]] %in% res[["SNP"]]))
}

