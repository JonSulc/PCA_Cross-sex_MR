require( data.table )
require( tidyverse )
require( TwoSampleMR )
require( ieugwasr )
require( parallel )
source( 'ld_clump.R' )

THRESHOLD       =  5e-8
RANK_THRESHOLD  =  5e-8
PROCESSED_PATH  =  'data/processed/%s/%s/%s.rds'
MIN_IVS         =  10

MAX_TRIES  =  50
CLUMP_KB   =  1e4
CLUMP_R2   =  0.01
CLUMP_P    =  1

SNP_CHUNK_SIZE = 3e4 # Required for clumping larger numbers of SNPs online

if (!exists( 'COMMON_SNPS' ))
    COMMON_SNPS  =  read_rds( 'data/common_snps.rds' )

PLINK_PATH  =  genetics.binaRies::get_plink_binary()
LD_PATH  =  'data/uk10k/uk10k_full'
# VCFFILE_PATTERN  =  '/data/sgg2/data/UK10K/for_SSimp/uk10k.chr%i.vcf.gz'

get_variables_from_all  =  function( category_name,
                                     all_phenotypes,
                                     variables      = c( 'beta', 'se', 'pval' ),
                                     snps_to_keep   = NULL,
                                     variant_list   = NULL,
                                     sex            = c( 'both_sexes', 'female', 'male' ),
                                     processed_path = PROCESSED_PATH,
                                     common_snps    = COMMON_SNPS,
                                     threshold      = THRESHOLD,
                                     return_names   = length( variables ) < 2,
                                     ... ) {
    sex = match.arg( sex )
    print( category_name )
    
    if (is.null( snps_to_keep ) & !is.null( variant_list )) {
        snps_to_keep  =  common_snps$variant %chin% variant_list
    }
    
    if (length( category_name ) > 1) {
        # category_stats  =  mclapply( category_name,
        category_stats  =  lapply( category_name,
                  get_variables_from_all,
                  all_phenotypes = all_phenotypes,
                  variables      = variables,
                  snps_to_keep   = snps_to_keep,
                  sex            = sex,
                  processed_path = processed_path,
                  common_snps    = common_snps,
                  threshold      = threshold,
                  ... )
        return(  setNames( nm = variables ) %>%
                     lapply( function( variable ){
                        lapply( category_stats, function( x ) x[[ variable ]] ) %>%
                            reduce( full_join )
                    } )
        )
    }
    
    phenotypes  =  all_phenotypes %>%
        filter( category == category_name ) %>%
        '$'( 'phenotype' )
    
    phenotypes  =  phenotypes %>%
        '['( sprintf( processed_path, sex, category_name, . ) %>% file.exists )
    
    stats  =  phenotypes %>%
        sprintf( processed_path,
                 sex,
                 category_name,
                 . ) %>%
        # mclapply( .get_variables_from_single_file,
        lapply( .get_variables_from_single_file,
                  variables    = variables,
                  snps_to_keep = snps_to_keep,
                  all_variants = common_snps$variant,
                  threshold    = threshold,
                  return_names = return_names )
    
    if (return_names) {
        return( unlist( stats ) %>%
                    unique )
    }
    if (!is.null( snps_to_keep )) {
        common_snps = common_snps[ snps_to_keep, ]
    }
    
    setNames( nm = variables) %>%
        lapply( function( variable ){
            lapply( stats, function( x ) x[[ variable ]] ) %>%
                bind_cols( common_snps[ , c( 'variant', 'rsid' ) ], . )
        } )
}

get_variables_from_each  =  function( category_name,
                                      all_phenotypes,
                                      variables      = c( 'beta', 'se', 'pval' ),
                                      snps_to_keep   = NULL,
                                      variant_list   = NULL,
                                      sex            = c( 'both_sexes', 'female', 'male' ),
                                      processed_path = PROCESSED_PATH,
                                      common_snps    = COMMON_SNPS,
                                      threshold      = THRESHOLD,
                                      return_names   = length( variables ) < 2,
                                      clump          = TRUE,
                                      ... ) {
    sex = match.arg( sex )
    
    if (is.null( snps_to_keep ) & !is.null( variant_list )) {
        snps_to_keep  =  common_snps$variant %chin% variant_list
    }
    
    if (length( category_name ) > 1) {
        # category_stats  =  mclapply( category_name,
        category_stats  =  lapply( category_name,
                                   get_variables_from_each,
                                   all_phenotypes = all_phenotypes,
                                   variables      = variables,
                                   snps_to_keep   = snps_to_keep,
                                   variant_list   = NULL,
                                   sex            = sex,
                                   processed_path = processed_path,
                                   common_snps    = common_snps,
                                   threshold      = threshold,
                                   return_names   = return_names,
                                   ... )
        return(  category_stats )
    }
    
    phenotypes  =  all_phenotypes %>%
        filter( category == category_name ) %>%
        '$'( 'phenotype' ) %>%
        '['( sprintf( processed_path, sex, category_name, . ) %>% file.exists )
    
    stats  =  phenotypes %>%
        sprintf( processed_path,
                 sex,
                 category_name,
                 . ) %>%
        setNames( phenotypes ) %>%
        # mclapply( .get_variables_from_single_file,
        lapply( .get_variables_from_single_file,
                  variables    = variables,
                  snps_to_keep = snps_to_keep,
                  all_variants = common_snps$variant,
                  threshold    = threshold,
                  return_names = return_names,
                  clump        = clump )
    return( stats )
}

get_clumped_from_each  =  function( clumped_snps_path,
                                    processed_path,
                                    common_snps = COMMON_SNPS,
                                    threshold = 1 ) {
    phenotype  =  str_match( clumped_snps_path,
                             'clump/(/)?(.+)[.]gwas' )[ , 3 ]
    if (length(clumped_snps_path) > 1) {
        return( clumped_snps_path %>%
                    setNames( phenotype ) %>%
                    lapply( get_clumped_from_each,
                            processed_path = processed_path,
                            common_snps = common_snps,
                            threshold = threshold ) )
    }
    all_snp_names = list.files( clumped_snps_path,
                                pattern = '.+[.]clumped$',
                                full.names = TRUE ) %>%
        lapply( function( filename ) {
            read_delim( filename,
                        delim = ' ',
                        trim_ws = TRUE )$SNP
        } ) %>%
        reduce( c )
    snps_to_keep  =  common_snps$rsid %chin% all_snp_names
    
    sprintf( processed_path,
             phenotype ) %>%
        .get_variables_from_single_file( snps_to_keep = snps_to_keep,
                                         threshold = threshold ) %>%
        bind_cols
}

clump_each  =  function( all_snps_summary_stats,
                         threshold = THRESHOLD,
                         min_ivs = MIN_IVS,
                         return_names = FALSE,
                         ... ) {
    if (threshold & threshold < 1) {
        min_pvals  =  all_snps_summary_stats$pval %>%
            select_if( is.numeric ) %>%
            apply( 1, min, na.rm = TRUE )
        to_keep  =  min_pvals < threshold
        all_snps_summary_stats[ c( 'beta', 'se', 'pval' ) ]  =  all_snps_summary_stats[ c( 'beta', 'se', 'pval' ) ] %>%
            lapply( function( x ) x[ to_keep, ] )
    }
    clumped_snps  =  all_snps_summary_stats$pval %>%
        select_if( is.numeric ) %>%
        colnames %>%
        setNames( nm = . ) %>%
        lapply( function( name ){
            all_snps_summary_stats$pval %>%
                filter( get(name) < threshold ) %>%
                (function( x ){
                    if (nrow(x) >= min_ivs) {
                        tibble( rsid      = x$rsid,
                                variant   = x$variant,
                                pval      = x[ , name ],
                                phenotype = name,
                                chr       = str_match( variant, '^[0-9X]+' )[ , 1 ],
                                position  = str_match( variant, '^[0-9X]+:([0-9]+):' )[ , 2 ] )
                    } else {
                        tibble( rsid      = character(),
                                variant   = character(),
                                pval      = numeric(),
                                phenotype = character(),
                                chr       = numeric(),
                                position  = numeric() )
                    }
                })
        } ) %>%
        lapply( clump_genome, ... )
    if (return_names) {
        return( clumped_snps )
    }
    
    if (!is.null( min_ivs )) {
        to_keep  =  sapply( clumped_snps, length ) >= min_ivs
        clumped_snps  =  clumped_snps[ to_keep ]
    }
    clumped_stats  =  clumped_snps %>%
        names %>%
        setNames( nm = . ) %>%
        lapply( function( pc_name ){
            lapply( c( 'beta', 'se', 'pval' ), function( stat ){
                all_snps_summary_stats[[ stat ]] %>%
                    select( variant, rsid, !!stat := !!pc_name ) %>%
                    left_join( tibble( variant = clumped_snps[[ pc_name ]] ), . )
            } ) %>%
                reduce( full_join )
        } )
    
    if ('pca' %in% names( all_snps_summary_stats )) {
        clumped_stats$pca  =  all_snps_summary_stats$pca[ c( 'rotation', 'correlation_matrix' ) ]
        if (!is.null( min_ivs )) {
            clumped_stats$pca$rotation  =  all_snps_summary_stats$pca$rotation[ , to_keep, drop = FALSE ]
        }
    }
    clumped_stats
}

.get_variables_from_single_file  =  function( filename,
                                              variables    = c( 'beta', 'se', 'pval' ),
                                              snps_to_keep = NULL,
                                              common_snps  = COMMON_SNPS,
                                              threshold    = THRESHOLD,
                                              return_names = length( variables ) < 2,
                                              clump        = FALSE,
                                              ... ) {
    if (!file.exists( filename )) return()
    stats = read_rds( filename ) %>%
        mutate( pval = 2 * pnorm( -abs(beta / se) ) ) %>%
        select( !!variables )
    
    if (!is.null( snps_to_keep )){
        stats  =  filter( stats, snps_to_keep )
        if (threshold & threshold < 1) {
            stats  =  filter( stats, pval < threshold )
        }
    } else if (clump) {
        return( common_snps %>%
            cbind( stats ) %>%
            filter( pval < threshold ) %>%
            clump_genome( ... ) %>%
            tibble( variant = . ) %>%
            left_join( cbind( common_snps[ , c( 'variant', 'rsid' ) ],
                              stats ) ) )
    }
    if (return_names) {
        return( common_snps$variant[ stats$pval < threshold & !is.na( stats$pval ) ] )
    }
    
    phenotype = str_match( filename, '([^/]+)[.]rds' )[ , 2 ]
    stats %>%
        lapply( function( x ) tibble( !!phenotype := x ) )
}

get_clumped_from_all  =  function( category_name,
                                   all_phenotypes,
                                   variables      = c( 'beta', 'se', 'pval' ),
                                   threshold      = THRESHOLD,
                                   ... ){
    get_variables_from_all( category_name,
                            all_phenotypes = all_phenotypes,
                            variables      = 'pval',
                            snps_to_keep   = NULL,
                            ... ) %>%
        get_variables_from_all( category_name,
                                all_phenotypes = all_phenotypes,
                                variables      = 'pval',
                                snps_to_keep   = .,
                                threshold      = FALSE,
                                ... ) %>%
        get_ranks( ... ) %>%
        clump_genome( ... ) %>%
        get_variables_from_all( category_name,
                                all_phenotypes = all_phenotypes,
                                variables      = variables,
                                variant_list   = .,
                                threshold      = FALSE,
                                ... )
}

get_ranks  =  function( pvalues,
                        phenotype = 'phenotype',
                        rank_threshold = RANK_THRESHOLD,
                        ... ){
    pvalues %>%
        transmute_if( is.numeric, function( x ) replace( x, x > rank_threshold, NA ) ) %>%
        transmute_all( rank, na.last = 'keep' ) %>%
        apply( 1, min, na.rm = TRUE ) %>%
        '/'( max( ., na.rm = TRUE ) ) %>%
        tibble( pval = . ) %>%
        bind_cols( pvalues %>% select_if( function(x) !is.numeric(x) ) ) %>%
        mutate( phenotype = phenotype,
                chr       = str_match( variant, '^[0-9X]+' )[ , 1 ],
                position  = str_match( variant, '^[0-9]+:([0-9]+):' )[ , 2 ] )
}

snp_bin  =  function( snp_ranks,
                      chunk_size = SNP_CHUNK_SIZE ){
    if (nrow( snp_ranks ) == 0) {
        return()
    }
    
    max_chr  =  snp_ranks$chr %>%
        table %>%
        cumsum %>%
        (function(x) x < chunk_size) %>%
        (function(x) names(x)[ max(which(x)) ] ) %>%
        as.numeric
    if (is.na( max_chr )) {
        max_chr = min( snp_ranks$chr )
    }
    
    bin = snp_ranks %>%
                filter( chr <= max_chr ) %>%
                list
    return( c( bin,
               snp_bin( snp_ranks[ snp_ranks$chr > max_chr, ],
                        chunk_size ) ) )
}

clump_genome  =  function( snp_ranks,
                           chunk_size = SNP_CHUNK_SIZE,
                           ... ){
    if (is.null( chunk_size )) {
        .clump_chromosome( snp_ranks, ... ) %>%
            ( function( rsids ) snp_ranks$variant[ snp_ranks$rsid %chin% rsids ] )
    } else {
        snp_ranks %>%
            snp_bin( chunk_size ) %>%
            lapply( .clump_chromosome, ... ) %>%
            unlist %>%
            unname %>%
            ( function( rsids ) snp_ranks$variant[ snp_ranks$rsid %chin% rsids ] )
    }
}

.clump_chromosome  =  function( chr_snp_ranks,
                                current_try = 1,
                                max_tries = MAX_TRIES,
                                clump_kb = CLUMP_KB,
                                clump_r2 = CLUMP_R2,
                                clump_p = CLUMP_P,
                                plink_path = PLINK_PATH,
                                ld_path = LD_PATH ){
    if (nrow( chr_snp_ranks ) == 0){
        print( paste( 'Skipping empty chromosome' ) )
        return()
    }
    # if (use_uk10k) {
        return(
            chr_snp_ranks %>%
                as_tibble %>%
                select( rsid, pval ) %>%
                ld_clump( plink_bin = plink_path,
                          bfile = ld_path,
                          clump_kb = clump_kb,
                          clump_r2 = clump_r2,
                          clump_p = clump_p ) %>%
                '$'( 'rsid' )
        )
    # }
    # tryCatch({
    #     # print( paste( 'Submitting chromosome', chr_snp_ranks$chr[1] ) )
    #     return(
    #         chr_snp_ranks %>%
    #             as_tibble %>%
    #             select( rsid, pval, id = chr ) %>%
    #             ld_clump( plink_bin = plink_path,
    #                       bfile = ld_path,
    #                       clump_kb = clump_kb,
    #                       clump_r2 = clump_r2,
    #                       clump_p = clump_p ) %>%
    #             '$'( 'rsid' )
    #         # chr_snp_ranks %>%
    #         #     format_data( snp_col = 'rsid',
    #         #                  phenotype_col = 'phenotype',
    #         #                  pval_col = 'pval' ) %>%
    #         #     clump_data( clump_kb = clump_kb,
    #         #                 clump_r2 = clump_r2 ) %>%
    #         #     '$'( 'SNP' )
    #     )
    # }, error = function( e ){
    #     print( "Warning: none of the SNPs were present in the reference panel." )
    #     return()
    # #     return()
    # # #     print( paste( 'Failed try', current_try ) )
    # # #     if (current_try < max_tries){
    # # #         return( .clump_chromosome( chr_snp_ranks,
    # # #                                    current_try = current_try + 1,
    # # #                                    max_tries = max_tries,
    # # #                                    clump_kb = clump_kb,
    # # #                                    clump_r2 = clump_r2 ) )
    # # #     } else {
    # # #         stop( paste( 'Failed, giving up on chromosome', chr_snp_ranks$chr[1] ) )
    # # #         # failed_snps <<- chr_snp_ranks
    # # #         # return()
    # # #     }
    # })
}













