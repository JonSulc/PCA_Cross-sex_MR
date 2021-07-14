require( tidyverse,   quietly = TRUE, warn.conflicts = FALSE )
require( TwoSampleMR, quietly = TRUE )

rankNorm  =  function( x ){
    qnorm( (rank( x, na.last = "keep" ) - 0.5) / sum( !is.na( x ) ) )
}

source( 'all_snps_summary_stats.R' )

# Resolve conflicts
filter  =  dplyr::filter

# Default values
THRESHOLD  =  5e-8

COVARS  =  c( '21022' = 'age', # Age at recruitment
              '31'    = 'sex' )
# COVARS  =  c()

get_category_pca_stats  =  function( all_snps_summary_stats,
                                     all_phenotypes,
                                     category_name,
                                     category_snps,
                                     ... ){
    .select_all_category_stats( all_snps_summary_stats,
                                all_phenotypes,
                                category_name,
                                category_snps ) %>%
        get_pca_summary_stats( all_phenotypes = all_phenotypes,
                               category_name = category_name,
                               ... )
}

get_pca_summary_stats  =  function( all_snps_summary_stats,
                                    all_phenotypes,
                                    line_up_positive = FALSE,
                                    npcs = 10,
                                    scaling = FALSE,
                                    category_name = FALSE,
                                    correlation_matrix = NULL,
                                    sex = 'both_sexes',
                                    ... ) {
    if (is.null( correlation_matrix )) {
        correlation_matrix  =  cor( all_snps_summary_stats$beta, use = 'pairwise' )
        removed_traits  =  FALSE
        while (any( eigen(correlation_matrix)$values < 0 )) {
            removed_traits  =  TRUE
            if (all(eigen(correlation_matrix * (1-1e-5) + diag(ncol(correlation_matrix)) * 1e-5)$values >= 0)) {
                print('Regularization works')
            }
            max_cor  =  correlation_matrix[ correlation_matrix != 1 ] %>%
                abs %>%
                max
            max_cor  =  correlation_matrix[ abs(correlation_matrix) == max_cor ][ 1 ]
            to_remove  =  which( abs(correlation_matrix) == abs( max_cor ), arr.ind = TRUE ) %>%
                rownames
            to_remove  =  correlation_matrix[ , to_remove ] %>%
                abs %>%
                colSums %>%
                which.max %>%
                names
            print( sprintf( 'Removing %s, correlation is %s', to_remove, max_cor ) )
            to_keep  =  setdiff( colnames(correlation_matrix), to_remove )
            correlation_matrix  =  correlation_matrix[ to_keep, ][ , to_keep ]
        }
        if (removed_traits) {
            all_snps_summary_stats  =  all_snps_summary_stats %>%
                lapply( function(snp_summary_stats){
                    snp_summary_stats[ , colnames( correlation_matrix ) ]
                } )
            all_phenotypes  =  all_phenotypes %>%
                filter( phenotype %in% colnames( correlation_matrix ) )
        }
    }
    
    all_snps_summary_stats$beta[ is.na( all_snps_summary_stats$beta ) ]  =  0
    all_snps_summary_stats$se[   is.na( all_snps_summary_stats$se   ) ]  =  0
    all_snps_summary_stats$pval[ is.na( all_snps_summary_stats$pval ) ]  =  1
    
    all_snps_summary_stats   =  all_snps_summary_stats  %>%
        lapply( as.matrix )
    
    pca_summary_stats   =  prcomp( all_snps_summary_stats$beta,
                                   center = FALSE,
                                   scale = scaling )
    if ( line_up_positive != FALSE ) {
        pca_summary_stats[ c( 'rotation', 'x' ) ]  =  pca_summary_stats[ c( 'rotation', 'x' ) ] %>%
            lapply( sweep,
                    MARGIN = 2,
                    sign( pca_summary_stats$rotation[ line_up_positive, ] ),
                    '*' )
    }
    
    pca_summary_stats$rotation  =  ( t( pca_summary_stats$rotation )
                                     %*% correlation_matrix
                                     %*% pca_summary_stats$rotation ) %>%
        diag %>%
        sapply( function( x ) max( x, 0 ) ) %>%
        sqrt %>%
        sweep( pca_summary_stats$rotation, MARGIN = 2, ., '/' )
    
    pca_summary_stats$x  =  all_snps_summary_stats$beta %*% pca_summary_stats$rotation
    
    pca_beta  =  pca_summary_stats$x[ , 1:min( npcs, ncol( pca_summary_stats$x ) ) ]
    
    npcs = min( npcs, ncol( pca_summary_stats$rotation ) )
    
    pca_se  =  matrix( nrow     = dim( pca_beta )[ 1 ],
                       ncol     = dim( pca_beta )[ 2 ],
                       dimnames = dimnames( pca_beta ) )
    
    for (pc in 1:npcs) {
        partial_var  =  sweep( all_snps_summary_stats$se,
                               MARGIN = 2,
                               pca_summary_stats$rotation[ , pc ],
                               '*' )
        
        pca_se[ , pc ]  =  ( ( partial_var %*% correlation_matrix ) * partial_var ) %>%
            rowSums %>%
            sqrt
    }
    
    pca_pval  =  2 * pnorm( -abs( pca_beta / pca_se ) )
    colnames( pca_pval )  =  colnames( pca_beta )
    pca_summary_stats$correlation_matrix  =  correlation_matrix
    
    list( beta = pca_beta, se = pca_se, pval = pca_pval, pca = pca_summary_stats )
}

calculate_pc_from_loadings  =  function( all_snps_summary_stats,
                                         pca_summary_stats = NULL,
                                         npcs = 10,
                                         restrict_npcs = TRUE,
                                         return_variants = FALSE ){
    
    if (return_variants) {
        variants = all_snps_summary_stats$beta %>%
            select( one_of( c( 'variant', 'rsid' ) ) )
    }
    
    all_snps_summary_stats$beta[ is.na( all_snps_summary_stats$beta ) ]  =  0
    all_snps_summary_stats$se[   is.na( all_snps_summary_stats$se   ) ]  =  0
    all_snps_summary_stats$pval[ is.na( all_snps_summary_stats$pval ) ]  =  1
    
    all_snps_summary_stats   =  all_snps_summary_stats  %>%
        # lapply( as.matrix )
        lapply( function(x) select_if( x, is.numeric ) %>% as.matrix )
    
    ydim  =  ncol( pca_summary_stats$rotation )
    
    pca_summary_stats$rotation  =  pca_summary_stats$rotation[ ,
                                                               1:ydim,
                                                               drop = FALSE ]
    correlation_matrix = pca_summary_stats$correlation_matrix
    rownames( correlation_matrix ) = colnames( correlation_matrix ) = colnames(correlation_matrix) %>%
        str_remove( '_irnt' )
    if (restrict_npcs) {
        npcs = min( npcs, ncol( pca_summary_stats$rotation ) )
    } else {
        npcs = ncol( pca_summary_stats$rotation )
    }
    
    pca_summary_stats$rotation  =  pca_summary_stats$rotation[ , 1:npcs, drop = FALSE ] %>%
        apply( 2, function( x ) x / sqrt(sum(x^2)) )
    pca_summary_stats$rotation  =  ( t( pca_summary_stats$rotation )
                                     %*% correlation_matrix
                                     %*% pca_summary_stats$rotation ) %>%
        diag %>%
        sapply( function( x ) max( x, 0 ) ) %>%
        sqrt %>%
        sweep( pca_summary_stats$rotation, MARGIN = 2, ., '/' )
    
    if (!is.null( pca_summary_stats$sdev )) {
        cnames  =  ( pca_summary_stats$sdev[ 1:npcs ]^2 * 100 / sum( pca_summary_stats$sdev^2 ) ) %>%
            signif( digits = 3 ) %>%
            paste0( 'PC ', 1:npcs, ' (', ., '%)' )
        colnames( pca_summary_stats$rotation )  =  cnames
    }
    
    pca_beta  =  all_snps_summary_stats$beta %*% pca_summary_stats$rotation
    
    pca_se  =  matrix( nrow     = dim( pca_beta )[ 1 ],
                       ncol     = dim( pca_beta )[ 2 ],
                       dimnames = dimnames( pca_beta ) )
    
    for (pc in 1:npcs) {
        partial_var  =  sweep( all_snps_summary_stats$se,
                               MARGIN = 2,
                               pca_summary_stats$rotation[ , pc ],
                               '*' )
        pca_se[ , pc ]  =  ( ( partial_var %*% correlation_matrix ) * partial_var ) %>%
            rowSums %>%
            sqrt
    }
    
    pca_pval  =  2 * pnorm( -abs( pca_beta / pca_se ) )
    colnames( pca_pval )  =  colnames( pca_beta )
    
    if (return_variants) {
        pca_beta  =  cbind( variants, pca_beta )
        pca_se    =  cbind( variants, pca_se   )
        pca_pval  =  cbind( variants, pca_pval )
    }
    
    list( beta = pca_beta, se = pca_se, pval = pca_pval, pca = pca_summary_stats )
}










