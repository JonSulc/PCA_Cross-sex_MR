require( tidyverse,   quietly = TRUE, warn.conflicts = FALSE )
require( TwoSampleMR, quietly = TRUE )
require( RNOmni )

source( 'pca.R' )

# Resolve conflicts
filter  =  dplyr::filter

# Default values
THRESHOLD  =  5e-8
MIN_IVS_BEFORE  =  10
MIN_IVS_AFTER   =  5
HET_THRESHOLD   =  1e-4
REVERSE_T_THRESHOLD  =  qnorm( 5e-2 ) # p-value of 5e-2 one-sided

mr_trait_vs_trait  =  function( all_snps_summary_stats,
                                all_phenotypes,
                                exposure_category,
                                exposure_snps,
                                outcome_category,
                                ... ){
    exposure  =  .select_all_category_stats( all_snps_summary_stats,
                                             all_phenotypes,
                                             exposure_category,
                                             exposure_snps )
    outcome   =  .select_all_category_stats( all_snps_summary_stats,
                                             all_phenotypes,
                                             outcome_category,
                                             exposure_snps )
    
    .perform_all_mr( exposure,
                     outcome,
                     exposure_snps,
                     ... )
}

mr_pc_vs_trait  =  function( all_snps_summary_stats,
                             all_phenotypes,
                             exposure_pca,
                             exposure_snps,
                             outcome_category,
                             ... ){

    outcome  =  .select_all_category_stats( all_snps_summary_stats,
                                            all_phenotypes,
                                            outcome_category,
                                            exposure_snps )

    .perform_all_mr( exposure_pca[ c( 'beta', 'pval', 'se' ) ],
                     outcome,
                     exposure_snps,
                     ... )
}

mr_trait_vs_pc  =  function( all_snps_summary_stats,
                             all_phenotypes,
                             exposure_category,
                             exposure_snps,
                             outcome_pca,
                             outcome_category,
                             npcs = 10,
                             restrict_npcs = TRUE,
                             all_phenotypes_outcome = all_phenotypes,
                             ... ){
    all_phenotypes  =  all_phenotypes %>%
        inner_join( all_phenotypes_outcome )
    exposure  =  .select_all_category_stats( all_snps_summary_stats,
                                             all_phenotypes,
                                             exposure_category,
                                             exposure_snps ) %>%
        lapply( as.matrix )
    
    outcome_pca  =  .select_all_category_stats( all_snps_summary_stats,
                                                all_phenotypes,
                                                outcome_category,
                                                exposure_snps ) %>%
        calculate_pc_from_loadings( outcome_pca$pca,
                                    npcs = npcs,
                                    restrict_npcs = restrict_npcs )
    
    # outcome_pca  =  get_category_pca_stats( all_snps_summary_stats,
    #                                         all_phenotypes,
    #                                         outcome_category,
    #                                         exposure_snps,
    #                                         pca_summary_stats = outcome_pca$pca,
    #                                         npcs = npcs )
    
    .perform_all_mr( exposure,
                     outcome_pca,
                     exposure_snps,
                     ... )
}


mr_pc_vs_pc  =  function( all_snps_summary_stats,
                          all_phenotypes,
                          exposure_pca,
                          exposure_snps,
                          outcome_pca,
                          outcome_category,
                          npcs = 10,
                          restrict_npcs = TRUE,
                          ... ){
    
    outcome_pca  =  .select_all_category_stats( all_snps_summary_stats,
                                                all_phenotypes,
                                                outcome_category,
                                                exposure_snps ) %>%
        calculate_pc_from_loadings( outcome_pca$pca,
                                    npcs = npcs,
                                    restrict_npcs = restrict_npcs )
    
    # outcome_pca  =  get_category_pca_stats( all_snps_summary_stats,
    #                                         all_phenotypes,
    #                                         outcome_category,
    #                                         exposure_snps,
    #                                         pca_summary_stats = outcome_pca$pca )
    
    .perform_all_mr( exposure_pca[ c( 'beta', 'pval', 'se' ) ],
                     outcome_pca[  c( 'beta', 'pval', 'se' ) ],
                     exposure_snps,
                     ... )
}


.perform_all_mr  =  function( exposure,
                              outcome,
                              exposure_snps,
                              threshold = THRESHOLD,
                              ... ){
    
    exposure  =  exposure %>%
        lapply( as.matrix )
    outcome   =  outcome  %>%
        lapply( as.matrix )
    
    mr_beta  =  data.frame( matrix( NA, nrow = ncol( exposure$beta ), ncol = ncol( outcome$beta ) ) )
    colnames( mr_beta )  =  colnames( outcome$beta  )
    rownames( mr_beta )  =  colnames( exposure$beta )
    
    mr_results  =  list( b                  = mr_beta,
                         se                 = mr_beta,
                         pval               = mr_beta,
                         egger_pval         = mr_beta,
                         heterogeneity_pval = mr_beta,
                         snps_start         = mr_beta,
                         snps_end           = mr_beta )
    mr_results$snps_start[]  =  colSums( exposure$pval < threshold )
    
    for (trait_name in colnames( exposure$beta )) {
        trait  =  list( beta = exposure$beta[ , trait_name, drop = FALSE ],
                        se   = exposure$se[   , trait_name, drop = FALSE ],
                        pval = exposure$pval[ , trait_name, drop = FALSE ] )
        for (outcome_category in colnames( outcome$beta )) {
            mr_result  =  single_mr( trait,
                                     list( beta   = outcome$beta[ , outcome_category, drop = FALSE ],
                                           se     = outcome$se[   , outcome_category, drop = FALSE ],
                                           pval   = outcome$pval[ , outcome_category, drop = FALSE ] ),
                                     exposure_snps$rsid,
                                     exposure_snps$alt,
                                     exposure_snps$ref,
                                     exposure_snps$AF,
                                     threshold = threshold,
                                     ... )
            # Save results and begin anew
            for (stat in c( 'b', 'se', 'pval' )) {
                mr_results[[ stat ]][ trait_name, outcome_category ]  =  mr_result$ivw[ 1, stat ]
            }
            
            if (is.null( mr_result$egger )) {
                mr_results$egger_pval[ trait_name, outcome_category ]  =  NA
            } else {
                mr_results$egger_pval[ trait_name, outcome_category ]  =  mr_result$egger
            }
            
            if (is.null( mr_result$heterogeneity )) {
                mr_results$heterogeneity_pval[ trait_name, outcome_category ]  =  NA
            } else {
                mr_results$heterogeneity_pval[ trait_name, outcome_category ]  =  mr_result$heterogeneity
            }
            
            mr_results$snps_end[ trait_name, outcome_category ]  =  mr_result$snps_end
        }
    }
    
    mr_results
}


single_mr  =  function( single_exposure,
                        single_outcome,
                        snp_names,
                        ea,
                        oa,
                        eaf,
                        threshold = THRESHOLD,
                        min_ivs_before = MIN_IVS_BEFORE,
                        min_ivs_after = MIN_IVS_AFTER,
                        het_threshold = HET_THRESHOLD,
                        reverse_t_threshold = REVERSE_T_THRESHOLD,
                        ... ){
    exposure_data  =  data.frame( snp_names,
                                  single_exposure_pval = unname( single_exposure$pval ),
                                  single_exposure_beta = unname( single_exposure$beta ),
                                  single_exposure_se   = unname( single_exposure$se   ),
                                  ea,
                                  oa,
                                  eaf,
                                  exposure_phenotype = 'exposure',
                                  stringsAsFactors = FALSE ) %>%
        filter( single_exposure_pval < threshold )
    if (nrow( exposure_data ) < min_ivs_before) {
        return( list( mr            = array( NA,
                                             dim = c( 1, 3 ),
                                             dimnames = list( NULL, c( 'b', 'se', 'pval' ) ) ),
                      egger         = NA,
                      heterogeneity = NA,
                      snps_end      = 0 ) )
    }
    exposure_data  =  exposure_data %>%
        format_data( type              = 'exposure',
                     phenotype_col     = 'exposure_phenotype',
                     snp_col           = 'snp_names',
                     beta_col          = 'single_exposure_beta',
                     se_col            = 'single_exposure_se',
                     eaf_col           = 'eaf',
                     effect_allele_col = 'ea',
                     other_allele_col  = 'oa',
                     pval_col          = 'single_exposure_pval' )
    
    mr_result  =  data.frame( snp_names,
                              outcome_pval = unname( single_outcome$pval ),
                              outcome_beta = unname( single_outcome$beta ),
                              outcome_se   = unname( single_outcome$se   ),
                              ea,
                              oa,
                              eaf,
                              outcome_phenotype = 'outcome',
                              stringsAsFactors = FALSE ) %>%
        format_data( type              = 'outcome',
                     phenotype_col     = 'outcome_phenotype',
                     snp_col           = 'snp_names',
                     beta_col          = 'outcome_beta',
                     se_col            = 'outcome_se',
                     eaf_col           = 'eaf',
                     effect_allele_col = 'ea',
                     other_allele_col  = 'oa',
                     pval_col          = 'outcome_pval' ) %>%
        harmonise_data( exposure_dat = exposure_data,
                        outcome_dat  = .,
                        action = 1 ) %>%
        filter( ( abs( beta.exposure ) - abs( beta.outcome ) ) / sqrt( se.exposure^2 + se.outcome^2 ) > reverse_t_threshold ) %>%
        simple_mr( min_ivs_before,
                   min_ivs_after,
                   het_threshold,
                   ... )
}

mr_each_trait_vs_trait  =  function( all_exposure_stats,
                                     all_phenotypes,
                                     exposure_category,
                                     exposure_snps,
                                     outcome_category,
                                     all_outcome_stats = all_exposure_stats,
                                     ... ){
    exposure  =  .select_each_category_stats( all_exposure_stats,
                                              all_phenotypes,
                                              exposure_category,
                                              exposure_snps )
    
    outcome   =  exposure_snps %>%
        lapply( function( trait ) trait$variant ) %>%
        unlist %>%
        unique %>%
        .select_all_category_stats( all_outcome_stats,
                                    all_phenotypes,
                                    outcome_category,
                                    .,
                                    return_variants = TRUE )
    
    .perform_each_mr( exposure,
                      outcome,
                      exposure_snps,
                      ... )
}

mr_each_pc_vs_trait  =  function( all_outcome_stats,
                                  all_phenotypes,
                                  exposure_pcs,
                                  exposure_snps,
                                  outcome_category,
                                  ... ){
    exposure  =  exposure_pcs[ names(exposure_pcs) != 'pca' ]
    if (length(exposure) < 1) {
        return()
    }
    
    variants  =  exposure_snps %>%
        lapply( function( pc ) pc$variant ) %>%
        unlist %>%
        unique
    outcome  =  .select_all_category_stats( all_outcome_stats,
                                            all_phenotypes,
                                            outcome_category,
                                            variants,
                                            return_variants = TRUE )
    # outcome %>%
    #     lapply( function( stat ){
    #         stat %>%
    #             mutate( variant = variants ) %>%
    #             select( variant, everything() )
    #     } ) %>%
    .perform_each_mr( exposure,
                      outcome,
                      exposure_snps,
                      ... )
}

mr_each_trait_vs_pc  =  function( all_exposure_stats,
                                  all_phenotypes,
                                  exposure_category,
                                  exposure_snps,
                                  outcome_pca,
                                  outcome_category,
                                  npcs = 10,
                                  restrict_npcs = TRUE,
                                  all_phenotypes_outcome = all_phenotypes,
                                  all_outcome_stats = all_exposure_stats,
                                  ... ){
    if (ncol( outcome_pca$pca$rotation ) == 0) {
        return()
    }
    
    all_phenotypes  =  all_phenotypes %>%
        inner_join( all_phenotypes_outcome )
    exposure  =  .select_each_category_stats( all_exposure_stats,
                                              all_phenotypes,
                                              exposure_category,
                                              exposure_snps )
    
    variants  =  exposure_snps %>%
        lapply( function( pc ) pc$variant ) %>%
        unlist %>%
        unique
    outcome_pca  =  variants %>%
        .select_all_category_stats( all_outcome_stats,
                                    all_phenotypes,
                                    outcome_category,
                                    .,
                                    return_variants = TRUE ) %>%
        calculate_pc_from_loadings( outcome_pca$pca,
                                    npcs = npcs,
                                    restrict_npcs = restrict_npcs )
    outcome  =  outcome_pca[  c( 'beta', 'pval', 'se' ) ] %>%
        lapply( function( x ) {
            tibble( variant = variants ) %>%
                cbind( x )
        } )
    
    .perform_each_mr( exposure,
                      outcome,
                      exposure_snps,
                      ... )
}


mr_each_pc_vs_pc  =  function( all_outcome_stats,
                               all_phenotypes,
                               exposure_pcs,
                               exposure_snps,
                               outcome_pca,
                               outcome_category,
                               npcs = 10,
                               restrict_npcs = TRUE,
                               ... ){
    if (ncol( outcome_pca$pca$rotation ) == 0) {
        return()
    }
    
    exposure  =  exposure_pcs[ names(exposure_pcs) != 'pca' ]
    if (length(exposure) < 1) {
        return()
    }
    
    variants  =  exposure_snps %>%
        lapply( function( pc ) pc$variant ) %>%
        unlist %>%
        unique
    outcome_pca  =  variants %>%
        .select_all_category_stats( all_outcome_stats,
                                    all_phenotypes,
                                    outcome_category,
                                    .,
                                    return_variants = TRUE ) %>%
        calculate_pc_from_loadings( outcome_pca$pca,
                                    npcs = npcs,
                                    restrict_npcs = restrict_npcs )
    outcome  =  outcome_pca[  c( 'beta', 'pval', 'se' ) ] %>%
        lapply( function( x ) {
            tibble( variant = variants ) %>%
                cbind( x )
        } )
    
    .perform_each_mr( exposure,
                      outcome,
                      exposure_snps,
                      ... )
}


.perform_each_mr  =  function( exposure,
                               outcome,
                               exposure_snps,
                               threshold = THRESHOLD,
                               ... ){
    
    mr_beta  =  data.frame( matrix( NA,
                                    nrow = length( exposure ),
                                    ncol = ncol( outcome$beta )-1 ) )
    colnames( mr_beta )  =  colnames( outcome$beta  )[ -1 ]
    rownames( mr_beta )  =     names( exposure )
    
    mr_results  =  list( b                  = mr_beta,
                         se                 = mr_beta,
                         pval               = mr_beta,
                         egger_pval         = mr_beta,
                         heterogeneity_pval = mr_beta,
                         snps_start         = mr_beta,
                         snps_end           = mr_beta )
    mr_results$snps_start[]  =  sapply( exposure, function(x) sum( x$pval < threshold ) )
    
    for (exposure_name in names( exposure )) {
        exposure_stats  =  exposure[[ exposure_name ]]
        all_outcome_stats  =  c( 'beta', 'se', 'pval' ) %>%
            setNames( nm = . ) %>%
            lapply( function( stat ){
                exposure_snps[[ exposure_name ]] %>%
                    select( variant ) %>%
                    left_join( outcome[[ stat ]] )
                # outcome[[ stat ]] %>%
                #     left_join( select( exposure_snps[[ exposure_name ]], variant ), . )
            } )
        for (outcome_name in colnames( all_outcome_stats$beta )[ -1 ]) {
            outcome_stats  =  c( 'beta', 'se', 'pval' ) %>%
                setNames( nm = . ) %>%
                lapply( function( stat ){
                    all_outcome_stats[[ stat ]] %>%
                        select( !!stat := !!outcome_name )
                } ) %>%
                bind_cols#( all_outcome_stats$beta %>% select( variant, rsid ) )
            mr_result  =  single_mr( single_exposure = exposure_stats,
                                     single_outcome  = outcome_stats,
                                     snp_names       = exposure_snps[[ exposure_name ]]$rsid,
                                     ea              = exposure_snps[[ exposure_name ]]$alt,
                                     oa              = exposure_snps[[ exposure_name ]]$ref,
                                     eaf             = exposure_snps[[ exposure_name ]]$AF,
                                     threshold = threshold,
                                     ... )
            # Save results and begin anew
            for (stat in c( 'b', 'se', 'pval' )) {
                mr_results[[ stat ]][ exposure_name, outcome_name ]  =  mr_result$mr[ 1, stat ]
            }
            
            if (is.null( mr_result$egger )) {
                mr_results$egger_pval[ exposure_name, outcome_name ]  =  NA
            } else {
                mr_results$egger_pval[ exposure_name, outcome_name ]  =  mr_result$egger
            }
            
            if (is.null( mr_result$heterogeneity )) {
                mr_results$heterogeneity_pval[ exposure_name, outcome_name ]  =  NA
            } else {
                mr_results$heterogeneity_pval[ exposure_name, outcome_name ]  =  mr_result$heterogeneity
            }
            
            mr_results$snps_end[ exposure_name, outcome_name ]  =  mr_result$snps_end
        }
    }
    
    mr_results
}

simple_mr  =  function( harmonised_data,
                        min_ivs_before,
                        min_ivs_after,
                        het_threshold = HET_THRESHOLD,
                        mr_method = 'mr_ivw' ){
    while (TRUE) {
        if (nrow( harmonised_data ) < min_ivs_before) {
            return( list( mr            = array( NA,
                                                 dim = c( 1, 3 ),
                                                 dimnames = list( NULL, c( 'b', 'se', 'pval' ) ) ),
                          egger         = NA,
                          heterogeneity = NA,
                          snps_end      = 0 ) )
        }
        mr_results  =  list( mr = suppressMessages( mr( harmonised_data, method_list = c( mr_method ) ) ),
                             egger = suppressMessages( mr_pleiotropy_test( harmonised_data )[ 1, 'pval' ] ),
                             snps_end      = nrow( harmonised_data ) )
        if (mr_method %in% mr_method_list()$obj[mr_method_list()$heterogeneity_test]) {
            mr_results$heterogeneity = suppressMessages( mr_heterogeneity( harmonised_data,
                                                                           method_list = c( mr_method ) )[ 1, 'Q_pval' ] )
        }
        
        d  =  harmonised_data$beta.outcome - mr_results$mr$b * harmonised_data$beta.exposure
        
        var_d  =  harmonised_data$se.outcome^2 +
            harmonised_data$beta.exposure^2 * mr_results$mr$se^2 +
            harmonised_data$se.exposure^2   * mr_results$mr$b^2  +
            harmonised_data$se.exposure^2   * mr_results$mr$se^2
        
        z  =  d^2 / var_d
        if (min( 1 - pchisq( z, 1 ) ) > het_threshold) {
            return( mr_results )
        }
        harmonised_data  =  harmonised_data[ z != max( z ), ]
        min_ivs_before   =  min_ivs_after
    }
}

get_category_fa_stats  =  function( all_snps_summary_stats,
                                    all_phenotypes,
                                    category_name,
                                    category_snps,
                                    phenotype_file = 'pheno/ukb21067.csv',
                                    line_up_positive = FALSE,
                                    fa_summary_stats = NULL,
                                    nfactors = 10 ){
    .select_all_category_stats( all_snps_summary_stats,
                                all_phenotypes,
                                category_name,
                                category_snps ) %>%
        get_fa_summary_stats( phenotype_file = phenotype_file,
                              line_up_positive = line_up_positive,
                              fa_summary_stats = fa_summary_stats,
                              nfactors = nfactors ) %>%
        return
}

get_fa_summary_stats  =  function( all_snps_summary_stats,
                                   phenotype_file = NULL,
                                   line_up_positive = FALSE,
                                   fa_summary_stats = NULL,
                                   nfactors = 10,
                                   factor_function = factanal,
                                   ... ){
    
    all_snps_summary_stats$beta[ is.na( all_snps_summary_stats$beta ) ]  =  0
    all_snps_summary_stats$se[   is.na( all_snps_summary_stats$se   ) ]  =  0
    all_snps_summary_stats$pval[ is.na( all_snps_summary_stats$pval ) ]  =  1
    
    all_snps_summary_stats   =  all_snps_summary_stats  %>%
        lapply( as.matrix )
    
    # Pick up modifications here
    if ( is.null( fa_summary_stats ) ){
        fa_summary_stats   =  factor_function( all_snps_summary_stats$beta,
                                               nfactors,
                                               ... )
        # return( fa_summary_stats )
        if ( line_up_positive != FALSE ) {
            fa_summary_stats[ c( 'rotmat', 'loadings' ) ]  =  fa_summary_stats[ c( 'rotmat', 'loadings' ) ] %>%
                lapply( sweep,
                        MARGIN = 2,
                        sign( fa_summary_stats$rotmat[ line_up_positive, ] ),
                        '*' )
        }
        
        fa_beta             =  fa_summary_stats$loadings[ , 1:nfactors ]
        correlation_matrix  =  cor( all_snps_summary_stats$beta )
    } else {
        fa_beta  =  all_snps_summary_stats$beta %*% fa_summary_stats$rotmat[ , 1:nfactors ]
        correlation_matrix  =  cor( fa_summary_stats$loadings %*% t(fa_summary_stats$rotmat^(-1)) )
    }
    
    
    # colnames( fa_beta )  =  ( fa_summary_stats$sdev[ 1:nfactors ]^2 * 100 / sum( fa_summary_stats$sdev^2 ) ) %>%
    #     signif( digits = 3 ) %>%
    #     paste0( 'Fa ', 1:nfactors, ' (', ., '%)' )
    
    fa_se  =  matrix( nrow     = dim( fa_beta )[ 1 ],
                      ncol     = dim( fa_beta )[ 2 ],
                      dimnames = dimnames( fa_beta ) )
    
    for (factor in 1:nfactors) {
        partial_var  =  sweep( all_snps_summary_stats$se,
                               MARGIN = 2,
                               fa_summary_stats$rotmat[ , factor ],
                               '*' )
        for (snp in 1:dim( fa_beta )[ 1 ]) {
            fa_se[ snp, factor ]  =  ( partial_var[ snp, , drop = FALSE ]
                                    %*% correlation_matrix
                                    %*% t( partial_var[ snp, , drop = FALSE ] ) ) %>%
                sqrt
        }
    }
    
    fa_pval      =  2 * pnorm( -abs( fa_beta / fa_se ) )
    colnames( fa_pval )  =  colnames( fa_beta )
    
    list( beta = fa_beta, se = fa_se, pval = fa_pval, fa = fa_summary_stats )
}

.convert_names  =  function( phenotypes,
                             all_phenotypes,
                             minlength = 20 ){
    phenotypes %>%
        tibble( phenotype = . ) %>%
        left_join( all_phenotypes ) %>%
        '$'( description ) %>%
        coalesce( phenotypes ) %>%
        abbreviate( minlength = minlength )
}









