require( tidyverse,   quietly = TRUE, warn.conflicts = FALSE )

.select_all_category_stats  =  function( all_snps_summary_stats,
                                         all_phenotypes,
                                         selected_category,
                                         variant_list,
                                         ... ){
    lapply( all_snps_summary_stats,
            .select_category_stats,
            all_phenotypes = all_phenotypes,
            selected_category = selected_category,
            variant_list = variant_list,
            ... )
}

.select_category_stats  =  function( summary_stats,
                                     all_phenotypes,
                                     selected_category,
                                     variant_list,
                                     return_variants = FALSE ){
    
    if (!is.null( dim( variant_list ) ))
        variant_list  =  variant_list$variant
    
    phenotypes  =  all_phenotypes %>%
        dplyr::filter( category == selected_category ) %>%
        '$'( phenotype )
    
    if (return_variants)
        phenotypes = c( 'variant', phenotypes )
    
    summary_stats %>%
        dplyr::left_join( tibble( variant = variant_list ), . ) %>%
        select( one_of( phenotypes ) )
}

.select_each_category_stats  =  function( all_snps_summary_stats,
                                          all_phenotypes,
                                          selected_category,
                                          variant_lists,
                                          ... ){
    phenotypes  =  all_phenotypes %>%
        dplyr::filter( category == selected_category ) %>%
        dplyr::filter( phenotype %in% names( variant_lists ) ) %>%
        '$'( 'phenotype' ) %>%
        setNames( nm = . )

    lapply( phenotypes,
            .select_single_trait_stats,
            all_snps_summary_stats,
            variant_lists )
}

.select_each_pc_stats  =  function( all_snps_summary_stats,
                                    all_phenotypes,
                                    selected_category,
                                    variant_lists,
                                    pca_summary_stats,
                                    ... ){
    
    phenotypes  =  all_phenotypes %>%
        dplyr::filter( category == selected_category ) %>%
        '$'( phenotype )
    
    names( variant_lists ) %>%
        setNames( nm = . ) %>%
        lapply( function( pc_name ) {
            pca = pca_summary_stats
            pca$rotation = pca$rotation[ , pc_name, drop = FALSE ]
            to_keep  =  all_snps_summary_stats$beta$variant %in% variant_lists[[ pc_name ]]$variant
            all_snps_summary_stats %>%
                lapply( function( summary_stats ){
                    summary_stats %>%
                        dplyr::filter( to_keep ) %>%
                        select( one_of( phenotypes ) )
                } ) %>%
                calculate_pc_from_loadings( pca, npcs = 1 ) %>%
                '['( c( 'beta', 'se', 'pval' ) ) %>%
                lapply( unname ) %>%
                lapply( c ) %>%
                bind_cols
        })
    
    # variant_lists %>%
    #     lapply( function( pc_variants ) {
    #         all_snps_summary_stats %>%
    #             lapply( function( summary_stats ){
    #                 summary_stats %>%
    #                     filter( variant %in% pc_variants$variant ) %>%
    #                     select( one_of( phenotypes ) )
    #             } ) %>%
    #             calculate_pc_from_loadings( pca_summary_stats )
    #     } )
}

.select_single_trait_stats  =  function( trait,
                                         all_snps_summary_stats,
                                         variants,
                                         ... ){
    if (!(trait %in% colnames( all_snps_summary_stats$beta ))) {
        return()
    }
    
    if (is_list(variants))
        variants  =  variants[[ trait ]]
    
    if (!is.null( dim( variants ) ))
        variants  =  variants$variant
    
    c( 'beta', 'se', 'pval' ) %>%
        lapply( function( stat ){
            all_snps_summary_stats[[ stat ]] %>%
                # select( variant, rsid, one_of( !!trait ) ) %>%
                # rename( list( !!stat = !!trait ) ) %>%
                select( variant, rsid, !!stat := !!trait ) %>%
                left_join( tibble( variant = variants ), . )
        } ) %>%
    reduce( full_join )
}

.convert_names  =  function( phenotypes,
                             all_phenotypes,
                             min_length = 50 ){
    phenotypes %>%
        tibble( phenotype = . ) %>%
        left_join( all_phenotypes ) %>%
        '$'( description ) %>%
        coalesce( phenotypes ) %>%
        abbreviate( minlength = min_length )
}









