require( tidyverse, quietly = TRUE )
require( data.table )

source( 'get_snp_stats.R' )
source( 'pca.R'  )
source( 'download_neale_files.R' )

# Parameters
SNP_THRESHOLD  =  5e-8  # P-value threshold for SNP selection
RAW            =  FALSE # Use raw phenotypes or inverse-rank normal transformed
NPCS           =  4     # Number of PCs to calculate
HETEROGENEITY       = 1e-3
MIN_IVS_BEFORE      = 10
MIN_IVS_AFTER       = 5
REVERSE_T_THRESHOLD = qnorm( .05 )

sexes = c( 'both_sexes', 'female', 'male' )
exposure_phenotypes  =  read_tsv( 'data/exposure_phenotypes.tsv' )
outcome_phenotypes   =  read_tsv( 'data/outcome_phenotypes.tsv'  )

# Downloading the summary statistics ####
# Will ignore any existing files in the appropriate folders
timeout = getOption( 'timeout' )
options( timeout = 1e3 )
download_neale_files( exposure_phenotypes$phenotype,
                      chunk_size = 3,
                      category = 'body' )
download_neale_files( outcome_phenotypes$phenotype,
                      chunk_size = 3,
                      category = 'disease' )
if (!file.exists('data/neale_files/variants.tsv.gz')) {
    download.file( 'https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/annotations/variants.tsv.bgz',
                   'data/neale_files/variants.tsv.gz' )
}

for (sex in sexes) {
    if (!file.exists( 'data/neale_files/phenotypes_%s.tsv.gz' %>%
                      sprintf( sex ) )) {
        download.file( 'https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/annotations/phenotypes.%s.tsv.bgz' %>%
                           sprintf( sex ),
                       'data/neale_files/phenotypes_%s.tsv.gz' %>%
                           sprintf( sex ) )
    }
}

options( timeout = timeout )

# Formatting the files for faster access ####
all_phenotypes  =  bind_rows( exposure_phenotypes %>%
                                  filter( phenotype != 'whr' ),
                              outcome_phenotypes )
snp_positions  =  fread( 'zcat data/neale_files/variants.tsv.gz', select = c( 'variant' ) ) %>%
    mutate( keep = variant %chin% COMMON_SNPS$variant )

for (sex in sexes) {
    dir.create( sprintf( 'data/%s', sex ),
                recursive = TRUE,
                showWarnings = FALSE )
    walk( 1:(nrow(all_phenotypes)),
          function( i ){
              neale_filename  =  'data/neale_files/%s/%s' %>%
                  sprintf( sex,
                           all_phenotypes$category[i] ) %>%
                  list.files( pattern = sprintf( '%s(_irnt)?.gwas.imputed_v3.%s.tsv.gz',
                                                 all_phenotypes$phenotype[i],
                                                 sex ),
                              full.names = TRUE ) %>%
                  '['(1)
              
              neale_output_path  =  sprintf( 'data/processed/%s/%s',
                                             sex,
                                             all_phenotypes$category[i] )
              dir.create( neale_output_path, showWarnings = FALSE, recursive = TRUE )
              neale_output_filename  =  sprintf( '%s/%s.rds',
                                                 neale_output_path,
                                                 all_phenotypes$phenotype[i] )
              
              
              if ( !file.exists( neale_filename )) {
                  stop( sprintf( 'File missing: %s',
                                 neale_filename ) )
              }
              stats = fread( cmd = sprintf( 'zcat %s', neale_filename ),
                             select = c( 'variant', 'low_confidence_variant', 'beta', 'se', 'pval' ) )
              stats  =  stats %>%
                  select( -pval ) %>%
                  left_join( snp_positions[ snp_positions$keep, 'variant', drop = FALSE ], . ) %>%
                  transmute_if( is.numeric,
                                function( stat ){
                                    replace( stat,
                                             str_detect( .$low_confidence_variant, regex( 'true', ignore_case = TRUE ) ),
                                             NA )
                                } )
              
              # Divide betas of binary phenotypes by the phenotype SD to obtain
              # beta estimates on a standardized scale (variance = 1)
              if (all_phenotypes$category[i] == 'disease') {
                  neale_phenotypes  =  fread( cmd = sprintf( 'zcat data/neale_files/phenotypes_%s.tsv.gz',
                                                             sex ) )
                  disease  =  which( all_phenotypes$phenotype[i] == neale_phenotypes$phenotype ) %>%
                      '['( neale_phenotypes, . )
                  stats  =  stats/ sqrt( disease$n_cases/disease$n_non_missing * disease$n_controls/disease$n_non_missing )
              }
              if (any( !is.na( stats ) )) {
                  saveRDS( stats %>% as.data.frame,
                           neale_output_filename )
              }
          })
}

# # Rank and clump SNPs, compute PCs ####
# 
# # The clumping requires a refernce LD structure.
# # UK10K (EGAD00001000740, EGAD00001000741) was used for the analyses presented
# # in the manuscript, however the 1000 Genomes data could be used
# # (e.g. via the TwoSampleMR package).
# # The code is here for reference but the subsequent files of clumped SNPs are
# # included with the other files and this section can be skipped.
# 
# all_snps_summary_stats  =  get_variables_from_all( 'body',
#                                                    all_phenotypes = exposure_phenotypes,
#                                                    variables = c( 'beta', 'se', 'pval' ),
#                                                    sex = 'both_sexes',
#                                                    threshold = FALSE )
# 
# traits_reclumped_snps  =  all_snps_summary_stats %>%
#     clump_each( threshold = SNP_THRESHOLD ) %>%
#     lapply( function( trait ) {
#         trait %>%
#             select( variant ) %>%
#             left_join( COMMON_SNPS )
#     } )
# write_rds( traits_reclumped_snps,
#            'data/clumped_snps/both_sexes_body_traits_clumped_snps.rds' )
# 
# all_snps_summary_stats  =  all_snps_summary_stats$pval[ , -c(1:2) ] %>%
#     is.na %>%
#     rowSums %>%
#     ( function(x) x < ncol( all_snps_summary_stats$pval ) - 2 ) %>%
#     lapply( all_snps_summary_stats, filter, . )
# 
# traits_clumped_snps  =  all_snps_summary_stats$pval %>%
#     '['( apply( .[ , -c(1:2) ],
#                 1,
#                 min, na.rm = TRUE ) < SNP_THRESHOLD, ) %>%
#     get_ranks( phenotype      = 'body',
#                rank_threshold = SNP_THRESHOLD ) %>%
#     clump_genome
# 
# correlation_matrix  =  read_rds( 'data/processed/both_sexes/body/correlation_matrix.rds' )
# 
# category_pcs  =
#     get_category_pca_stats( all_snps_summary_stats,
#                             exposure_phenotypes,
#                             'body',
#                             traits_clumped_snps,
#                             line_up_positive = 'whr',
#                             correlation_matrix = correlation_matrix,
#                             adjust = TRUE,
#                             sex = 'both_sexes',
#                             npcs = NPCS )
# 
# # Lining up PCs to match the manuscript orientation
# category_pcs$pca$rotation[,4] = -category_pcs$pca$rotation[,4]
# 
# full_category_pcs  =  calculate_pc_from_loadings( all_snps_summary_stats,
#                                                   category_pcs$pca,
#                                                   npcs = 4 )
# 
# full_category_pcs[ c( 'beta', 'se', 'pval' ) ]  =  full_category_pcs[ c( 'beta', 'se', 'pval' ) ] %>%
#     lapply( function( x ) as_tibble( x ) %>%
#                 bind_cols( all_snps_summary_stats$beta[ , c( 'rsid', 'variant' ) ], . ) )
# 
# reclumped_pcs  =  clump_each( full_category_pcs,
#                               threshold = SNP_THRESHOLD )
# write_rds( reclumped_pcs, 'data/both_sexes/body_pcs.rds' )

# Saving SNP lists
pc_snps  =  reclumped_pcs[ names( reclumped_pcs ) != 'pca' ] %>%
    lapply( function(x) x$variant ) %>%
    unlist %>%
    unique %>%
    tibble( variant = . ) %>%
    left_join( COMMON_SNPS )

pcs_reclumped_snps  =  reclumped_pcs[ names( reclumped_pcs ) != 'pca' ] %>%
    lapply( function( pc ) {
        pc %>%
            select( variant ) %>%
            left_join( pc_snps )
    } )
write_rds( pcs_reclumped_snps,
           'data/clumped_snps/both_sexes_body_pcs_clumped_snps.rds' )

# Repeating SNP selection for disease (diabetes) as exposure
all_snps_summary_stats  =  get_variables_from_all( 'disease',
                                                   all_phenotypes = outcome_phenotypes,
                                                   variables = c( 'beta', 'se', 'pval' ),
                                                   sex = 'both_sexes',
                                                   threshold = FALSE )

traits_reclumped_snps  =  all_snps_summary_stats %>%
    clump_each( threshold = SNP_THRESHOLD ) %>%
    lapply( function( trait ) {
        trait %>%
            select( variant ) %>%
            left_join( COMMON_SNPS )
    } )
write_rds( traits_reclumped_snps,
           'data/clumped_snps/both_sexes_disease_traits_clumped_snps.rds' )


# Combine stats ####

categories  =  c( 'body', 'disease' )
sexes  =  c( 'female', 'male' )

all_snps  =  list.files( 'data/clumped_snps',
                         pattern = 'both_sexes_(body|disease)_(pcs|traits)_clumped_snps[.]rds',
                         full.names = TRUE ) %>%
    lapply( read_rds ) %>%
    '['( lapply( ., length ) > 0 ) %>%
    lapply( function( clumped_snps ) {
        clumped_snps %>%
            '['( lapply( ., length ) > 0 ) %>%
            lapply( function( x ) x$variant ) %>%
            unlist
    } ) %>%
    unlist %>%
    unique

for (sex in sexes) {
    all_snps_summary_stats  =  bind_rows( exposure_phenotypes,
                                          outcome_phenotypes ) %>%
        get_variables_from_all( category_name = categories,
                                all_phenotypes = .,
                                variables = c( 'beta', 'se', 'pval' ),
                                variant_list = all_snps,
                                sex = sex,
                                threshold = FALSE )
    
    write_rds( all_snps_summary_stats,
               sprintf( 'data/all_snps_summary_stats_%s.rds',
                        sex ) )
}

# Cross-sex MR analysis ####

source( 'mr_analysis.R' )

for (sex_exposure in sexes) {
    sex_outcome = setdiff( sexes, sex_exposure )
    
    all_exposure_stats  =  sprintf( 'data/all_snps_summary_stats_%s.rds',
                                    sex_exposure ) %>%
        read_rds
    
    all_outcome_stats  =  sprintf( 'data/all_snps_summary_stats_%s.rds',
                                   sex_outcome ) %>%
        read_rds
    
    exposure_pcs_ref  =  read_rds( 'data/both_sexes/body_pcs.rds' )
    
    exposure_pcs_ref$pca$correlation_matrix  =  sprintf( 'data/processed/%s/body/correlation_matrix.rds',
                                                         sex_exposure ) %>%
        read_rds
    
    pc_snps  =  read_rds( 'data/clumped_snps/both_sexes_body_pcs_clumped_snps.rds' )
    
    exposure_pcs  =  .select_each_pc_stats( all_exposure_stats,
                                            exposure_phenotypes,
                                            'body',
                                            pc_snps,
                                            exposure_pcs_ref$pca )
    
    dir_name  =  sprintf( 'results/%s_exp_%s_out_b/',
                          sex_exposure,
                          sex_outcome )
    dir.create( dir_name,
                showWarnings = FALSE,
                recursive    = TRUE )
    
    all_phenotypes  =  bind_rows( exposure_phenotypes, outcome_phenotypes )
    mr_results  =  mr_each_pc_vs_trait( all_outcome_stats,
                                        all_phenotypes      = all_phenotypes,
                                        exposure_pcs        = exposure_pcs,
                                        exposure_snps       = pc_snps,
                                        outcome_category    = 'disease',
                                        threshold           = SNP_THRESHOLD,
                                        het_threshold       = HETEROGENEITY,
                                        min_ivs_before      = MIN_IVS_BEFORE,
                                        min_ivs_after       = MIN_IVS_AFTER,
                                        reverse_t_threshold = REVERSE_T_THRESHOLD )
    write_rds( mr_results,
               sprintf( '%s/body_disease_pvt_mr.rds',
                        dir_name ) )
    
    dir_name  =  sprintf( 'results/%s_exp_%s_out_b/',
                          sex_outcome,
                          sex_exposure )
    dir.create( dir_name,
                showWarnings = FALSE,
                recursive    = TRUE )
    
    mr_results  =  mr_each_trait_vs_pc( all_outcome_stats,
                                        all_phenotypes      = all_phenotypes,
                                        exposure_category   = 'disease',
                                        exposure_snps       = read_rds('data/clumped_snps/both_sexes_disease_traits_clumped_snps.rds'),
                                        outcome_pca         = exposure_pcs_ref,
                                        outcome_category    = 'body',
                                        all_outcome_stats   = all_exposure_stats,
                                        threshold           = SNP_THRESHOLD,
                                        het_threshold       = HETEROGENEITY,
                                        min_ivs_before      = MIN_IVS_BEFORE,
                                        min_ivs_after       = MIN_IVS_AFTER,
                                        reverse_t_threshold = REVERSE_T_THRESHOLD )
    write_rds( mr_results,
               sprintf( '%s/disease_body_tvp_mr.rds',
                        dir_name ) )
}

# Meta-analysis of sex-specific results ####
require( meta )

source( 'meta_analysis.R' )
dir.create( 'results/both_sexes',
            showWarnings = FALSE,
            recursive    = TRUE )

list.files( 'results/female_exp_male_out_b',
            pattern = '[tp]v[tp]_mr.rds$' ) %>%
    walk( function( filename ){
        mr_results_fm  =  sprintf( 'results/female_exp_male_out_b/%s',
                                   filename ) %>%
            read_rds
        mr_results_mf  =  sprintf( 'results/male_exp_female_out_b/%s',
                                   filename ) %>%
            read_rds
        meta_analyze( mr_results_fm,
                      mr_results_mf,
                      'female_exp_male_out',
                      'male_exp_female_out' ) %>%
            write_rds( sprintf( 'results/both_sexes/%s',
                                filename ) )
    } )

# Print (forward) results ####
read_rds( 'results/both_sexes/body_disease_pvt_mr.rds' )[ c('b', 'se', 'pval') ] %>%
    print



































