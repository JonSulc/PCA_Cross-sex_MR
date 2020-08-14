suppressMessages( require( tidyverse ) )
require( gridExtra )
require( RColorBrewer )
require( ggplot2 )
require( plotly )
require( reshape2 )
require( scales )
require( igraph )
require( deming )

PREFIX = '.'

THRESHOLD  =  5e-2
GW_THRESHOLD  =  5e-8
BONFERRONI_CORRECTION  =  TRUE
MF_PREVALENCE_THRESHOLD  =  .01
PATH  =  'thr5e-08het0.001rev-1.644854b10'
DATA_PATH  =  'data'
FUMA_DATA_PATH  =  '%s/fuma/results'
GTEX_TISSUES  =  'fuma/gtex_descriptions.rds'

COLUMN_NAMES = c( 'PC1',
                  'Weight',
                  'PC2',
                  'BMI',
                  'WHR',
                  'PC3',
                  'PC4' )

DEPICT_DESCRIPTION = 'pascal/depict/merged_height_allMAF_nonSynSplice_RestrictedGenes_hg19_final_WithClusterAnnotations_genesetenrichment_4-19-16_v2'
PROCESSED_DESCRIPTION = 'pascal/depict/processed_descriptions.rds'
PASCAL_PATH = 'pascal/'

# Filename blanks: path, sex_pca initial, 'pc' or '', pc number or trait, sex
PASCAL_FILENAME = '%s/body_%s%s%s_%s.PathwaySet--depict_discretized_cutoff3.2--sum.txt'

if (!file.exists( PROCESSED_DESCRIPTION )) {
    sort_depict_groups  =  function( description ) {
        if (length( description ) > 1) {
            return( sapply( description, sort_depict_groups ) )
        }
        description  =  str_to_lower( description )
        if (str_detect( description, 'short|bone|osteo' )) {
            return( 'Bone' )
        }
        if (str_detect( description, 'morpho|embryo|exencephaly|organ|lethal|patterning|(?:eye|system|tube) development|cleft|regionalization|transformation|otic' )) {
            return( 'embryo-/morphogenesis' )
        }
        if (str_detect( description, 'hormone|testosterone|sexual|corpus luteum' )) {
            return( 'hormonal/sexual characteristics' )
        }
        if (str_detect( description,
                        'axon|synap|learning|hippocampus|neuro|glutamate receptor|behavior|receptor|ataxia|cerebellum|myelin' )) {
            return( 'brain/nerves' )
        }
        if (str_detect( description, 'lip[io]|insulin|glucose|obesity' )) {
            return( 'adipose/energy homeostasis' )
        }
        if (str_detect( description, 'erythro|arter|heart|blood|pr interval|myeloid' )) {
            return( 'blood/heart' )
        }
        return( 'other' )
    }
    descriptions  =  read_tsv( DEPICT_DESCRIPTION,
                               col_types  =  cols_only( ID = 'c',
                                                        Desc = 'c',
                                                        MetaCluster = 'c',
                                                        MetaCluster_Desc = 'c' ) ) %>%
        rename( Name = ID ) %>%
        filter( str_detect(Name, '^MP:') ) %>%
        mutate( group = sort_depict_groups( MetaCluster_Desc ) )
    write_rds( descriptions, PROCESSED_DESCRIPTION )
}

sex_name  =  function( sex_exposure,
                        sex_ivs = sex_exposure,
                        sex_pca = sex_ivs,
                        sex_outcome = sex_exposure ){
    if (sex_exposure != sex_outcome) {
        sex  =  sex_ivs %>%
            str_sub( 1, 1 ) %>%
            sprintf( '%s_exp_%s_out_%s',
                     sex_exposure,
                     sex_outcome,
                     . )
    } else if (sex_exposure != sex_ivs){
        sex  =  sex_ivs %>%
            str_sub( 1, 1 ) %>%
            paste0( sex_exposure, '_', . )
    } else if (sex_pca != sex_ivs) {
        return( sprintf( '%s_%sivs_%spca',
                         sex_exposure,
                         str_sub( sex_ivs, 1, 1 ),
                         str_sub( sex_pca, 1, 1 ) ) )
    } else {
        return( sex_exposure )
    }
    
    if (sex_pca != sex_ivs) {
        return( sprintf( '%sivs_%spca',
                         sex,
                         str_sub( sex_pca, 1, 1 ) ) )
    }
    sex
}


add_a_transform  =  function( transforms_list,
                              groups,
                              suffixes,
                              new_marker_additions,
                              base_marker_additions = list() ){
    
    if (is.null( names( suffixes ))){
        suffixes  =  setNames( c(suffixes), suffixes )
        new_marker_additions  =  setNames( list( new_marker_additions ), suffixes )
    }
    transforms_list[[1]]$groups  =  groups
    transforms_styles  =  list()
    for (suffix in suffixes) {
        transforms_styles  =  transforms_list[[1]]$styles %>%
            lapply( function( x ){
                x$target  =  paste0( x$target, suffix )
                x$value$marker   =  c( new_marker_additions[[ suffix ]],
                                       x$value$marker )
                x
            } ) %>%
            c( transforms_styles )
    }
    if (length( base_marker_additions ) > 0) {
        transforms_list[[1]]$styles  =  c(
            lapply( transforms_list[[1]]$styles, function( x ){
                x$value$marker   =  c( base_marker_additions,
                                       x$value$marker )
                x
            } ),
            transforms_styles
        )
    } else {
        transforms_list[[1]]$styles  =  transforms_styles
    }
    
    transforms_list
}


line_coordinates  =  function( x_range,
                               y_range,
                               x_range_2 = NULL,
                               y_range_2 = NULL,
                               slope = 1,
                               intercept = 0,
                               include_0 = TRUE ) {
    
    x_range  =  range( x_range, x_range_2 )
    x_range  =  ((x_range[2] - x_range[1])/20) %>%
        c( -., . ) %>%
        '+'( x_range )
    y_range  =  range( y_range, y_range_2 )
    y_range  =  ((y_range[2] - y_range[1])/20) %>%
        c( -., . ) %>%
        '+'( y_range )
    
    if (slope > 0) {
        x  =  c( x_range[ 2 ],
                 (y_range[ 2 ] - intercept) / slope ) %>%
            min
        x  =  c( x_range[ 1 ],
                 (y_range[ 1 ] - intercept) / slope ) %>%
            max %>%
            c( x )
        y  =  c( x_range[ 2 ] * slope + intercept,
                 y_range[ 2 ] ) %>%
            min
        y  =  c( x_range[ 1 ] * slope + intercept,
                 y_range[ 1 ] ) %>%
            max %>%
            c( y )
    } else {
        x  =  c( x_range[ 2 ],
                 (y_range[ 1 ] - intercept) / slope ) %>%
            min
        x  =  c( x_range[ 1 ],
                 (y_range[ 2 ] - intercept) / slope ) %>%
            max %>%
            c( x )
        y  =  c( x_range[ 2 ] * slope + intercept,
                 y_range[ 1 ] ) %>%
            max
        y  =  c( x_range[ 1 ] * slope + intercept,
                 y_range[ 2 ] ) %>%
            min %>%
            c( y )
    }
    if (include_0) {
        if (prod(x) > 0) {
            x[ which.min( abs(x) ) ] = 0
        }
        if (prod(y) > 0) {
            y[ which.min( abs(y) ) ] = 0
        }
    }
    
    list( x = x, y = y )
}


interactive_mr_plot  =  function( mrxb,
                                  mrxp,
                                  mryb,
                                  mryp,
                                  outcome_descriptions,
                                  x_label,
                                  y_label,
                                  main_title,
                                  threshold_x = THRESHOLD,
                                  threshold_y = THRESHOLD,
                                  threshold_xy = THRESHOLD,
                                  annot = list(),
                                  separate_xy = FALSE,
                                  mrxs = NULL,
                                  mrys = NULL,
                                  show_subcategories = FALSE,
                                  showlegend = TRUE,
                                  selected = c(),
                                  selected_text = NULL,
                                  selected_group = c(),
                                  flip_line = FALSE,
                                  export_image = FALSE,
                                  add_regression_line = FALSE,
                                  leave_labels = FALSE,
                                  error_bar_size = qnorm( 0.975 )){
    data_to_plot  =  tibble( mrxb          = mrxb,
                             mrxs          = mrxs,
                             mrxp          = mrxp,
                             mryb          = mryb,
                             mrys          = mrys,
                             mryp          = mryp,
                             x_significant = ( mrxp < threshold_x ) + 1,
                             y_significant = ( mryp < threshold_y ) + 1,
                             separate_p    = 2*pnorm( -abs( mrxb - mryb ) / sqrt( mrxs^2 + mrys^2 ) ),
                             separate_sig  = ( separate_p < threshold_xy
                                               & x_significant + y_significant > 2 ) + 1,
                             separate_neg_p = 2*pnorm( -abs( mrxb + mryb ) / sqrt( mrxs^2 + mrys^2 ) ),
                             separate_neg_sig = ( separate_neg_p < threshold_xy
                                                  & x_significant + y_significant > 2 ) + 1,
                             group         = c( 'Not significant', x_label, y_label, 'Both' )[ x_significant + 2*y_significant - 2],
                             subcategory   = outcome_descriptions$group,
                             key           = names( mrxb ),
                             toggle        = key %in% selected )# %>%
    # mutate( mrxe = mrxs * (-qnorm( threshold_x/2 )),
    #         mrye = mrys * (-qnorm( threshold_y/2 )) )
    
    # data_to_plot <<- data_to_plot
    
    if (!is.null( outcome_descriptions$description )) {
        data_to_plot  =  data_to_plot %>%
            mutate_at( c( 'mrxp', 'mryp', 'separate_p', 'separate_neg_p' ),
                       function( x ) formatC( x, format = 'e', digits = 2 ) ) %>%
            mutate( xb = formatC( mrxb, digits = 2 ),
                    yb = formatC( mryb, digits = 2 ),
                    phenotype_group = outcome_descriptions$group )
        if (flip_line) {
            data_to_plot  =  data_to_plot %>%
                mutate( description = sprintf( '%s
%s b = %s, p = %s
%s b = %s, p = %s
p_diff = %s',
                                               outcome_descriptions$description %>% str_to_title,
                                               x_label, xb, mrxp,
                                               y_label, yb, mryp,
                                               separate_neg_p ) )
        } else {
            data_to_plot  =  data_to_plot  %>%
                mutate( description = sprintf( '%s
%s b = %s, p = %s
%s b = %s, p = %s
p_diff = %s',
                                               outcome_descriptions$description %>% str_to_title,
                                               x_label, xb, mrxp,
                                               y_label, yb, mryp,
                                               separate_p ) )
        }
    } else {
        data_to_plot$description  =  outcome_descriptions$hover_text
    }
    
    data_to_plot  =  data_to_plot %>%
        filter( !is.na( mrxb ) & !is.na( mryb ) )
    
    if (show_subcategories) {
        data_to_plot  =  data_to_plot %>%
            mutate( group = data_to_plot$subcategory %>% str_to_sentence %>% as.factor )
        n_groups  =  data_to_plot$group %>% nlevels
        group_palette = brewer_pal( 'qual',
                                    ifelse( TRUE, #n_groups > 9 | export_image,
                                            'Paired',
                                            'Set1' ) )( n_groups )
        # 'Paired' )( n_groups )
        
        names( group_palette )  =  levels( data_to_plot$group )
        group_styles  =  data_to_plot$group %>%
            levels %>%
            lapply( function( x ) list(
                target = x,
                name   = x,
                value  =  list( marker = list( color = group_palette[ x ] %>% unname ) )
            ) )
        transforms_list  =  list(
            list(
                type = 'groupby',
                groups = data_to_plot$group,
                styles = group_styles
            )
        )
        data_to_plot  =  data_to_plot %>%
            mutate( group = paste0( group, c( '', ' (x)', ' (y)', ' (b)' )[ x_significant + 2*( y_significant - 1 ) ] ) )
        transforms_list  =  transforms_list %>%
            add_a_transform( data_to_plot$group,
                             c( ' (x)', ' (y)', ' (b)' ) %>% setNames( ., . ),
                             list( ' (x)' = list( symbol = 'triangle-down',
                                                  size   = 10 ),
                                   ' (y)' = list( symbol = 'triangle-left',
                                                  size   = 10 ),
                                   ' (b)' = list( symbol = 'circle',
                                                  size   = 10 ) ),
                             base_marker_additions = list( opacity =.42 ) )
    } else {
        transforms_list  =  list(
            list(
                type = 'groupby',
                groups = data_to_plot$group,
                styles = list(
                    list( target = 'Not significant',
                          name   = 'Not significant',
                          value = list( marker = list( color  = '#848484',
                                                       symbol = 'circle-open' ) ) ),
                    list( target = x_label,
                          name   = x_label,
                          value = list( marker = list( color  = 'blue',
                                                       symbol = 'triangle-down' ) ) ),
                    list( target = y_label,
                          name   = y_label,
                          value = list( marker = list( color  = 'red',
                                                       symbol = 'triangle-left' ) ) ),
                    list( target = 'Both',
                          name   = 'Both',
                          value = list( marker = list( color  = 'black',
                                                       symbol = 'circle' ) ) )
                )
            )
        )
    }
    
    if (separate_xy) {
        if (flip_line) {
            data_to_plot  =  data_to_plot %>%
                mutate( group = paste0( group, c( '', ' (s)' )[ separate_neg_sig ] ) )
        } else {
            data_to_plot  =  data_to_plot %>%
                mutate( group = paste0( group, c( '', ' (s)' )[ separate_sig ] ) )
        }
        transforms_list  =  add_a_transform( transforms_list,
                                             data_to_plot$group,
                                             ' (s)',
                                             new_marker_additions = list( opacity = 1,
                                                                          size    = 10 ),
                                             base_marker_additions = list( opacity = .42,
                                                                           size    = 5 ) )
    }
    xrange = range( data_to_plot$mrxb )
    xrange = c( min( xrange[ 1 ], 0 ),
                max( xrange[ 2 ], 0 ) )
    xsize  = xrange[ 2 ] - xrange[ 1 ]
    xrange = xrange + c( -xsize, xsize )/20
    yrange = range( data_to_plot$mryb )
    yrange = c( min( yrange[ 1 ], 0 ),
                max( yrange[ 2 ], 0 ) )
    ysize  = yrange[ 2 ] - yrange[ 1 ]
    yrange = yrange + c( -ysize, ysize )/20
    
    if (export_image) {
        width = 1000 + 150*showlegend
        height = 800
        if (!leave_labels) {
            x_label = ''
            y_label = ''
            main_title = NULL
        }
        
        # if (showlegend) {
        #     m  =  list( l   = 10,
        #                 r   = 250,
        #                 b   = 10,
        #                 t   = 10,
        #                 pad = 0 )
        # }
    } else {
        width = height = NULL
    }
    
    identity  =  line_coordinates( data_to_plot$mrxb,
                                   data_to_plot$mryb,
                                   slope = ifelse( flip_line,
                                                   -1,
                                                   1 ) )
    
    all_lines  =  list( list( type = 'line',
                              line = list( color = 'lightgray',
                                           dash  = 'dash' ),
                              x0 = identity$x[ 1 ],
                              x1 = identity$x[ 2 ],
                              y0 = identity$y[ 1 ],
                              y1 = identity$y[ 2 ] ) )
    
    
    
    # Accepts RGB or RGBA
    error_bar_color  =  '#00000011'
    
    m = list( t = 50 )
    
    p  =  plot_ly( data_to_plot, type = 'scatter', mode = 'none',
                   width = width,
                   height = height
    ) %>%
        add_markers( x = ~mrxb, y = ~mryb, type = 'scatter', mode = 'markers',
                     name = ' ',
                     error_x = ~list( array = mrxs * error_bar_size,
                                      color = error_bar_color ),
                     error_y = ~list( array = mrys * error_bar_size,
                                      color = error_bar_color ),
                     transforms = transforms_list,
                     hoverinfo = 'text',
                     text = ~description ) %>%
        plotly::layout( title = main_title,
                        xaxis  = list( title = x_label,
                                       range = xrange ),
                        yaxis  = list( title = y_label,
                                       range = yrange ),
                        # shapes = all_lines,
                        margin = m,
                        showlegend = showlegend,
                        annotations = annot,
                        legend = list( title = list( text = '<b>Significant in</b>')))
    
    if (!is_empty(selected)) {
        if (is.null( selected_text )) {
            selected_text = outcome_descriptions$description[ data_to_plot$toggle ] %>%
                str_to_sentence
        }
        
        p  =  data_to_plot %>%
            filter( toggle ) %>%
            (function( data ){
                add_annotations( p,
                                 x = data$mrxb,
                                 y = data$mryb,
                                 text = selected_text,
                                 xref = "x",
                                 yref = "y",
                                 showarrow = TRUE,
                                 arrowhead = 0,
                                 arrowsize = 1,
                                 standoff = 5,
                                 ax = -20,
                                 ay = -20,
                                 # xanchor = 'left',
                                 font = list( size = 16 ))
            })
        
    }
    if (add_regression_line) {
        p  =  p %>%
            add_deming_regression_line( data_to_plot$mrxb,
                                        data_to_plot$mryb,
                                        data_to_plot$mrxs,
                                        data_to_plot$mrys,
                                        all_lines,
                                        flip_line = flip_line )
    } else {
        p  =  p %>%
            plotly::layout( shapes = all_lines )
    }
    if (is_empty( annot )) {
        p  =  p %>%
            add_annotations(
                x  =  data_to_plot$mrxb,
                y  =  data_to_plot$mryb,
                text = outcome_descriptions$description %>% str_to_sentence,
                visible = FALSE,
                xref = 'x',
                yref = 'y',
                showarrow = TRUE,
                arrowhead = 0,
                clicktoshow = 'onoff'
            )
    }
    p
}


add_deming_regression_line  =  function( p,
                                         x,
                                         y,
                                         xs,
                                         ys,
                                         all_lines,
                                         conf_int = TRUE,
                                         flip_line = FALSE,
                                         ... ) {
    slope  =  deming( y ~ x-1,
                      xstd = xs, ystd = ys,
                      jackknife = FALSE )$coefficients[ 2 ]
    slope_angles  =  c( 1:length(x) ) %>%
        sapply( function( index ){
            deming( y[-index] ~ x[-index]-1,
                    xstd = xs[-index], ystd = ys[-index],
                    jackknife = FALSE )$coefficients[ 2 ] %>%
                atan
        } )
    jackknife_se  =  (length(x) * atan(slope) - (length(x)-1) * slope_angles) %>%
        sd %>%
        '/'( sqrt( length(x) ) )
    slope_p  =  2*pnorm( -abs( atan(slope) / jackknife_se ) )
    slope_not_1_p  =  2*pnorm( -abs( ( atan(slope) - ifelse(flip_line, -pi/4, pi/4) ) / jackknife_se ) )
    
    regression_line  =  line_coordinates( x,
                                          y,
                                          slope = slope) #,
                                          # intercept = regression[ 'Intercept', 'EST' ] )
    all_lines  =  all_lines %>%
        c( list( list( type = 'line',
                       line = list( color = 'black',
                                    dash  = 'solid' ),
                       x0 = regression_line$x[ 1 ],
                       x1 = regression_line$x[ 2 ],
                       y0 = regression_line$y[ 1 ],
                       y1 = regression_line$y[ 2 ]) ) )
    
    if (conf_int) {
        upper  =  (atan( slope ) + qnorm( 0.975 ) * jackknife_se) %>%
            tan %>%
            line_coordinates( c(-1,1), c(-1,1), slope = . )
        lower  =  (atan( slope ) - qnorm( 0.975 ) * jackknife_se) %>%
            tan %>%
            line_coordinates( c(-1,1), c(-1,1), slope = . )
        all_lines  =  list( list( type = 'path',
                                  path = sprintf( 'M0 0 L%f %f L%f %f L0 0 L%f %f L%f %f Z',
                                                  upper$x[ 1 ],
                                                  upper$y[ 1 ],
                                                  lower$x[ 1 ],
                                                  lower$y[ 1 ],
                                                  upper$x[ 2 ],
                                                  upper$y[ 2 ],
                                                  lower$x[ 2 ],
                                                  lower$y[ 2 ] ),
                                  line = list( width = 0 ),
                                  fillcolor = 'rgba(128,128,128,.1)' ) ) %>%
            c( all_lines )
    }
    
    # p <<- p
    p %>%
        add_annotations(
            x= 0.05,
            y= 0.95,
            xref = 'paper',
            yref = 'paper',
            align = 'left',
            text = sprintf( "<b>Regression</b>
slope = %.3g
p<sub>slope %s 1</sub> = %.2g
<b>Correlation</b>
r = %.2g
p = %.2g",
# intercept = %.2g, p = %.2g",
                            slope, 
                            # slope_p,
                            # 2*pnorm( -abs( regression[ 'Slope', 'EST' ]/regression[ 'Slope', 'SE' ] ) ),
                            '\U2260',
                            slope_not_1_p,
                            cor( x, y ),
                            cor.test( x, y )$p.value ),
                            # 2*pnorm( -abs( (regression[ 'Slope', 'EST' ]-1)/regression[ 'Slope', 'SE' ] ) ),
                            # regression[ 'Intercept', 'EST' ], 
                            # 2*pnorm( -abs( regression[ 'Intercept', 'EST' ]/regression[ 'Intercept', 'SE' ] ) ) ),
            showarrow = F
        ) %>%
        plotly::layout( shapes = all_lines )
}


compare_2_exposures  =  function( exposure_category_1,
                                  exposure_category_2,
                                  outcome_category,
                                  sex_exposure  = 'both_sexes',
                                  sex_ivs  = sex_exposure,
                                  sex_pca  = sex_ivs,
                                  sex_outcome  = sex_exposure,
                                  sexy_exposure = sex_exposure,
                                  sexy_ivs = sex_ivs,
                                  sexy_pca = sexy_ivs,
                                  sexy_outcome = sexy_exposure,
                                  exposure_type_1 = 'p',
                                  exposure_type_2 = 'p',
                                  outcome_type = 't',
                                  exp_trait_1 = 1,
                                  exp_trait_2 = 2,
                                  x_path = PATH,
                                  x_data_path = DATA_PATH,
                                  y_path = PATH,
                                  y_data_path = DATA_PATH,
                                  threshold = THRESHOLD,
                                  bonferroni_correction = BONFERRONI_CORRECTION,
                                  all_phenotypes = NULL,
                                  x_label = NULL,
                                  y_label = NULL,
                                  main_title = NULL,
                                  leave_labels = !is.null(c( x_label,
                                                             y_label,
                                                             main_title )),
                                  ... ){
    
    if (is.null( all_phenotypes )) {
        all_phenotypes  =  sprintf( '%s/all_phenotypes_%s.rds',
                                    x_data_path,
                                    sex_exposure ) %>%
            read_rds %>%
            full_join( read_rds( sprintf( '%s/all_phenotypes_%s.rds',
                                          y_data_path,
                                          sexy_exposure ) ),
                       by = c( 'phenotype', 'category', 'group', 'description' ) )
    }
    
    sex  =  sex_name( sex_exposure,
                      sex_ivs = sex_ivs,
                      sex_pca = sex_pca,
                      sex_outcome = sex_outcome )
    sexy  =  sex_name( sexy_exposure,
                       sex_ivs = sexy_ivs,
                       sex_pca = sexy_pca,
                       sex_outcome = sexy_outcome )
    
    mrx  =  sprintf( '%s/%s/%s_%s_%sv%s_mr.rds',
                     x_path,
                     sex,
                     exposure_category_1,
                     outcome_category,
                     exposure_type_1,
                     outcome_type ) %>%
        read_rds %>%
        '['( c( 'b', 'se', 'pval' ) )
    
    mry  =  sprintf( '%s/%s/%s_%s_%sv%s_mr.rds',
                     y_path,
                     sexy,
                     exposure_category_2,
                     outcome_category,
                     exposure_type_2,
                     outcome_type ) %>%
        read_rds %>%
        '['( c( 'b', 'se', 'pval' ) )
    
    if (bonferroni_correction) {
        threshold_x  =  threshold / prod( dim( mrx$b ) )
        threshold_y  =  threshold / prod( dim( mry$b ) )
    } else {
        threshold_x = threshold_y = threshold
    }
    if (exposure_type_1 == 't')
        exp_trait_1  =  .convert_name( exp_trait_1, mrx, all_phenotypes )
    if (exposure_type_2 == 't')
        exp_trait_2  =  .convert_name( exp_trait_2, mry, all_phenotypes )
    
    common_outcomes  =  list( mrx$b, mry$b ) %>%
        lapply( colnames ) %>%
        do.call( intersect, . )
    
    mrx  =  mrx %>%
        lapply( function( x ) x[ exp_trait_1, common_outcomes ] )
    mry  =  mry %>%
        lapply( function( x ) x[ exp_trait_2, common_outcomes ] )
    
    if (outcome_type == 't') {
        outcome_descriptions  =  suppressMessages(
            tibble( phenotype = common_outcomes ) %>%
                left_join( all_phenotypes[ , c( 'phenotype', 'description', 'group' ) ] )
        )
    } else {
        outcome_descriptions  =  tibble( phenotype = common_outcomes,
                                         description = common_outcomes,
                                         group = 'PCs' )
    }
    
    if (exposure_type_1 == 'p'){
        names( exp_trait_1 )  =  paste( exposure_category_1 %>% str_to_title,
                                        rownames( mrx$b ) )
    }
    if (exposure_type_2 == 'p'){
        names( exp_trait_2 )  =  paste( exposure_category_2 %>% str_to_title,
                                        rownames( mry$b ) )
    }
    
    if (is.null( leave_labels )) {
        leave_labels = !(is.null( x_label ) & is.null( y_label ))
    }
    
    if (is.null( x_label )) {
        x_label  =  sprintf( '%s -> %s',
                             names( exp_trait_1 ),
                             outcome_category %>% str_to_title )
    }
    if (is.null( y_label )) {
        y_label  =  sprintf( '%s -> %s',
                             names( exp_trait_2 ),
                             outcome_category %>% str_to_title )
    }
    if (is.null( main_title )) {
        main_title  =  paste( 'Effects of', names( exp_trait_1 ),
                              '&', names( exp_trait_2 ), '<br>on', outcome_category )
    }
    
    interactive_mr_plot( unlist( mrx$b ),
                         unlist( mrx$p ),
                         unlist( mry$b ),
                         unlist( mry$p ),
                         outcome_descriptions,
                         x_label,
                         y_label,
                         main_title = main_title,
                         threshold_x = threshold_x,
                         threshold_y = threshold_y,
                         threshold_xy = threshold / sum( (mrx$p<threshold) + (mry$p<threshold),
                                                         na.rm = 'TRUE' ),
                         mrxs = unlist( mrx$se ),
                         mrys = unlist( mry$se ),
                         leave_labels = leave_labels,
                         ... )
}


compare_2_outcomes  =  function( exposure_category,
                                 outcome_category_1,
                                 outcome_category_2,
                                 sex_exposure = 'both_sexes',
                                 sex_ivs      = sex_exposure,
                                 sex_pca      = sex_ivs,
                                 sex_outcome  = sex_exposure,
                                 sexy_exposure = sex_exposure,
                                 sexy_ivs     = sex_ivs,
                                 sexy_pca     = sexy_ivs,
                                 sexy_outcome = sexy_exposure,
                                 exposure_type = 'p',
                                 outcome_type_1 = 't',
                                 outcome_type_2 = 't',
                                 out_trait_1  = 1,
                                 out_trait_2  = 2,
                                 x_path = PATH,
                                 x_data_path = DATA_PATH,
                                 y_path = PATH,
                                 y_data_path = DATA_PATH,
                                 threshold = THRESHOLD,
                                 bonferroni_correction = BONFERRONI_CORRECTION,
                                 all_phenotypes = NULL,
                                 ... ){
    
    if (is.null( all_phenotypes )) {
        all_phenotypes  =  sprintf( '%s/all_phenotypes_%s.rds',
                                    x_data_path,
                                    sex_exposure ) %>%
            read_rds %>%
            full_join( read_rds( sprintf( '%s/all_phenotypes_%s.rds',
                                          y_data_path,
                                          sexy_exposure ) ),
                       by = c( 'phenotype', 'category', 'group', 'description' ) )
    }
    
    sex  =  sex_name( sex_exposure,
                      sex_ivs = sex_ivs,
                      sex_pca = sex_pca,
                      sex_outcome = sex_outcome )
    sexy  =  sex_name( sexy_exposure,
                       sex_ivs = sexy_ivs,
                       sex_pca = sexy_pca,
                       sex_outcome = sexy_outcome )
    
    mrx  =  sprintf( '%s/%s/%s_%s_%sv%s_mr.rds',
                     x_path,
                     sex,
                     exposure_category,
                     outcome_category_1,
                     exposure_type,
                     outcome_type_1 ) %>%
        read_rds
    
    mry  =  sprintf( '%s/%s/%s_%s_%sv%s_mr.rds',
                     y_path,
                     sexy,
                     exposure_category,
                     outcome_category_2,
                     exposure_type,
                     outcome_type_2 ) %>%
        read_rds
    
    if (bonferroni_correction) {
        threshold_x  =  threshold / prod( dim( mrx$b ) )
        threshold_y  =  threshold / prod( dim( mry$b ) )
    } else {
        threshold_x = threshold_y = threshold
    }
    if (outcome_type_1 == 't')
        out_trait_1  =  .convert_name( out_trait_1,
                                       mrx,
                                       all_phenotypes,
                                       name_function = colnames )
    if (outcome_type_2 == 't')
        out_trait_2  =  .convert_name( out_trait_2,
                                       mry,
                                       all_phenotypes,
                                       name_function = colnames )
    
    common_exposures  =  list( mrx$b, mry$b ) %>%
        lapply( rownames ) %>%
        do.call( intersect, . )
    
    mrx  =  mrx %>%
        lapply( function( x ) x[ common_exposures, out_trait_1, drop = FALSE ] )
    mry  =  mry %>%
        lapply( function( x ) x[ common_exposures, out_trait_2, drop = FALSE ] )
    
    if (exposure_type == 't') {
        exposure_descriptions  =  suppressMessages(
            tibble( phenotype = common_exposures ) %>%
                left_join( all_phenotypes[ , c( 'phenotype', 'description', 'group' ) ] )
        )
    } else {
        exposure_descriptions  =  tibble( phenotype = common_exposures,
                                          description = common_exposures,
                                          group = 'PCs' )
    }
    
    if (outcome_type_1 == 'p')
        names( out_trait_1 )  =  paste( outcome_category_1, colnames( mrx$b ) )
    if (outcome_type_2 == 'p')
        names( out_trait_2 )  =  paste( outcome_category_2, colnames( mry$b ) )
    
    x_label  =  paste( exposure_category %>% str_to_sentence, '->', names( out_trait_1 ) )
    y_label  =  paste( exposure_category %>% str_to_sentence, '->', names( out_trait_2 ) )
    
    main_title  =  paste( 'Effects of', exposure_category,
                          '<br>on', names( out_trait_1 ), '&', names( out_trait_2 ) )
    
    interactive_mr_plot( unlist( mrx$b ),
                         unlist( mrx$p ),
                         unlist( mry$b ),
                         unlist( mry$p ),
                         exposure_descriptions,
                         x_label,
                         y_label,
                         main_title = main_title,
                         threshold_x = threshold_x,
                         threshold_y = threshold_y,
                         threshold_xy = threshold / sum( (mrx$p<threshold) + (mry$p<threshold),
                                                         na.rm = 'TRUE' ),
                         mrxs = unlist( mrx$se ),
                         mrys = unlist( mry$se ),
                         ... )
}

.convert_name  =  function( trait,
                            mr_results = NULL,
                            all_phenotypes = stop( 'Requires all_phenotypes' ),
                            pheno_names = NULL,
                            name_function = rownames ){
    if (is.numeric( trait )) {
        if (is.null( pheno_names ))
            pheno_names = mr_results$b %>%
                name_function
        trait  =  all_phenotypes %>%
            filter( phenotype == pheno_names[ trait ] )
        trait  =  setNames( trait$phenotype, trait$description )
    } else if (str_detect( all_phenotypes$phenotype, paste0( '^', trait, '(_irnt)?$' ) ) %>% any ) {
        trait  =  all_phenotypes %>%
            filter( str_detect( phenotype, paste0( '^', trait, '(_irnt)?$' ) ) )
        trait  =  setNames( trait$phenotype, trait$description )
    } else {
        names( trait )  =  trait
    }
    trait
}

barplot_pcs  =  function( category,
                          pc_number,
                          all_phenotypes = read_rds( paste0( data_path, '/all_phenotypes_combined.rds' ) ),
                          threshold = NULL,
                          data_path = DATA_PATH,
                          sex = 'both_sexes',
                          sex_pca = '',
                          show_sexes = c( 'both_sexes', 'female', 'male' ),
                          single = length( show_sexes ) == 1,
                          keep_order = FALSE ){
    if (category == 'brain') {
        show_sexes  =  sex  =  'both_sexes'
    }
    if (single) {
        if (!(sex %in% show_sexes)) {
            stop( 'Main sex not displayed' )
        }
        show_sexes  =  sex
    }
    if (length( pc_number ) > 1) {
        return( lapply( pc_number, function( x ) barplot_pcs( category, x, all_phenotypes, threshold ) ) %>%
                    do.call( subplot, . ) )
    }
    
    pcs  =  c( 'both_sexes', 'female', 'male' ) %>%
        setNames( nm = . ) %>%
        lapply( function( sex ){
            filename = sprintf( '%s/%s/%s_%spcs.rds',
                                data_path,
                                sex,
                                category,
                                ifelse( sex == sex_pca,
                                        '',
                                        str_sub( sex_pca, 1, 1 ) ) )
            if (file.exists( filename )) {
                filename %>%
                    read_rds %>%
                    '$'( 'pca' ) %>%
                    '$'( 'rotation' ) %>%
                    apply( 2, function( x ) x / sqrt( sum( x^2, na.rm = TRUE ) ) )
            } else {
                list()
            }
        } )
    
    numbers  =  rep( pc_number, 3 )
    names( numbers )  =  c( 'both_sexes', 'female', 'male' )
    
    pc_names  =  names( numbers ) %>%
        setNames( nm = . ) %>%
        lapply( function( sex ){
            pcs[[ sex ]] %>%
                colnames %>%
                str_detect( sprintf( '^PC[ ]?%s(?: |$)', numbers[ sex ] ) ) %>%
                '['( colnames( pcs[[ sex ]] ), . )
        } )
    
    # Change show_sexes to reflect availability
    pcs_available  =  pc_names %>%
        sapply( function( x ) !is.null( x ) & !identical( x, character(0) ) )
    if (all( !pcs_available )) {
        stop( 'PC non-significant' )
    }
    pc_names  =  pc_names[ pcs_available ]
    show_sexes  =  pcs_available[ pcs_available ] %>%
        names %>%
        intersect( show_sexes )
    
    p = pcs %>%
        sapply( rownames ) %>%
        unlist %>%
        c %>%
        unique %>%
        tibble( phenotype = . ) %>%
        left_join( all_phenotypes ) %>%
        select( phenotype, trait = description, group )
    
    for (sex in show_sexes) {
        p  =  p %>%
            left_join( tibble( phenotype = rownames( pcs[[ sex ]] ),
                               !!sex := pcs[[ sex ]][ , pc_names[[ sex ]] ] ) )
    }
    
    p  =  select( p, -phenotype )
    
    if (is.null( threshold )) {
        threshold  =  p[ , show_sexes ] %>%
            abs %>%
            max( na.rm = TRUE ) %>%
            '/'( 3 )
    }
    
    sex_color_scale  =  brewer.pal( 3, 'Set1' )
    sex_color_scale[3]  =  '#424242'
    names( sex_color_scale ) = c( 'female', 'male', 'both_sexes' )
    sex_color_scale  =  scale_fill_manual( name = 'variable', values = sex_color_scale )
    
    p  =  p %>%
        filter( rowSums( abs( .[ , show_sexes, drop = FALSE ] ) > threshold, na.rm = TRUE ) > 0 )
    if (!keep_order) {
        p  =  select( p, sex ) %>%
            unlist %>%
            abs %>%
            order( decreasing = TRUE ) %>%
            '['( p, ., )
    } 
    p  =  p %>%
        mutate( trait = factor( trait, levels = trait ) ) %>%
        dplyr::select( c( 'trait', 'group', show_sexes ) )
    
    showlegend  =  TRUE
    
    p %>%
        split( p$group %>% str_to_sentence ) %>%
        map2( names( . ), function( mydata, myname ) {
            plot_ly( data = mydata, x = ~droplevels( trait ), showlegend = showlegend ) %>%
                {
                    showlegend  <<-  FALSE
                    if ('female' %in% show_sexes)
                        add_trace( ., y = ~female, type = 'bar',
                                   name = paste( 'Female', pc_names[ 'female' ] ),
                                   marker = list( color = ifelse( single,
                                                                  '#424242',
                                                                  '#E41A1C' ) ),
                                   width = ifelse( single,
                                                   .8,
                                                   .2 ),
                                   opacity = ifelse( single,
                                                     1,
                                                     .5 ) )
                    else
                        .
                } %>%
                {
                    if ('both_sexes' %in% show_sexes)
                        add_trace( ., y = ~both_sexes, type = 'bar',
                                   name = paste( 'Both', pc_names[ 'both_sexes' ] ),
                                   marker = list( color = '#424242' ),
                                   width = ifelse( single,
                                                   .8,
                                                   .5 ) )
                    else
                        .
                } %>%
                {
                    if ('male' %in% show_sexes)
                        add_trace( ., y = ~male, type = 'bar',
                                   name = paste( 'Male', pc_names[ 'male' ] ),
                                   marker = list( color = ifelse( single,
                                                                  '#424242',
                                                                  '#377EB8' ) ),
                                   width = ifelse( single,
                                                   .8,
                                                   .2 ),
                                   opacity = ifelse( single,
                                                     1,
                                                     .5 ) )
                    else
                        .
                } %>%
                plotly::layout( annotations = list(
                    text = myname,
                    xref = "paper",
                    yref = "paper",
                    yanchor = "bottom",
                    xanchor = "center",
                    align = "center",
                    x = 0.5,
                    y = 1,
                    showarrow = FALSE),
                    xaxis = list( tickangle = 45 ) )
        } ) %>%
        subplot( margin = .01, shareY = TRUE,
                 widths = p$group %>% table %>% '/'( sum( . ) ) %>% unname ) %>%
        plotly::layout( yaxis = list( title = 'Loadings' ))
}


bidirectional_mr_plot  =  function( category_1,
                                    category_2,
                                    sex = 'both_sexes',
                                    sex_ivs = sex,
                                    sex_pca = sex_ivs,
                                    type_1 = 'p',
                                    type_2 = 't',
                                    trait = NULL,
                                    path = PATH,
                                    threshold = THRESHOLD,
                                    bonferroni_correction = BONFERRONI_CORRECTION,
                                    data_path = DATA_PATH,
                                    all_phenotypes = sprintf( '%s/all_phenotypes_%s.rds',
                                                              data_path,
                                                              sex ) %>%
                                        read_rds,
                                    ... ){
    
    sex_full = sex_name( sex,
                         sex_ivs = sex_ivs,
                         sex_pca = sex_pca )
    
    direction_1  =  sprintf( '%s/%s/%s_%s_%sv%s_mr.rds',
                             path,
                             sex_full,
                             category_1,
                             category_2,
                             type_1,
                             type_2 ) %>%
        read_rds %>%
        '['( c( 'b', 'se', 'pval' ) )
    
    direction_2  =  sprintf( '%s/%s/%s_%s_%sv%s_mr.rds',
                             path,
                             sex_full,
                             category_2,
                             category_1,
                             type_2,
                             type_1 ) %>%
        read_rds %>%
        '['( c( 'b', 'se', 'pval' ) )
    
    if (bonferroni_correction) {
        threshold_x  =  threshold / ( dim( direction_1$b ) %>% prod )
        threshold_y  =  threshold / ( dim( direction_2$b ) %>% prod )
    } else {
        threshold_x = threshold_y = threshold
    }
    
    if (type_1 == 'p') {
        common_traits = rownames( direction_1$pval ) %>%
            str_match( '^PC[ ]?([0-9]+)' ) %>%
            '['( , 2 ) %>%
            intersect( str_match( colnames( direction_2$pval ),
                                  '^PC[ ]?([0-9]+)' )[ , 2 ] )
        to_keep  =  common_traits %>%
            paste( collapse = '|' ) %>%
            sprintf( '^PC[ ]?(%s)( [(0-9.%%)]+)?$', . ) %>%
            str_detect( colnames( direction_2$pval ), . )
        direction_2  =  direction_2 %>%
            lapply( function(x) x[ , to_keep ] )
    } else {
        common_traits  =  intersect( rownames( direction_1$pval ),
                                     colnames( direction_2$pval ) )
        direction_1  =  direction_1 %>%
            lapply( function(x) x[ common_traits,  ] )
        direction_2  =  direction_2 %>%
            lapply( function(x) x[ , common_traits ] )
    }
    
    if (type_2 == 'p') {
        common_traits = rownames( direction_2$pval ) %>%
            str_match( '^PC[ ]?([0-9]+)' ) %>%
            '['( , 2 ) %>%
            intersect( str_match( colnames( direction_1$pval ),
                                  '^PC[ ]?([0-9]+)' )[ , 2 ] )
        to_keep  =  common_traits %>%
            paste( collapse = '|' ) %>%
            sprintf( '^PC[ ]?(%s)( [(0-9.%%)]+)?$', . ) %>%
            str_detect( colnames( direction_1$pval ), . )
        direction_1  =  direction_1 %>%
            lapply( function(x) x[ , to_keep ] )
    } else {
        common_traits  =  intersect( rownames( direction_2$pval ),
                                     colnames( direction_1$pval ) )
        direction_1  =  direction_1 %>%
            lapply( function(x) x[ , common_traits ] )
        direction_2  =  direction_2 %>%
            lapply( function(x) x[ common_traits,  ] )
    }
    
    if (!is.null( trait )) {
        if(type_1 == 't'){
            trait   =  .convert_name( trait, direction_1, all_phenotypes )
            direction_1  =  direction_1 %>%
                lapply( function( x ) x[ trait, , drop = FALSE ] )
            direction_2  =  direction_2 %>%
                lapply( function( x ) x[ , trait, drop = FALSE  ] )
            rownames( direction_1$b )  =  trait
        } else {
            names( trait )  =  rownames( direction_1$b )[ trait ]
            direction_1  =  direction_1 %>%
                lapply( function( x ) x[ trait, , drop = FALSE ] )
            direction_2  =  direction_2 %>%
                lapply( function( x ) x[ , trait, drop = FALSE  ] )
        }
    }
    
    mrx  =  direction_1[ c( 'b', 'pval', 'se' ) ] %>%
        lapply( function( x ) x %>% unlist %>% as.numeric )
    
    mry  =  direction_2[ c( 'b', 'pval', 'se' ) ] %>%
        lapply( function( x ) x %>% t %>% as.data.frame %>% unlist %>% as.numeric )
    
    
    if (type_2 == 't') {
        outcomes  =  suppressMessages(
            colnames( direction_1$b ) %>%
                tibble( phenotype = . ) %>%
                left_join( all_phenotypes[ , c( 'phenotype', 'description', 'group' ) ] )
        )
        colnames( direction_1$b )  =  outcomes$description
        outcome_groups  =  outcomes$group
    } else {
        outcomes  =  tibble( group = 'PCs' )
    }
    
    trait_1  =  rep( rownames( direction_1$b ), times = ncol( direction_1$b ) )
    trait_2  =  rep( colnames( direction_1$b ), each  = nrow( direction_1$b ) )
    
    if (!is.null( trait )) {
        main_title  =  paste0( 'Bidirectional MR<br>', names( trait ), ' vs ', category_2 )
        x_label  =  paste( names( trait ), '->', category_2 )
        y_label  =  paste( category_2,     '->', names( trait ) )
        outcome_descriptions  =  paste0( trait_2, '<br>',
                                         '  v p = ', formatC( mrx$pval, format = 'e', digits = 2 ), '<br>',
                                         '  < p = ', formatC( mry$pval, format = 'e', digits = 2 ) ) %>%
            tibble( hover_text = .,
                    group = outcomes$group,
                    description = trait_2 )
    } else {
        main_title  =  paste0( 'Bidirectional MR<br>', category_1, ' vs ', category_2 )
        x_label  =  paste( category_1, '->', category_2 )
        y_label  =  paste( category_2, '->', category_1 )
        outcome_descriptions  =  paste0( trait_1, ' - ', trait_2, '<br>',
                                         '  v p = ', formatC( mrx$pval, format = 'e', digits = 2 ), '<br>',
                                         '  < p = ', formatC( mry$pval, format = 'e', digits = 2 ) ) %>%
            tibble( hover_text = .,
                    group = outcomes$group,
                    description = paste( trait_1, '-', trait_2 ) )
    }
    
    interactive_mr_plot( mrxb = setNames( mrx$b, trait ),
                         mrxp = mrx$pval,
                         mrxs = mrx$se,
                         mryb = mry$b,
                         mryp = mry$pval,
                         mrys = mry$se,
                         outcome_descriptions = outcome_descriptions,
                         x_label = x_label,
                         y_label = y_label,
                         main_title = main_title,
                         threshold_x = threshold_x,
                         threshold_y = threshold_y,
                         threshold_xy = threshold / sum( (mrx$pval<threshold) + (mry$pval<threshold),
                                                         na.rm = 'TRUE' ),
                         ... )
}


compare_instruments  =  function( exposure_category_1,
                                  exposure_category_2,
                                  outcome_category,
                                  sex = 'both_sexes',
                                  sex_ivs = sex,
                                  sex_pca = sex_ivs,
                                  sex_2 = sex,
                                  sex_ivs_2 = sex_2,
                                  sex_pca_2 = sex_ivs_2,
                                  exposure_type_1 = 'p',
                                  exposure_type_2 = 't',
                                  outcome_type = 't',
                                  sex_outcome = sex,
                                  sex_outcome_2 = sex_2,
                                  exp_trait_1 = 1,
                                  exp_trait_2 = 1,
                                  out_trait = 1,
                                  prefix_1 = '..',
                                  path_1 = list.files( prefix_1,
                                                       pattern = '^thr',
                                                       full.names = TRUE )[1],
                                  data_path_1 = paste0( prefix_1, '/data' ),
                                  prefix_2 = '..',
                                  path_2 = list.files( prefix_2,
                                                       pattern = '^thr',
                                                       full.names = TRUE )[1],
                                  data_path_2 = paste0( prefix_2, '/data' ),
                                  threshold = GW_THRESHOLD,
                                  all_phenotypes = NULL,
                                  show_mr_estimates = TRUE,
                                  ... ){
    
    if (is.null( all_phenotypes )) {
        all_phenotypes  =  sprintf( '%s/all_phenotypes_%s.rds',
                                    data_path_1,
                                    sex ) %>%
            read_rds %>%
            select( phenotype ) %>%
            full_join( read_rds( sprintf( '%s/all_phenotypes_%s.rds',
                                          data_path_2,
                                          sex ) ) )
    }
    
    if (exposure_type_1 == 'p') {
        type_name  =  'pcs'
    } else if (exposure_type_1 == 't') {
        type_name  =  'traits'
    } else {
        type_name  =  'dxa'
    }
    clumped_snps_1  =  sprintf( '%s/clumped_snps/%s_%s_%s%s_clumped_snps.rds',
                                data_path_1,
                                sex_ivs,
                                exposure_category_1,
                                ifelse( sex_ivs == sex_pca | exposure_type_1 != 'p', '', str_sub( sex_pca, 1, 1 ) ),
                                type_name ) %>%
        read_rds
    
    if (exposure_type_1 == 't') {
        exposure_1  =  sprintf( '%s/all_snps_summary_stats_%s.rds',
                                data_path_1,
                                sex ) %>%
            read_rds %>%
            .select_each_category_stats( all_phenotypes,
                                         exposure_category_1,
                                         clumped_snps_1[ exp_trait_1 ] ) %>%
            '[['( 1 ) %>%
            select( -one_of( c( 'variant', 'rsid' ) ) )
        exp_trait_1  =  .convert_name( exp_trait_1,
                                       pheno_names = names( clumped_snps_1 ),
                                       all_phenotypes = all_phenotypes )
    } else {
        exposure_1  =  sprintf( '%s/%s/%s_%spcs.rds',
                                data_path_1,
                                sex,
                                exposure_category_1,
                                ifelse( sex_ivs == sex_pca, '', str_sub( sex_pca, 1, 1 ) ) ) %>%
            read_rds %>%
            '[['( exp_trait_1 ) %>%
            select( -one_of( c( 'variant', 'rsid' ) ) )
        names( exp_trait_1 )  =  paste( exposure_category_1 %>% str_to_sentence,
                                        names( clumped_snps_1 )[ exp_trait_1 ] )
    }
    
    if (exposure_type_2 == 'p') {
        type_name  =  'pcs'
    } else if (exposure_type_2 == 't') {
        type_name  =  'traits'
    } else {
        type_name  =  'dxa'
    }
    clumped_snps_2  =  sprintf( '%s/clumped_snps/%s_%s_%s%s_clumped_snps.rds',
                                data_path_2,
                                sex_ivs_2,
                                exposure_category_2,
                                ifelse( sex_ivs == sex_pca | exposure_type_2 != 'p', '', str_sub( sex_pca, 1, 1 ) ),
                                type_name ) %>%
        read_rds
    if (exposure_type_2 == 't') {
        exposure_2  =  sprintf( '%s/all_snps_summary_stats_%s.rds',
                                data_path_2,
                                sex_2 ) %>%
            read_rds %>%
            .select_each_category_stats( all_phenotypes,
                                         exposure_category_2,
                                         clumped_snps_2[ exp_trait_2 ] ) %>%
            '[['( 1 ) %>%
            select( -one_of( c( 'variant', 'rsid' ) ) )
        exp_trait_2  =  .convert_name( exp_trait_2,
                                       pheno_names = names( clumped_snps_2 ),
                                       all_phenotypes = all_phenotypes )
    } else {
        exposure_2  =  sprintf( '%s/%s/%s_%spcs.rds',
                                data_path_2,
                                sex_2,
                                exposure_category_2,
                                ifelse( sex_ivs_2 == sex_pca_2,
                                        '',
                                        str_sub( sex_pca_2, 1, 1 ) ) ) %>%
            read_rds %>%
            '[['( exp_trait_2 ) %>%
            select( -one_of( c( 'variant', 'rsid' ) ) )
        names( exp_trait_2 )  =  paste( exposure_category_2 %>% str_to_sentence,
                                        names( clumped_snps_2 )[ exp_trait_2 ] )
    }
    
    if (outcome_type == 't') {
        out_trait    =  all_phenotypes %>%
            filter( category == outcome_category ) %>%
            .convert_name( out_trait,
                           pheno_names = .$phenotype,
                           all_phenotypes = . )
        
        outcome_1  =  sprintf( '%s/all_snps_summary_stats_%s.rds',
                               data_path_1,
                               sex ) %>%
            read_rds %>%
            .select_each_category_stats( all_phenotypes %>% filter( phenotype == out_trait ),
                                         outcome_category,
                                         setNames( clumped_snps_1[ exp_trait_1 ], out_trait ) ) %>%
            '[['( 1 ) %>%
            select( -one_of( c( 'variant', 'rsid' ) ) )
        
        outcome_2  =  sprintf( '%s/all_snps_summary_stats_%s.rds',
                               data_path_2,
                               sex ) %>%
            read_rds %>%
            .select_each_category_stats( all_phenotypes %>% filter( phenotype == out_trait ),
                                         outcome_category,
                                         setNames( clumped_snps_2[ exp_trait_2 ], out_trait ) ) %>%
            '[['( 1 ) %>%
            select( -one_of( c( 'variant', 'rsid' ) ) )
    } else {
        outcome_1_pc  =  sprintf( '%s/%s/%s_%spcs.rds',
                                  data_path_1,
                                  sex,
                                  outcome_category,
                                  ifelse( sex_ivs == sex_pca,
                                          '',
                                          str_sub( sex_pca, 1, 1 ) ) ) %>%
            read_rds
        if (is.integer( out_trait )) {
            out_name  =  (outcome_1_pc %>% names)[ out_trait ]
        } else {
            out_name  =  out_trait
        }
        names( out_trait )  =  out_name
        
        outcome_1  =  sprintf( '%s/all_snps_summary_stats_%s.rds',
                               data_path_1,
                               sex ) %>%
            read_rds %>%
            .select_each_pc_stats( all_phenotypes,
                                   outcome_category,
                                   clumped_snps_1[ exp_trait_1 ] %>% setNames( out_name ),
                                   outcome_1_pc$pca )
        
        outcome_2_pc  =  sprintf( '%s/%s/%s_%spcs.rds',
                                  data_path_2,
                                  sex_2,
                                  outcome_category,
                                  ifelse( sex_ivs_2 == sex_pca_2,
                                          '',
                                          str_sub( sex_pca_2, 1, 1 ) ) ) %>%
            read_rds
        
        outcome_2  =  sprintf( '%s/all_snps_summary_stats_%s.rds',
                               data_path_2,
                               sex_2 ) %>%
            read_rds %>%
            .select_each_pc_stats( all_phenotypes,
                                   outcome_category,
                                   clumped_snps_2[ exp_trait_2 ] %>% setNames( out_name ),
                                   outcome_2_pc$pca )
        
    }
    
    x_label  =  paste0( names( exp_trait_1 ), ', ', names( exp_trait_2 ) )
    y_label  =  names( out_trait )
    
    main_title  =  paste( 'Effect size comparison of', names( exp_trait_1 ),
                          '&', names( exp_trait_2 ), '<br>on', names( out_trait ) )
    
    data_to_plot  =  tibble( exp_1    = exposure_1$beta,
                             out_1    = outcome_1$beta,
                             exp_1_se = exposure_1$se,
                             out_1_se = outcome_1$se )
    data_to_plot_2  =  tibble( exp_2    = exposure_2$beta,
                               out_2    = outcome_2$beta,
                               exp_2_se = exposure_2$se,
                               out_2_se = outcome_2$se )
    
    if (show_mr_estimates) {
        mr_sex_1  =  sex_name( sex,
                               sex_ivs = sex_ivs,
                               sex_pca = sex_pca,
                               sex_outcome = sex_outcome )
        mr_1  =  sprintf( '%s/%s/%s_%s_%sv%s_mr.rds',
                          path_1,
                          mr_sex_1,
                          exposure_category_1,
                          outcome_category,
                          exposure_type_1,
                          outcome_type ) %>%
            read_rds
        line_type_1  =  (mr_1$pval[ exp_trait_1, out_trait ] < .05/prod( dim( mr_1$b ) )) %>%
            ifelse( 'solid', 'dot' )
        
        line_1  =  line_coordinates( x_range = exposure_1$beta,
                                     y_range = outcome_1$beta,
                                     x_range_2 = exposure_2$beta,
                                     y_range_2 = outcome_2$beta,
                                     slope = mr_1$b[ exp_trait_1, out_trait ] )
        
        mr_sex_2  =  sex_name( sex_2,
                               sex_ivs = sex_ivs_2,
                               sex_pca = sex_pca_2,
                               sex_outcome = sex_outcome_2 )
        
        mr_2  =  sprintf( '%s/%s/%s_%s_%sv%s_mr.rds',
                          path_2,
                          mr_sex_2,
                          exposure_category_2,
                          outcome_category,
                          exposure_type_2,
                          outcome_type ) %>%
            read_rds
        line_type_2  =  (mr_2$pval[ exp_trait_2, out_trait ] < .05/prod( dim( mr_2$b ) )) %>%
            ifelse( 'solid', 'dot' )
        
        line_2  =  line_coordinates( x_range = exposure_1$beta,
                                     y_range = outcome_1$beta,
                                     x_range_2 = exposure_2$beta,
                                     y_range_2 = outcome_2$beta,
                                     slope = mr_2$b[ exp_trait_2, out_trait ] )
    }
    
    plot_ly( data_to_plot, type = 'scatter', mode = 'none' ) %>%
        plotly::layout( title = main_title,
                        xaxis = list( title = x_label ),
                        yaxis = list( title = y_label ),
                        shapes = list( list( type = 'line',
                                             line = list( color = 'lightgray',
                                                          dash  = line_type_1 ),
                                             x0 = line_1$x[ 1 ],
                                             x1 = line_1$x[ 2 ],
                                             y0 = line_1$y[ 1 ],
                                             y1 = line_1$y[ 2 ] ),
                                       list( type = 'line',
                                             line = list( color = 'lightgray',
                                                          dash  = line_type_2 ),
                                             x0 = line_2$x[ 1 ],
                                             x1 = line_2$x[ 2 ],
                                             y0 = line_2$y[ 1 ],
                                             y1 = line_2$y[ 2 ] ) ),
                        showlegend = TRUE ) %>%
        add_trace( x = ~exp_1,
                   y = ~out_1,
                   mode = 'markers', name = names( exp_trait_1 ),
                   error_x = list(
                       type = 'data',
                       array = ~exp_1_se,
                       thickness = .5,
                       width = 0
                   ),
                   error_y = list(
                       type = 'data',
                       array = ~out_1_se,
                       thickness = .5,
                       width = 0
                   ) ) %>%
        add_trace( x = data_to_plot_2$exp_2,
                   y = data_to_plot_2$out_2,
                   mode = 'markers', name = names( exp_trait_2 ),
                   error_x = list(
                       type = 'data',
                       array = data_to_plot_2$exp_2_se,
                       thickness = .5,
                       width = 0
                   ),
                   error_y = list(
                       type = 'data',
                       array = data_to_plot_2$out_2_se,
                       thickness = .5,
                       width = 0
                   ) )
}

barplot_fuma  =  function( phenotype,
                           dataset,
                           threshold = NULL,
                           data_path = FUMA_DATA_PATH,
                           sex_pca = 'female',
                           magma_results_filename = sprintf( '%s/body_%s%s%s_%s/magma_exp_%s.gsa.out',
                                                             data_path,
                                                             '%s',
                                                             ifelse( is.numeric( phenotype ),
                                                                     'pc',
                                                                     '' ),
                                                             phenotype,
                                                             '%s',
                                                             dataset ),
                           tissue_descriptions = sprintf( '%s/%s.rds', data_path, dataset ) %>%
                               read_rds,
                           show_sexes = c( 'both_sexes', 'female', 'male' ),
                           single = length(show_sexes) == 1,
                           title_text = sprintf( ifelse( is.numeric( phenotype ),
                                                         'PC%s tissue enrichment',
                                                         '%s tissue enrichment' ),
                                                 phenotype ),
                           ordered = TRUE ){
    
    pcs  =  show_sexes %>%
        setNames( nm = . ) %>%
        lapply( function( sex ){
            filename  =  sprintf( magma_results_filename,
                                  ifelse( sex == sex_pca | !is.numeric( phenotype ),
                                          '',
                                          str_sub( sex_pca, 1, 1 ) ),
                                  sex )
            if (file.exists( filename )) {
                read_delim( filename,
                            delim = ' ',
                            comment = '#',
                            trim_ws = TRUE ) %>%
                    mutate( log10p = -log10( P ) )
            } else {
                NULL
            }
        } )
    
    show_sexes  =  show_sexes[ !sapply( pcs, is.null ) ]
    pcs  =  pcs[ show_sexes ]
    
    p  =  pcs %>%
        lapply( function(x) x$log10p ) %>%
        as_tibble %>%
        bind_cols( pcs[[ 1 ]][ 'VARIABLE' ], . ) %>%
        left_join( tissue_descriptions ) %>%
        mutate( tissue_group = factor( tissue_group, levels = unique( tissue_group ) ) )
    
    if (ordered) {
        p$VARIABLE  =  factor( p$VARIABLE,
                               ordered = TRUE,
                               levels  = p$VARIABLE )
    }
    
    # p  =  p %>%
    #     dplyr::select( c( 'VARIABLE', !!show_sexes ) )
    # 
    # if (!ordered) {
    #     p$VARIABLE  =  factor( p$VARIABLE,
    #                            ordered = TRUE,
    #                            levels  = p$VARIABLE )
    # }
    
    if (is.null( threshold )) {
        threshold = -log10( .05/nrow( p ) )
    }
    
    significance_line  =  list(
        type = 'line',
        line = list( color = 'grey' ),
        xref = 'paper',
        yref = 'y',
        x0   = 0,
        x1   = 1,
        y0   = threshold,
        y1   = threshold
    )
    
    showlegend  =  TRUE
    
    margin = 0.005
    widths = p$tissue_group %>%
        table %>%
        '/'(sum(.)) %>%
        unname %>%
        '-'( c( margin/2,
                rep(margin, length(unique(p$tissue_group))-2),
                margin/2 ) )
    
    p %>%
        split( p$tissue_group ) %>%
        map2( names( . ), function( my_data, my_name ){
            plot_ly( data = my_data, x = ~droplevels( VARIABLE ), showlegend = showlegend ) %>%
                {
                    showlegend  <<-  FALSE
                    if ('female' %in% show_sexes)
                        add_trace( ., y = ~female, type = 'bar',
                                   name = 'Female',
                                   marker = list( color = ifelse( single,
                                                                  '#424242',
                                                                  '#E41A1C' ) ),
                                   width = ifelse( single,
                                                   1,
                                                   .5 ),
                                   opacity = ifelse( single,
                                                     1,
                                                     .5 ) )
                    else
                        .
                } %>%
                {
                    if ('both_sexes' %in% show_sexes)
                        add_trace( ., y = ~both_sexes, type = 'bar',
                                   name = 'Both sexes',
                                   marker = list( color = '#424242' ),
                                   width = ifelse( single,
                                                   1,
                                                   .5 ) )
                    else
                        .
                } %>%
                # Hasn't been done in male
                # {
                #     if ('male' %in% show_sexes)
                #         add_trace( ., y = ~male, type = 'bar',
                #                    name = 'Male',
                #                    marker = list( color = '#377EB8' ),
                #                    width = .2, opacity = .5 )
                #     else
                #         .
                # } %>%
                plotly::layout( annotations = list(
                    text = my_name,
                    xref = "paper",
                    yref = "paper",
                    yanchor = "bottom",
                    xanchor = "center",
                    align = "center",
                    x = 0.5,
                    y = 1,
                    showarrow = FALSE),
                    xaxis = list( tickangle = 45 ) )
        } ) %>%
        # '['(3:10) %>%
        subplot(
            margin = margin, shareY = TRUE,
            widths = widths
            ) %>%
        plotly::layout( yaxis  = list( title = '-log 10 p-value' ),
                        xaxis  = list( title = '' ),
                        shapes = significance_line )
}


scatterplot_fuma  =  function( phenotype_x,
                               phenotype_y,
                               dataset,
                               threshold_x = NULL,
                               threshold_y = NULL,
                               data_path = FUMA_DATA_PATH,
                               sex_x = 'both_sexes',
                               sex_y = 'both_sexes',
                               sex_pca_x = 'female',
                               sex_pca_y = 'female',
                               x_label = ifelse( is.numeric(phenotype_x),
                                                 paste0( 'PC', phenotype_x ),
                                                 phenotype_x ),
                               y_label = ifelse( is.numeric(phenotype_y),
                                                 paste0( 'PC', phenotype_y ),
                                                 phenotype_y ),
                               x_results_filename = sprintf( '%s/body_%s%s%s_%s/magma_exp_%s.gsa.out',
                                                             data_path,
                                                             ifelse( sex_x == sex_pca_x | !is.numeric( phenotype_x ),
                                                                     '',
                                                                     str_sub(sex_pca_x, 1, 1) ),
                                                             ifelse( is.numeric( phenotype_x ),
                                                                     'pc',
                                                                     '' ),
                                                             phenotype_x,
                                                             sex_x,
                                                             dataset ),
                               y_results_filename = sprintf( '%s/body_%s%s%s_%s/magma_exp_%s.gsa.out',
                                                             data_path,
                                                             ifelse( sex_y == sex_pca_y | !is.numeric( phenotype_y ),
                                                                     '',
                                                                     str_sub(sex_pca_y, 1, 1) ),
                                                             ifelse( is.numeric( phenotype_y ),
                                                                     'pc',
                                                                     '' ),
                                                             phenotype_y,
                                                             sex_y,
                                                             dataset ),
                               tissue_descriptions = sprintf( '%s/%s.rds',
                                                              data_path,
                                                              dataset ) %>%
                                   read_rds,
                               title_text = sprintf( 'Comparison of tissue enrichment for %s and %s',
                                                     x_label,
                                                     y_label ),
                               show_subcategories = FALSE,
                               separate_xy = FALSE,
                               showlegend = FALSE,
                               export_image = FALSE,
                               flip_line = FALSE,
                               annot = list(),
                               add_regression_line = FALSE ){
    
    data_to_plot = c( x_results_filename, y_results_filename ) %>%
        lapply( function( filename ){
            read_delim( filename,
                        delim = ' ',
                        comment = '#',
                        trim_ws = TRUE )
        } )
    threshold_x  =  ifelse( is.null( threshold_x ),
                            .05/nrow( data_to_plot[[1]] ),
                            threshold_x )
    threshold_y  =  ifelse( is.null( threshold_y ),
                            .05/nrow( data_to_plot[[2]] ),
                            threshold_y )
    combine_by  =  data_to_plot %>%
        lapply( colnames ) %>%
        c( list( c( 'VARIABLE', 'FULL_NAME' ) ) ) %>%
        reduce( intersect )
    
    data_to_plot  =  data_to_plot %>%
        c( list( by = combine_by ) ) %>%
        do.call( full_join, args = . ) %>%
        left_join( tissue_descriptions ) %>%
        transmute( phenotype = VARIABLE,
                   full_name = FULL_NAME,
                   tissue_group = tissue_group,
                   xb = BETA.x,
                   xs = SE.x,
                   xp = P.x,
                   yb = BETA.y,
                   ys = SE.y,
                   yp = P.y,
                   x_significant = xp < threshold_x,
                   y_significant = yp < threshold_y,
                   separate_p    = 2*pnorm( -abs( xb - yb ) / sqrt( xs^2 + ys^2 ) ),
                   separate_sig  = (separate_p < .05/sum( x_significant | y_significant ) &
                                        ( x_significant | y_significant )),
                   separate_neg_p = 2*pnorm( -abs( xb + yb ) / sqrt( xs^2 + ys^2 ) ),
                   separate_neg_sig  = (separate_neg_p < .05/sum( x_significant | y_significant ) &
                                            ( x_significant | y_significant )),
                   group         = c( 'Not significant',
                                      x_label,
                                      y_label,
                                      'Both' )[x_significant + 2*y_significant + 1] ) %>%
        mutate_at( c( 'xp', 'yp', 'separate_p', 'separate_neg_p' ),
                   formatC, format = 'e', digits = 2 ) %>%
        # mutate_at( c( 'xb', 'yb' ),
        #            formatC, digits = 2 ) %>%
        mutate( description = sprintf( '%s
%s b = %s, p = %s
%s b = %s, p = %s
p_diff = %s',
                                       full_name,
                                       x_label, formatC( xb, digits = 2 ), xp,
                                       y_label, formatC( yb, digits = 2 ), yp,
                                       ifelse( flip_line,
                                               separate_neg_p,
                                               separate_p ) ) ) %>%
        filter( !is.na( xb ) & !is.na( yb ) )
    
    if (show_subcategories) {
        data_to_plot  =  data_to_plot %>%
            mutate( group = tissue_group %>% str_to_sentence %>% as.factor )
        n_groups  =  data_to_plot$group %>% nlevels
        group_palette = brewer_pal( 'qual','Paired' )( n_groups )
        
        names( group_palette )  =  levels( data_to_plot$group )
        group_styles  =  data_to_plot$group %>%
            levels %>%
            lapply( function( x ) list(
                target = x,
                name   = x,
                value  =  list( marker = list( color = group_palette[ x ] %>% unname ) )
            ) )
        transforms_list  =  list(
            list(
                type = 'groupby',
                groups = data_to_plot$group,
                styles = group_styles
            )
        )
        data_to_plot  =  data_to_plot %>%
            mutate( group = paste0( group, c( '', ' (x)', ' (y)', ' (b)' )[ x_significant + 2*y_significant + 1 ] ) )
        transforms_list  =  transforms_list %>%
            add_a_transform( data_to_plot$group,
                             c( ' (x)', ' (y)', ' (b)' ) %>% setNames( ., . ),
                             list( ' (x)' = list( symbol = 'triangle-down',
                                                  size   = 10 ),
                                   ' (y)' = list( symbol = 'triangle-left',
                                                  size   = 10 ),
                                   ' (b)' = list( symbol = 'circle',
                                                  size   = 10 ) ),
                             base_marker_additions = list( opacity =.42 ) )
    } else {
        transforms_list  =  list(
            list(
                type = 'groupby',
                groups = data_to_plot$group,
                styles = list(
                    list( target = 'Not significant',
                          name   = 'Not significant',
                          value = list( marker = list( color  = '#848484',
                                                       symbol = 'circle-open' ) ) ),
                    list( target = x_label,
                          name   = x_label,
                          value = list( marker = list( color  = 'blue',
                                                       symbol = 'triangle-down' ) ) ),
                    list( target = y_label,
                          name   = y_label,
                          value = list( marker = list( color  = 'red',
                                                       symbol = 'triangle-left' ) ) ),
                    list( target = 'Both',
                          name   = 'Both',
                          value = list( marker = list( color  = 'black',
                                                       symbol = 'circle' ) ) )
                )
            )
        )
    }
    
    if (separate_xy) {
        data_to_plot  =  data_to_plot %>%
            mutate( group = paste0( group, c( '', ' (s)' )[ separate_sig+1 ] ) )
        transforms_list  =  add_a_transform( transforms_list,
                                             data_to_plot$group,
                                             ' (s)',
                                             new_marker_additions = list( opacity = 1,
                                                                          size    = 10 ),
                                             base_marker_additions = list( opacity = .42,
                                                                           size    = 5 ) )
    }
    xrange = range( data_to_plot$xb )
    xsize  = xrange[ 2 ] - xrange[ 1 ]
    xrange = xrange + c( -xsize, xsize )/20
    xrange = c( min( xrange[ 1 ], 0 ),
                max( xrange[ 2 ], 0 ) )
    yrange = range( data_to_plot$yb )
    ysize  = yrange[ 2 ] - yrange[ 1 ]
    yrange = yrange + c( -ysize, ysize )/20
    yrange = c( min( yrange[ 1 ], 0 ),
                max( yrange[ 2 ], 0 ) )
    identity  =  line_coordinates( data_to_plot$xb,
                                   data_to_plot$yb,
                                   slope = ifelse( flip_line,
                                                   -1,
                                                   1 ) )
    
    all_lines  =  list( list( type = 'line',
                              line = list( color = 'lightgray',
                                           dash  = 'dash' ),
                              x0 = identity$x[ 1 ],
                              x1 = identity$x[ 2 ],
                              y0 = identity$y[ 1 ],
                              y1 = identity$y[ 2 ] ) )
    
    # if (add_regression_line) {
    #     regression  =  lm( yb~xb-1,
    #                        data = data_to_plot )
    #     regression_line  =  line_coordinates( data_to_plot$xb,
    #                                           data_to_plot$yb,
    #                                           slope = regression$coef %>%
    #                                               unname )
    #     all_lines  =  all_lines %>%
    #         c( list( list( type = 'line',
    #                        line = list( color = 'black',
    #                                     dash  = 'solid' ),
    #                        x0 = regression_line$x[ 1 ],
    #                        x1 = regression_line$x[ 2 ],
    #                        y0 = regression_line$y[ 1 ],
    #                        y1 = regression_line$y[ 2 ]) ) )
    # }
    
    if (export_image) {
        width = 1000 + 150*showlegend
        height = 800
        main_title = NULL
        x_label = ''
        y_label = ''
        if (showlegend) {
            m  =  list( l   = 10,
                        r   = 250,
                        b   = 10,
                        t   = 10,
                        pad = 0 )
        }
    } else {
        width = height = NULL
    }
    
    # p0  =  max( xrange[ 1 ], yrange[ 1 ] )
    # p1  =  min( xrange[ 2 ], yrange[ 2 ] )
    # flip_p0  =  max( xrange[ 1 ], -yrange[ 2 ] )
    # flip_p1  =  min( xrange[ 2 ], -yrange[ 1 ] )
    
    # Accepts RGB or RGBA
    error_bar_color  =  '#00000011'
    
    m = list()
    
    # if (!is.null( selected_group )) {
    #     data_to_plot  =  data_to_plot %>%
    #         filter( subcategory %in% selected_group )
    # }
    
    p  =  plot_ly( data_to_plot, type = 'scatter', mode = 'none',
                   width = width,
                   height = height ) %>%
        add_markers( x = ~xb, y = ~yb, type = 'scatter', mode = 'markers',
                     name = ' ',
                     error_x = ~list( array = xs,
                                      color = error_bar_color ),
                     error_y = ~list( array = ys,
                                      color = error_bar_color ),
                     transforms = transforms_list,
                     hoverinfo = 'text',
                     text = ~description ) %>%
        plotly::layout( title = title_text,
                        xaxis  = list( title = x_label,
                                       range = xrange ),
                        yaxis  = list( title = y_label,
                                       range = yrange ),
                        # shapes = all_lines,
                        margin = m,
                        showlegend = showlegend,
                        annotations = annot )
    
    if (add_regression_line) {
        p  =  p %>%
            add_deming_regression_line( data_to_plot$xb,
                                        data_to_plot$yb,
                                        data_to_plot$xs,
                                        data_to_plot$ys,
                                        all_lines )
    }
    p
}


get_pascal_results  =  function( phenotype = NULL,
                                 sex = 'both_sexes',
                                 sex_pca = 'female',
                                 results_filename = sprintf( PASCAL_FILENAME,
                                                             PASCAL_PATH,
                                                             ifelse( sex == sex_pca | !is.numeric( phenotype ),
                                                                     '',
                                                                     str_sub( sex_pca, 1, 1 ) ),
                                                             ifelse( is.numeric( phenotype ),
                                                                     'pc',
                                                                     '' ),
                                                             phenotype,
                                                             sex ),
                                 description_file = DEPICT_DESCRIPTION,
                                 mp_only = TRUE,
                                 fdr = TRUE,
                                 threshold = 0.05,
                                 bonferroni = !fdr,
                                 ... ){
    
    results = read_tsv( results_filename,
                        col_types = cols( Name = 'c',
                                          chi2Pvalue = 'n',
                                          empPvalue = 'n') )
    if (mp_only) {
        results = results %>%
            filter( str_detect( Name, '^MP:[0-9]+$' ) )
    }
    if (fdr) {
        results  =  results %>%
            mutate_if( is.numeric, p.adjust, 'fdr' ) %>%
            rename_all( str_replace, 'Pvalue', 'fdr' ) %>%
            full_join( results ) %>%
            filter( empfdr <= threshold | chi2fdr <= threshold )
    } else if (bonferroni) {
        threshold  =  threshold / nrow( results )
        results  =  results %>%
            filter_if( is.numeric, any_vars( . <= threshold ) )
    }
    
    read_tsv( description_file,
              col_types = cols_only( ID = 'c',
                                     Desc = 'c' ) ) %>%
        rename( Name = ID ) %>%
        right_join( results )
}


create_pascal_tree  =  function( phenotype,
                                 sex = 'both_sexes',
                                 sex_pca = '',
                                 prune = TRUE,
                                 description_file = DEPICT_DESCRIPTION,
                                 title = sprintf( 'Pathway enrichment for %s',
                                                  ifelse( is.numeric(phenotype),
                                                          paste0( 'PC', phenotype ),
                                                          phenotype ) ),
                                 ... ){
    results  =  get_pascal_results( phenotype,
                                    sex = sex,
                                    sex_pca = sex_pca,
                                    description_file = description_file,
                                    ... )
    descriptions  =  read_tsv( description_file,
                               col_types  =  cols_only( ID = 'c',
                                                        Desc = 'c',
                                                        MetaCluster = 'c',
                                                        MetaCluster_Desc = 'c' ) ) %>%
        rename( Name = ID )
    if (prune) {
        new_descriptions  =  descriptions %>%
            filter( Name %in% results$Name )
        while (any( !(new_descriptions$MetaCluster %in% new_descriptions$Name) )) {
            new_descriptions  =  descriptions %>%
                filter( Name %in% new_descriptions$MetaCluster ) %>%
                full_join( new_descriptions )
        }
    }
    
    my_graph = mapply( c,
                       new_descriptions$MetaCluster,
                       new_descriptions$Name,
                       SIMPLIFY = FALSE ) %>%
        unlist %>%
        unname %>%
        graph
    
    vertex_descriptions  =  my_graph %>%
        V %>%
        '$'( 'name' ) %>%
        tibble( Name = . ) %>%
        left_join( new_descriptions )
    
    vertex_descriptions  =  results %>%
        mutate( lowerp = apply( results[ , c( 'empPvalue', 'chi2Pvalue' ) ], 1, min ) ) %>%
        left_join( vertex_descriptions, . )
    
    V( my_graph )$label  =  vertex_descriptions$Desc
    
    graph_layout = layout.fruchterman.reingold( my_graph )
    
    pathways = V( my_graph )
    metacluster_edges = get.edgelist( my_graph ) %>%
        as.data.frame( stringsAsFactors = FALSE ) %>%
        mutate_all( factor, levels = pathways$name )
    
    n_pathways = length( pathways )
    n_edges = length( metacluster_edges[1]$V1 )
    
    Xn  =  graph_layout[ , 1 ]
    Yn  =  graph_layout[ , 2 ]
    
    pathways_network  =  plot_ly() %>%
        add_markers( x = ~Xn[ !is.na( vertex_descriptions$lowerp ) ],
                     y = ~Yn[ !is.na( vertex_descriptions$lowerp ) ],
                     mode = "markers",
                     marker = list( 
                         color = -log10(vertex_descriptions$lowerp[ !is.na( vertex_descriptions$lowerp ) ]),
                         colorscale = 'Viridis',
                         colorbar = list(
                             title = '-log10 p-value'
                         )
                     ),
                     text = pathways$label[ !is.na( vertex_descriptions$lowerp ) ],
                     hoverinfo = "text" ) %>%
        add_markers( x = ~Xn[ is.na( vertex_descriptions$lowerp ) ],
                     y = ~Yn[ is.na( vertex_descriptions$lowerp ) ],
                     mode = "markers",
                     marker = list(
                         color = 'lightgrey'
                     ),
                     text = pathways$label[ is.na( vertex_descriptions$lowerp ) ],
                     hoverinfo = "text" )
    
    edge_shapes = list()
    for(i in 1:n_edges) {
        v0 = metacluster_edges[ i, ]$V1
        v1 = metacluster_edges[ i, ]$V2
        
        edge_shape = list(
            type = "line",
            line = list( color = "#030303", width = 0.3 ),
            x0 = Xn[ v0 ],
            y0 = Yn[ v0 ],
            x1 = Xn[ v1 ],
            y1 = Yn[ v1 ]
        )
        
        edge_shapes[[ i ]] = edge_shape
    }
    
    axis  = list( title = "",
                  showgrid = FALSE,
                  showticklabels = FALSE,
                  zeroline = FALSE )
    
    p =  layout(
        pathways_network,
        title = title,
        shapes = edge_shapes,
        xaxis = axis,
        yaxis = axis
    )
    p
}


pascal_heatmap  =  function( phenotypes = NULL,
                             sex = 'both_sexes',
                             sex_pca = 'female',
                             prune = TRUE,
                             pascal_path = PASCAL_PATH,
                             pascal_pattern = PASCAL_FILENAME %>%
                                 str_match( '[^/]+$' ) %>%
                                 sprintf( ifelse( is.null( phenotypes ),
                                                  sprintf( '(?:%s(?:pc)?[1-4]|2100[12]_irnt|whr)',
                                                           ifelse( sex == sex_pca,
                                                                   '',
                                                                   str_sub( sex_pca, 1, 1 ) ) ),
                                                  paste( phenotypes, collapse = '|' ) %>%
                                                      paste0( '(?:', ., ')' ) ),
                                          '',
                                          '',
                                          sex ),
                             descriptions = read_rds( PROCESSED_DESCRIPTION ),
                             mp_only = TRUE,
                             fdr = TRUE,
                             threshold = 0.05,
                             bonferroni = !fdr,
                             column_order = c( 'fpc1',
                                               '21002_irnt',
                                               'fpc2',
                                               '21001_irnt',
                                               'whr',
                                               'fpc3',
                                               'fpc4' ),
                             column_names = COLUMN_NAMES,
                             # title = sprintf( 'Pathway enrichment for %s',
                             #                  ifelse( is.numeric(phenotype),
                             #                          paste0( 'PC', phenotype ),
                             #                          phenotype ) ),
                             ... ){
    pascal_results  =  list.files( pascal_path,
                                   pattern = pascal_pattern,
                                   full.names = TRUE ) %>%
        lapply( function( filename ){
            phenotype  =  filename %>%
                str_match( '/body_(.+)_both_sexes' ) %>%
                '['( 1, 2 )
            results  =  read_tsv( filename,
                                  col_types = cols( Name = 'c',
                                                    chi2Pvalue = 'n',
                                                    empPvalue = 'n') )
            if (mp_only) {
                results  =  results %>%
                    filter( str_detect( Name, '^MP:[0-9]+$' ) )
            }
            if (fdr) {
                results  =  results %>%
                    mutate_if( is.numeric, p.adjust, 'fdr' )
            }
            results %>%
                mutate_if( is.numeric, function(x) -log10(x) ) %>%
                transmute( Name = Name,
                           !!phenotype := apply( .[ , c( 'empPvalue', 'chi2Pvalue' ) ], 1, min ) )
        } ) %>%
        reduce( full_join ) %>%
        right_join( descriptions, . )
    
    if (bonferroni) {
        threshold  =  -log10( threshold / nrow( pascal_results ) )
    } else {
        threshold  =  -log10( threshold )
    }
    # pascal_results  =  pascal_results %>%
    #     filter( apply( pascal_results[ , -c(1:ncol(descriptions)) ], 1, function(x) any(x >= threshold) ) )
    pascal_results  =  pascal_results[ , -c(1:ncol(descriptions)) ] %>%
        apply( 1, function(x) any(x >= threshold) ) %>%
        filter( pascal_results, . )
    
    if (is.null( column_order )) {
        pascal_results  =  (hclust( pascal_results[ , -c(1:ncol(descriptions)) ] %>% t %>% dist )$order + ncol(descriptions)) %>%
            c( 1:ncol(descriptions), . ) %>%
            '['( pascal_results, , . )
    } else {
        pascal_results  =  pascal_results[ , c(colnames( descriptions ), column_order) ]
    }
    
    
    scale_colors  =  RColorBrewer::brewer.pal( length(unique(pascal_results$group)),'Dark2') %>%
        setNames( nm = pascal_results$group %>% unique ) %>%
        lapply( function(group_color){
            values  =  unique( scales::rescale( unlist( pascal_results[ , -c(1:ncol(descriptions)) ] %>% c( 0, . ) ) ) )
            o  =  order( values, decreasing = FALSE )
            cols  =  scales::col_numeric( colorRamp(c("#FFFFFF", group_color), interpolate="spline"), domain = NULL )( values )
            setNames( data.frame( values[ o ], cols[ o ] ), NULL )
        } )
    
    pascal_results %>%
        split( .$group ) %>%
        map2( names( . ), function( mydata, myname ) {
            mydata = hclust( mydata[ , -c(1:ncol(descriptions)) ] %>%
                                 (function(x){
                                     dist( x ) + dist( x < threshold )
                                 }) %>%
                                 dist )$order %>%
                '['( mydata, .,  )
            mydata %>%
                # mutate_if( is.numeric, function(x) replace(x, x < threshold, NA) ) %>%
                plot_ly( data = ., x = colnames( pascal_results )[ -c(1:ncol(descriptions)) ] ) %>%
                add_heatmap( name = myname,
                             y = ~Desc,
                             z = mydata[, -c(1:ncol(descriptions))] %>%
                                 as.matrix %>%
                                 replace( . < threshold, NA ),
                             text = mydata[, -c(1:ncol(descriptions))] %>%
                                 map2( names(.),
                                       function( values, name ) sprintf('%s\n%s\nfdr = %g',
                                                                        name,
                                                                        mydata$Desc,
                                                                        10^(-values))) %>%
                                 as_tibble %>%
                                 as.matrix,
                             hoverinfo = 'text',
                             colorscale = scale_colors[[ myname ]]
                ) %>%
                layout( yaxis = list( #title = myname,
                    tickmode = 'array',
                    tickvals = c(),
                    legend = list( title = list( text = myname ) ) ),
                    xaxis = list( 
                        tickmode = 'array',
                        tickvals = column_order,
                        ticktext = column_names ) ) %>%
                colorbar( title = myname )
        } ) %>%
        subplot( nrows = length(unique(pascal_results$group)),
                 margin = .001, shareX = TRUE,
                 heights = (table( pascal_results$group )%>%unname)/nrow(pascal_results) )
    
}


fuma_heatmap  =  function( phenotypes = c( 'fpc1',
                                           '21002_irnt',
                                           'fpc2',
                                           '21001_irnt',
                                           'whr',
                                           'fpc3',
                                           'fpc4' ),
                           dataset = 'gtex_v8_ts_avg_log2TPM',
                           data_path = FUMA_DATA_PATH %>%
                               sprintf( PREFIX ),
                           sex = 'both_sexes',
                           sex_pca = 'female',
                           fuma_filenames = sprintf( '%s/body_%s_%s/magma_exp_%s.gsa.out',
                                                     data_path,
                                                     phenotypes,
                                                     sex,
                                                     dataset ),
                           descriptions = sprintf( '%s/%s.rds',
                                                   data_path,
                                                   dataset ) %>%
                               read_rds %>%
                               select( tissue = VARIABLE,
                                       tissue_group ),
                           title_text = 'MAGMA tissue enrichment',
                           ordered = TRUE,
                           prune = TRUE,
                           fdr = FALSE,
                           threshold = 0.05,
                           bonferroni = !fdr,
                           column_order = c( 'fpc1',
                                             '21002_irnt',
                                             'fpc2',
                                             '21001_irnt',
                                             'whr',
                                             'fpc3',
                                             'fpc4' ),
                           column_names = COLUMN_NAMES,
                           group_order = c( 'Brain',
                                            'Hormonal',
                                            'Adipose',
                                            'Cells',
                                            'Digestive',
                                            'Female',
                                            'Other',
                                            'Vascular',
                                            'Male',
                                            'Skin',
                                            'Prenatal',
                                            'Infancy',
                                            'Childhood',
                                            'Adulthood' ),
                           ... ){
    fuma_results  =  fuma_filenames %>%
        setNames( phenotypes ) %>%
        lapply( read_delim,
                delim = ' ',
                comment = '#',
                trim_ws = TRUE )
    fuma_results  =  phenotypes %>%
        lapply( function( phenotype ){
            fuma_results[[ phenotype ]] %>%
                transmute( tissue = VARIABLE,
                           !!phenotype := -log10( P ) )
        } ) %>%
        reduce( full_join ) %>%
        right_join( descriptions, . ) %>%
        mutate( tissue = str_replace_all( tissue, '_', ' ' ) )
    descriptions  =  descriptions %>%
        mutate( tissue = str_replace_all( tissue, '_', ' ' ) )
    
    if (fdr) {
        fuma_results  =  fuma_results %>%
            mutate_if( is.numeric, function( x ){
                -(p.adjust( 10^(-x), 'fdr' ) %>%
                    log10)
            } )
    } else if (bonferroni) {
        threshold  =  -log10( threshold / nrow( fuma_results ) )
    } else {
        threshold  =  -log10( threshold )
    }
    
    if (prune) {
        to_keep  =  fuma_results$tissue_group %>%
            unique %>%
            sapply( function( tissue_group_name ){
                any( fuma_results[ fuma_results$tissue_group == tissue_group_name,
                                   -c(1:ncol(descriptions)) ] > threshold )
            } )
        to_keep  =  names(to_keep)[ to_keep ]
        fuma_results  =  fuma_results %>%
            filter( tissue_group %in% to_keep )
        descriptions  =  descriptions %>%
            filter( tissue_group %in% to_keep )
        # if (!is.null( group_order )) {
        #     group_order = group_order[ group_order %in% to_keep ]
        # }
    }
    
    group_order  =  group_order[ group_order %in% fuma_results$tissue_group ]
    
    if (is.null( column_order )) {
        fuma_results  =  (hclust( fuma_results[ , -c(1:ncol(descriptions)) ] %>% t %>% dist )$order + ncol(descriptions)) %>%
            c( 1:ncol(descriptions), . ) %>%
            '['( fuma_results, , . )
    } else {
        fuma_results  =  fuma_results[ , c(colnames( descriptions ), column_order) ]
    }
    
    
    scale_colors  =  RColorBrewer::brewer.pal( length(unique(fuma_results$tissue_group)),'Paired') %>%
        setNames( nm = fuma_results$tissue_group %>% unique ) %>%
        map2( names( . ), function( group_color, group_name ){
            values  =  unique( scales::rescale( unlist( fuma_results[ , -c(1:ncol(descriptions)) ] %>% c( 0, . ) ) ) )
            o  =  order( values, decreasing = FALSE )
            # if (all( fuma_results[ fuma_results$tissue_group == group_name,
            #                        c(1:ncol(descriptions)) ] < threshold )) {
            #     return( setNames( data.frame( values[ o ], '#FFFFFF' ), NULL ) )
            # } else {
            cols  =  scales::col_numeric( colorRamp(c("#FFFFFF", group_color), interpolate="spline"), domain = NULL )( values )
            setNames( data.frame( values[ o ], cols[ o ] ), NULL )
            # }
        } )
        # lapply( function(group_color){
        #     values  =  unique( scales::rescale( unlist( fuma_results[ , -c(1:ncol(descriptions)) ] %>% c( 0, . ) ) ) )
        #     o  =  order( values, decreasing = FALSE )
        #     cols  =  scales::col_numeric( colorRamp(c("#FFFFFF", group_color), interpolate="spline"), domain = NULL )( values )
        #     setNames( data.frame( values[ o ], cols[ o ] ), NULL )
        # } )
    
    fuma_results %>%
        split( .$tissue_group ) %>%
        '['( group_order ) %>%
        map2( names( . ), function( mydata, myname ) {
            if (nrow( mydata ) > 1){
                if (str_detect( 'gtex', dataset )) {
                    mydata = hclust( mydata[ , -c(1:ncol(descriptions)) ] %>%
                                         (function(x){
                                             dist( x ) + dist( x < threshold )
                                         }) %>%
                                         dist )$order %>%
                        '['( mydata, .,  )
                } else {
                    mydata  =  mydata[ nrow(mydata):1, ]
                }
            }
            
            mydata %>%
                plot_ly( data = ., x = colnames( fuma_results )[ -c(1:ncol(descriptions)) ] ) %>%
                add_heatmap( name = myname,
                             y = ~tissue,
                             z = mydata[, -c(1:ncol(descriptions))] %>%
                                 as.matrix %>%
                                 replace( . < threshold, 0 ),
                             text = mydata[, -c(1:ncol(descriptions))] %>%
                                 map2( names(.),
                                       function( values, name ) sprintf('%s\n%s\n%s = %g',
                                                                        name,
                                                                        mydata$tissue,
                                                                        ifelse( fdr,
                                                                                'fdr',
                                                                                'p-value' ),
                                                                        10^(-values))) %>%
                                 as_tibble %>%
                                 as.matrix,
                             hoverinfo = 'text',
                             colorscale = scale_colors[[ myname ]]
                ) %>%
                layout( yaxis = list( legend = list( title = list( text = myname ) ) ),
                        xaxis = list( 
                            tickmode = 'array',
                            tickvals = column_order,
                            ticktext = column_names ) ) %>%
                colorbar( title = myname )
        } ) %>%
        subplot( nrows = length(unique(fuma_results$tissue_group)),
                 margin = .001, shareX = TRUE,
                 heights = (table( fuma_results$tissue_group )[ group_order ]%>%unname)/nrow(fuma_results) )
    
}









