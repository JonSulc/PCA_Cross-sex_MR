require( meta )
require( tidyverse )

meta_analyze  =  function( mr_results_x,
                           mr_results_y,
                           name_x = '',
                           name_y = '',
                           ... ){
    common_rows  =  list( mr_results_x$b, mr_results_y$b ) %>%
        lapply( rownames ) %>%
        do.call( intersect, . )
    common_columns  =  list( mr_results_x$b, mr_results_y$b ) %>%
        lapply( colnames ) %>%
        do.call( intersect, . )
    mr_results_x  =  mr_results_x %>%
        lapply( '[', common_rows, common_columns, drop = FALSE )
    mr_results_y  =  mr_results_y %>%
        lapply( '[', common_rows, common_columns, drop = FALSE )
    
    mr_results  =  mr_results_x[ c( 'b', 'se', 'pval' ) ]
    for (exposure in common_rows) {
        for (outcome in common_columns) {
            meta_results  =  metagen( c( mr_results_x$b[  exposure, outcome ],
                                         mr_results_y$b[  exposure, outcome ] ),
                                      c( mr_results_x$se[ exposure, outcome ],
                                         mr_results_y$se[ exposure, outcome ] ),
                                      ... )
            mr_results$b[    exposure, outcome ]  =  meta_results$TE.fixed
            mr_results$se[   exposure, outcome ]  =  meta_results$seTE.fixed
            mr_results$pval[ exposure, outcome ]  =  meta_results$pval.fixed
        }
    }
    mr_results %>%
        c( list( mr_results_x,
                 mr_results_y ) %>%
               setNames( c( name_x, name_y ) ) )
}







