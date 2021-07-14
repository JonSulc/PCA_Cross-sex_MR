require( tidyverse )

MAX_REGEX      =  200
CHUNK_SIZE     =  10
MAX_DOWNLOADS  =  500
OUTPUT_PATH    =  'data/neale_files'

download_neale_files  =  function( phenotype_ids,
                                   irnt = TRUE,
                                   max_regex = MAX_REGEX,
                                   chunk_size = CHUNK_SIZE,
                                   max_downloads = MAX_DOWNLOADS,
                                   output_path = OUTPUT_PATH,
                                   reference_file = list.files( 
                                           path = output_path,
                                           pattern = 'UKBB.*Manifest.*[.]tsv',
                                           full.names = TRUE
                                       )[ 1 ],
                                   sex = c( 'both_sexes', 'female', 'male' ),
                                   category = 'unsorted' ){
    
    if (length( phenotype_ids ) > max_regex) {
        n_downloads  =  download_neale_files( phenotype_ids[ 1:max_regex ],
                                              irnt = irnt,
                                              max_regex = max_regex,
                                              chunk_size = chunk_size,
                                              max_downloads = max_downloads,
                                              output_path = output_path,
                                              reference_file = reference_file,
                                              sex = sex,
                                              category = category )
        
        download_neale_files( phenotype_ids[ (max_regex+1):length( phenotype_ids ) ],
                              irnt = irnt,
                              max_regex = max_regex,
                              chunk_size = chunk_size,
                              max_downloads = max_downloads - n_downloads,
                              output_path = output_path,
                              reference_file = reference_file,
                              sex = sex,
                              category = category )
    }
    if (length( sex ) > 1)
        sex = paste( sex, collapse = '|') %>%
            paste0( '(', ., ')' )
    
    phenotype_ids  =  paste0( '^', phenotype_ids, ifelse( irnt, '(_irnt|)$', '(_raw|)$' ) ) %>%
        paste( collapse = '|' )
    
    if (!file.exists( output_path )) {
      dir.create( output_path )
    }
    
    existing_files  =  list.files( output_path,
                                   recursive = TRUE,
                                   pattern = '[.]gz' ) %>%
        str_match( '[^/]+$' ) %>%
        c
    
    to_download  =  read_tsv( reference_file ) %>%
        filter( str_detect( .$'Phenotype Code', phenotype_ids ) ) %>%
        select( wget = 'wget command' ) %>%
        mutate( wget = str_replace( wget, '[.]bgz$', '.gz' ) ) %>%
        filter( str_detect( .$wget, paste0( '-O .*[.]', sex, '([.]v2)?[.]tsv[.]gz$' ) ) ) %>%
        transmute( url = str_match( wget, '^wget (http.*) -O' )[ , 2 ],
                   folder = paste( output_path,
                                   str_match( wget, '-O .*[.](both_sexes|female|male)([.]v2)?[.]tsv[.]gz$' )[ , 2 ],
                                   category,
                                   sep = '/' ),
                   filename = str_match( wget, '-O (.*[.]gz$)' )[ , 2 ] ) %>%
        filter( !( filename %in% existing_files ) ) %>%
        mutate( filename = paste0( folder,
                                   '/',
                                   filename ) )
    
    to_download$folder %>%
        unique %>%
        walk( dir.create, showWarnings = FALSE, recursive = TRUE )
    to_download <<- to_download
    if (nrow( to_download ) != 0) {
        for (chunk in 1:ceiling( nrow( to_download ) / chunk_size )) {
            if (chunk * chunk_size > max_downloads)
                stop( 'Maximum number of downloads reached.' )
            download.file( to_download[ (( chunk - 1 ) * chunk_size + 1 ):min( nrow( to_download ), (chunk * chunk_size) ), ]$url,
                           to_download[ (( chunk - 1 ) * chunk_size + 1 ):min( nrow( to_download ), (chunk * chunk_size) ), ]$filename,
                           method = 'libcurl')
        }
    }
    nrow( to_download )
}










