source( 'results_analysis.R' )

all_phenotypes  =  read_rds( './data/all_phenotypes_combined.rds' )

PATH            =  'thr5e-08het0.001rev-1.644854b10_%s'
DATA_PATH       =  'data/'
FUMA_DATA_PATH  =  'fuma/results'

server  =  function( input, output, session ) {
    output$mrplot  =  renderPlotly({
        
        if (input$mrtype == 'mvf') {
            if (input$exposure_type == 'p') {
                exposure = input$pc_number
            } else {
                exposure = input$trait
            }
            
            arguments <<- list( exposure_category_1 = input$exposure_category,
                                exposure_category_2 = input$exposure_category,
                                outcome_category = input$outcome_category,
                                sex_exposure = 'female',
                                sex_ivs = input$sex_ivs,
                                sex_pca = input$sex_pca,
                                sex_outcome = 'male',
                                sexy_exposure = 'male',
                                sexy_ivs = input$sex_ivs,
                                sexy_pca = input$sex_pca,
                                sexy_outcome = 'female',
                                exposure_type_1 = input$exposure_type,
                                exposure_type_2 = input$exposure_type,
                                outcome_type  = input$outcome_type,
                                exp_trait_1 = exposure,
                                exp_trait_2 = exposure,
                                separate_xy   = 'separate_xy' %in% input$plot_options,
                                show_subcategories = 'show_subcategories' %in% input$plot_options,
                                x_path = sprintf( PATH,
                                                  input$mr_method_x ),
                                x_data_path = DATA_PATH,
                                y_path = sprintf( PATH,
                                                  input$mr_method_y ),
                                y_data_path = DATA_PATH,
                                showlegend = 'showlegend' %in% input$plot_options,
                                export_image = 'export_image' %in% input$plot_options,
                                add_regression_line = 'add_regression_line' %in% input$plot_options,
                                x_label = 'Male-specific effect',
                                y_label = 'Female-specific effect' )
            myplot <<- do.call( compare_2_exposures, arguments )
        } else if (input$mrtype == 'bimr') {
            if (input$type_1 == 'p') {
                exposure = input$pc_number_1
            } else {
                exposure = input$trait_1
            }
            
            arguments <<- list( category_1 = input$category_1,
                                category_2 = input$category_2,
                                sex = input$sex,
                                sex_pca = input$sex_pca,
                                type_1 = input$type_1,
                                type_2 = input$type_2,
                                trait = exposure,
                                # separate_xy = 'separate_xy' %in% input$plot_options,
                                show_subcategories = 'show_subcategories' %in% input$plot_options,
                                path = sprintf( PATH,
                                                input$mr_method_x ),
                                data_path = DATA_PATH,
                                export_image = 'export_image' %in% input$plot_options,
                                showlegend = 'showlegend' %in% input$plot_options,
                                # add_regression_line = 'add_regression_line' %in% input$plot_options,
                                flip_line = 'flip_line' %in% input$plot_options )
            do.call( bidirectional_mr_plot, arguments )
        } else if (input$mrtype == '2exp') {
            if (input$type_1 == 'p') {
                exposure_1 = input$pc_number_1
            } else {
                exposure_1 = input$trait_1
            }
            if (input$type_2 == 'p') {
                exposure_2 = input$pc_number_2
            } else {
                exposure_2 = input$trait_2
            }
            
            if (input$sexx_mr == 'male') {
                sexx_exposure  =  'female'
            } else if (input$sexx_mr == 'female') {
                sexx_exposure  =  'male'
            } else {
                sexx_exposure  =  'both_sexes'
            }
            
            if (input$sexy_mr == 'male') {
                sexy_exposure  =  'female'
            } else if (input$sexy_mr == 'female') {
                sexy_exposure  =  'male'
            } else {
                sexy_exposure  =  'both_sexes'
            }
            
            arguments <<- list( exposure_category_1 = input$category_1,
                                exposure_category_2 = input$category_2,
                                outcome_category    = input$outcome_category,
                                sex_exposure = sexx_exposure,
                                sex_pca      = input$sexx_pca,
                                sex_outcome  = input$sexx_mr,
                                sexy_exposure = sexy_exposure,
                                sexy_pca     = input$sexy_pca,
                                sexy_outcome = input$sexy_mr,
                                exposure_type_1 = input$type_1,
                                exposure_type_2 = input$type_2,
                                outcome_type = input$outcome_type,
                                exp_trait_1  = exposure_1,
                                exp_trait_2  = exposure_2,
                                separate_xy  = 'separate_xy' %in% input$plot_options,
                                show_subcategories = 'show_subcategories' %in% input$plot_options,
                                x_path = sprintf( PATH,
                                                  input$mr_method_x ),
                                x_data_path = DATA_PATH,
                                y_path = sprintf( PATH,
                                                  input$mr_method_y ),
                                y_data_path = DATA_PATH,
                                showlegend = 'showlegend' %in% input$plot_options,
                                export_image = 'export_image' %in% input$plot_options,
                                add_regression_line = 'add_regression_line' %in% input$plot_options,
                                flip_line = 'flip_line' %in% input$plot_options )
            myplot <<- do.call( compare_2_exposures, arguments )
        } else if (input$mrtype == 'pcbar') {
            if ('filter_loadings' %in% input$pcbar_options) {
                threshold = NULL
            } else {
                threshold = 0
            }
            arguments  <<-  list( category   = 'body',
                                  pc_number  = input$bar_pc_number,
                                  sex        = 'both_sexes',
                                  sex_pca    = '',
                                  show_sexes = input$show_sexes,
                                  data_path  = DATA_PATH,
                                  threshold  = threshold,
                                  keep_order = !('sort' %in% input$pcbar_options) )
            do.call( barplot_pcs, arguments )
        } else if (input$mrtype == 'fuma_heatmap') {
            arguments <<- list( dataset = input$fuma_tissue_dataset,
                                data_path = FUMA_DATA_PATH,
                                prune = str_detect( input$fuma_tissue_dataset, 'gtex' ) )
            do.call( fuma_heatmap, arguments )
        } else if (input$mrtype == 'pascal_heatmap') {
            pascal_heatmap()
        }
    })
    
    output$type_1_select  =  renderUI({
        if (input$category_1 == 'body') {
            choices = list( 'PC'    = 'p',
                            'Trait' = 't' )
            selected = 'p'
        } else {
            choices = list( 'Trait' = 't' )
            selected = 't'
        }
        radioButtons(
            'type_1',
            'X type',
            choices = choices,
            selected = selected,
            inline = TRUE
        )
    })
    output$trait_select  =  renderUI({
        traits  =  all_phenotypes %>%
            filter( category == input$exposure_category )
        traits  =  setNames( traits$phenotype, traits$description )
        selectInput(
            'trait',
            'MR exposure trait',
            choices = traits
        )
    })
    output$trait_1_select  =  renderUI({
        traits_1  =  all_phenotypes %>%
            filter( category == input$category_1 )
        traits_1  =  setNames( traits_1$phenotype, traits_1$description )
        selectInput(
            'trait_1',
            'X trait',
            choices = traits_1
        )
    })
    
    output$type_2_select  =  renderUI({
        if (input$category_2 == 'body') {
            choices = list( 'PC'    = 'p',
                            'Trait' = 't' )
        } else {
            choices = list( 'Trait' = 't' )
        }
        radioButtons(
            'type_2',
            'X type',
            choices = choices,
            selected = 't',
            inline = TRUE
        )
    })
    output$trait_2_select  =  renderUI({
        traits_2  =  all_phenotypes %>%
            filter( category == input$category_2 )
        traits_2  =  setNames( traits_2$phenotype, traits_2$description )
        selectInput(
            'trait_2',
            'Y trait',
            choices = traits_2
        )
    })
    
    output$outcome_type_select  =  renderUI({
        if (input$outcome_category == 'body') {
            choices = list( 'PC'    = 'p',
                            'Trait' = 't' )
            selected = 'p'
        } else {
            choices = list( 'Trait' = 't' )
            selected = 't'
        }
        radioButtons(
            'outcome_type',
            'Outcome type',
            choices = choices,
            selected = selected,
            inline = TRUE
        )
    })
    output$outcome_select  =  renderUI({
        outcome  =  all_phenotypes %>%
            filter( category == input$outcome_category )
        outcome  =  setNames( outcome$phenotype, outcome$description )
        selectInput(
            'out_trait',
            'Outcome trait',
            choices = outcome
        )
    })
    
    output$out_trait_1_select  =  renderUI({
        out_traits_1  =  all_phenotypes %>%
            filter( category == input$outcome_category_1 )
        out_traits_1  =  setNames( out_traits_1$phenotype, out_traits_1$description )
        selectInput(
            'out_trait_1',
            'X outcome trait',
            choices = out_traits_1
        )
    })
    
    output$out_trait_2_select  =  renderUI({
        out_traits_2  =  all_phenotypes %>%
            filter( category == input$outcome_category_2 )
        out_traits_2  =  setNames( out_traits_2$phenotype, out_traits_2$description )
        selectInput(
            'out_trait_2',
            'Y outcome trait',
            choices = out_traits_2
        )
    })
    output$mr_method_x = renderUI({
        radioButtons( 'mr_method_x',
                      'X MR method',
                      inline = TRUE,
                      choices = list( 'IVW'    = 'mr_ivw',
                                      'Median' = 'mr_weighted_median' ),
                      selected = 'mr_ivw' )
    })
    output$mr_method_y = renderUI({
        radioButtons( 'mr_method_y',
                      'Y MR method',
                      inline = TRUE,
                      choices = list( 'IVW'    = 'mr_ivw',
                                      'Median' = 'mr_weighted_median' ),
                      selected = 'mr_ivw' )
    })
}
















