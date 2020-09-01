library( shiny )
library( shinydashboard )
source( 'results_analysis.R' )

# all_phenotypes  =  read_rds( '../data/all_phenotypes_combined.rds' )

ui  =  dashboardPage(
    dashboardHeader(
        title = 'MR results'
    ),
    dashboardSidebar(
        radioButtons( 'mr_method_x',
                      'X MR method',
                      inline = TRUE,
                      choices = list( 'IVW'    = 'mr_ivw',
                                      'Median' = 'mr_weighted_median' ),
                      selected = 'mr_ivw' ),
        radioButtons( 'mr_method_y',
                      'Y MR method',
                      inline = TRUE,
                      choices = list( 'IVW'    = 'mr_ivw',
                                      'Median' = 'mr_weighted_median' ),
                      selected = 'mr_ivw' ),
        radioButtons(
            'mrtype',
            'Plot type',
            choices = list('Compare exposures'         = '2exp',
                           'Male vs. Female'           = 'mvf',
                           'Bidirectional MR'          = 'bimr',
                           'PC loadings'               = 'pcbar',
                           'Tissue heatmap'            = 'fuma_heatmap',
                           'Pascal heatmap'            = 'pascal_heatmap' ),
            selected = '2exp'
        ),
        checkboxGroupInput(
            'plot_options',
            'Plot options',
            choices = list(
                'Difference of effects' = 'separate_xy',
                'Color by subcategory'  = 'show_subcategories',
                'Legend'                = 'showlegend',
                'Export image'          = 'export_image',
                'Flip 1:1 line'         = 'flip_line',
                'Add regression line'   = 'add_regression_line'
            ),
            selected = c( 'showlegend' )
        )
    ),
    dashboardBody(
        fluidRow(
            shinydashboard::box( width = 8,
                                 plotlyOutput( 'mrplot', height = '800px' ) ),
            shinydashboard::box( width = 4,
                                 title = NULL,
                                 conditionalPanel(
                                     condition = 'input.mrtype == "bimr"',
                                     fluidRow(
                                         column( 4,
                                                 radioButtons(
                                                     'sex',
                                                     h3( 'Sex (Betas)' ),
                                                     choices = list( 'X-sex meta'   = 'both_sexes',
                                                                     'Female' = 'female',
                                                                     'Male'   = 'male' ),
                                                     selected = 'both_sexes',
                                                     inline = FALSE
                                                 )
                                         ),
                                         column( 4,
                                                 radioButtons(
                                                     'sex_pca',
                                                     h3( 'Sex (PCA)' ),
                                                     choices = list( 'Both'   = 'both_sexes',
                                                                     'Female' = 'female',
                                                                     'Male'   = 'male' ),
                                                     selected = 'female',
                                                     inline = FALSE
                                                 )
                                         ),
                                         column( 4,
                                                 radioButtons(
                                                     'sex_ivs',
                                                     h3( 'Sex (IVs)' ),
                                                     choices = list( 'Both'   = 'both_sexes',
                                                                     'Female' = 'female',
                                                                     'Male'   = 'male' ),
                                                     selected = 'both_sexes',
                                                     inline = FALSE
                                                 )
                                         )
                                     )
                                 ),
                                 conditionalPanel(
                                     condition = 'input.mrtype == "mvf"',
                                     fluidRow(
                                         column( 7,
                                                 selectInput(
                                                     'exposure_category',
                                                     'Exposure category',
                                                     choices = list( 'Body'       = 'body',
                                                                     'Biomarkers' = 'disease_proxy',
                                                                     'Diet'       = 'diet',
                                                                     'Diseases'   = 'disease',
                                                                     'Lifestyle'  = 'lifestyle' ),
                                                     selected = 'body'
                                                 )
                                         ),
                                         column( 5,
                                                 radioButtons(
                                                     'exposure_type',
                                                     'Type',
                                                     choices = list( 'PC'    = 'p',
                                                                     'Trait' = 't',
                                                                     'DXA'   = 'dxa' ),
                                                     selected = 'p',
                                                     inline = TRUE
                                                 )
                                         )
                                     ),
                                     conditionalPanel(
                                         condition = 'input.exposure_type == "p"',
                                         numericInput(
                                             'pc_number',
                                             'PC number',
                                             value = 1
                                         )
                                     ),
                                     conditionalPanel(
                                         condition = 'input.exposure_type == "t"',
                                         uiOutput( 'trait_select' )
                                     ),
                                     conditionalPanel(
                                         condition = 'input.exposure_type == "dxa"',
                                         uiOutput( 'dxa_select' )
                                     )
                                 ),
                                 
                                 conditionalPanel(
                                     condition = 'input.mrtype == "2exp"',
                                     fluidRow(
                                         h4( 'X axis' ),
                                         column( 3,
                                                 radioButtons(
                                                     'sexx_exposure',
                                                     'Sex exposure',
                                                     choices = list( 'X-sex meta'   = 'both_sexes',
                                                                     'Female' = 'female',
                                                                     'Male'   = 'male' ),
                                                     selected = 'both_sexes'
                                                 )
                                         ),
                                         column( 3,
                                                 radioButtons(
                                                     'sexx_outcome',
                                                     'Sex outcome',
                                                     choices = list( 'X-sex meta'   = 'both_sexes',
                                                                     'Female' = 'female',
                                                                     'Male'   = 'male' ),
                                                     selected = 'both_sexes'
                                                 )
                                         ),
                                         column( 3,
                                                 radioButtons(
                                                     'sexx_pca',
                                                     'Sex (PCA)',
                                                     choices = list( 'Both'   = 'both_sexes',
                                                                     'Female' = 'female',
                                                                     'Male'   = 'male' ),
                                                     selected = 'female'
                                                 )
                                         ),
                                         column( 3,
                                                 radioButtons(
                                                     'sexx_ivs',
                                                     'Sex (IVs)',
                                                     choices = list( 'Both'   = 'both_sexes',
                                                                     'Female' = 'female',
                                                                     'Male'   = 'male' ),
                                                     selected = 'both_sexes'
                                                 )
                                         )
                                     ),
                                     fluidRow(
                                         h4( 'Y axis' ),
                                         column( 3,
                                                 radioButtons(
                                                     'sexy_exposure',
                                                     'Sex exposure',
                                                     choices = list( 'X-sex meta'   = 'both_sexes',
                                                                     'Female' = 'female',
                                                                     'Male'   = 'male' ),
                                                     selected = 'both_sexes'
                                                 )
                                         ),
                                         column( 3,
                                                 radioButtons(
                                                     'sexy_outcome',
                                                     'Sex outcome',
                                                     choices = list( 'X-sex meta'   = 'both_sexes',
                                                                     'Female' = 'female',
                                                                     'Male'   = 'male' ),
                                                     selected = 'both_sexes'
                                                 )
                                         ),
                                         column( 3,
                                                 radioButtons(
                                                     'sexy_pca',
                                                     'Sex (PCA)',
                                                     choices = list( 'Both'   = 'both_sexes',
                                                                     'Female' = 'female',
                                                                     'Male'   = 'male' ),
                                                     selected = 'female'
                                                 )
                                         ),
                                         column( 3,
                                                 radioButtons(
                                                     'sexy_ivs',
                                                     'Sex (IVs)',
                                                     choices = list( 'Both'   = 'both_sexes',
                                                                     'Female' = 'female',
                                                                     'Male'   = 'male' ),
                                                     selected = 'both_sexes'
                                                 )
                                         )
                                     )
                                 ),
                                 conditionalPanel(
                                     condition = 'input.mrtype == "bimr" | input.mrtype == "2exp"',
                                     fluidRow(
                                         column( 7,
                                                 selectInput(
                                                     'category_1',
                                                     'X category',
                                                     choices = list( 'Body'       = 'body',
                                                                     'Biomarkers' = 'disease_proxy',
                                                                     'Diet'       = 'diet',
                                                                     'Diseases'   = 'disease',
                                                                     'Lifestyle'  = 'lifestyle' ),
                                                     selected = 'body'
                                                 )
                                         ),
                                         column( 5,
                                                 uiOutput( 'type_1_select' )
                                         )
                                     ),
                                     conditionalPanel(
                                         condition = 'input.type_1 == "p"',
                                         numericInput(
                                             'pc_number_1',
                                             'X PC number',
                                             value = 1
                                         )
                                     ),
                                     conditionalPanel(
                                         condition = 'input.type_1 == "t"',
                                         uiOutput( 'trait_1_select')
                                     ),
                                     conditionalPanel(
                                         condition = 'input.type_1 == "dxa"',
                                         uiOutput( 'dxa_1_select')
                                     ),
                                     fluidRow(
                                         column( 7,
                                                 selectInput(
                                                     'category_2',
                                                     'Y category',
                                                     choices = list( 'Body'       = 'body',
                                                                     'Biomarkers' = 'disease_proxy',
                                                                     'Diet'       = 'diet',
                                                                     'Diseases'   = 'disease',
                                                                     'Lifestyle'  = 'lifestyle' ),
                                                     selected = 'body'
                                                 )
                                         ),
                                         column( 5,
                                                 uiOutput( 'type_2_select' )
                                         )
                                     ),
                                     conditionalPanel(
                                         condition = 'input.mrtype == "2exp"',
                                         conditionalPanel(
                                             condition = 'input.type_2 == "p"',
                                             numericInput(
                                                 'pc_number_2',
                                                 'Y PC number',
                                                 value = 1
                                             )
                                         ),
                                         conditionalPanel(
                                             condition = 'input.type_2 == "t"',
                                             uiOutput( 'trait_2_select')
                                         ),
                                         conditionalPanel(
                                             condition = 'input.type_2 == "dxa"',
                                             uiOutput( 'dxa_2_select')
                                         )
                                         
                                     )
                                 ),
                                 conditionalPanel(
                                     condition = 'input.mrtype == "mvf" | input.mrtype == "2exp"',
                                     fluidRow(
                                         column( 7,
                                                 selectInput(
                                                     'outcome_category',
                                                     'Outcome category',
                                                     choices = list( 'Body'       = 'body',
                                                                     'Biomarkers' = 'disease_proxy',
                                                                     'Diet'       = 'diet',
                                                                     'Diseases'   = 'disease',
                                                                     'Lifestyle'  = 'lifestyle' ),
                                                     selected = 'disease'
                                                 )
                                         )
                                     )
                                 ),
                                 conditionalPanel(
                                     condition = 'input.mrtype == "fuma_heatmap"',
                                     selectInput(
                                         'fuma_tissue_dataset',
                                         h3( 'FUMA tissue dataset' ),
                                         choices = list( 'GTEx v8 30 tissues' = 'gtex_v8_ts_general_avg_log2TPM',
                                                         'GTEx v8 54 tissues' = 'gtex_v8_ts_avg_log2TPM',
                                                         'Brain development'  = 'bs_dev_avg_log2RPKM',
                                                         'Brain ages'         = 'bs_age_avg_log2RPKM' ),
                                         selected = 'gtex_v8_ts_general_avg_log2TPM'
                                     )
                                 ),
                                 conditionalPanel(
                                     condition = 'input.mrtype == "pcbar"',
                                     numericInput(
                                         'bar_pc_number',
                                         h3( 'PC number' ),
                                         value = 1
                                     )
                                 ),
                                 conditionalPanel(
                                     condition = 'input.mrtype == "pcbar"',
                                     checkboxGroupInput(
                                         'show_sexes',
                                         'Show sexes',
                                         choices = list(
                                             'Both combined' = 'both_sexes',
                                             'Female'        = 'female',
                                             'Male'          = 'male'
                                         ),
                                         selected = c( 'both_sexes', 'female', 'male' ),
                                         inline = TRUE
                                     )
                                 ),
                                 conditionalPanel(
                                     condition = 'input.mrtype == "pcbar"',
                                     checkboxGroupInput(
                                         'pcbar_options',
                                         '',
                                         choices = list(
                                             'Hide weak loadings' = 'filter_loadings',
                                             'Sort' = 'sort'
                                         ),
                                         selected = c( 'sort' ),
                                         inline = TRUE
                                     )
                                 )
            )
        )
    )
)













