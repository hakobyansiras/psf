library(data.table)
library(DT)
library(miniUI)
library(shiny)
library(Seurat)
library(ggplot2)
library(psf)
library(magick)
library(shinyjs)
library(visNetwork)
library(plotly)

shinyUI(
  navbarPage("PSF Spatial Browser",
             tabPanel("App",
                      fluidPage(
                        tags$script('
        $(document).on("keydown", function (e) {
           if(e.which == 17) {
              Shiny.onInputChange("multiple_choice", "yes");
           }
        });
      '),
                        #### hover dialog design ####
                        tags$head(tags$style('
                           #hover_text {
                             position: absolute;
                             max-width: 600px;
                             min-width: 0px;
                             border-radius: 5px;
                             z-index: 100;
                             color: white;
                             background-color: rgba(0, 0, 0, 0.8);
                             }
                             ')),
                        #### hovering script ####
                        tags$script('
                  $(document).ready(function(){
                  // id of the plot
                  $("#pathway_image").mousemove(function(e){ 
                  
                  // ID of uiOutput
                  $("#hover_text").show();
                  $("#hover_text").css({
                  top: (e.pageY + 5) + "px",
                  left: (e.pageX + 5) + "px"
                  });     
                  });     
                  });
        '),
                        uiOutput("hover_text"),
                        # titlePanel(span("PSF spatial browser", style = "color:#1f992f"), windowTitle = "PSF SP"),
                        column(3,
                               wellPanel(
                                 # selectInput(
                                 #   inputId = "annotation",
                                 #   label = "Select spot Annotation",
                                 #   choices = "",
                                 #   selected = "",
                                 #   selectize = FALSE,
                                 #   width = "100%"
                                 # ),
                                 selectInput(
                                   inputId = "cluster",
                                   label = "Select cluster",
                                   choices = "all", # c("all", as.character(unique_clusters)),
                                   selected = "all", # as.character(unique_clusters)[1],
                                   selectize = FALSE,
                                   width = "100%"
                                 ),
                                 tabsetPanel(
                                   id = "app_mode",
                                   tabPanel(
                                     title = "Top features",
                                     selectizeInput(
                                       inputId = "selected_feature",
                                       label = "Feature to visualize",
                                       choices = NULL,
                                       selected = NULL,
                                       width = "100%"
                                     ),
                                     htmlOutput('feature_stat')
                                   ),
                                   tabPanel(
                                     title = "Pathway vis",
                                     selectizeInput("selected_pathway", label = "Select pathway", choices = NULL, options = list(
                                       placeholder = 'Please select a pathway',
                                       onInitialize = I('function() { this.setValue(""); }')
                                     )),
                                     radioButtons(
                                       "node_coloring_type",
                                       label = "Select mapping value",
                                       choices = list("PSF" = 1, "FC" = 2),
                                       selected = 1,
                                       inline = TRUE
                                     )
                                   )
                                 ),
                                 htmlOutput('error_text'),
                                 width = "100%"
                               )
                        ),
                        column(
                          9,
                          plotlyOutput("spatial_plot", height = "550px", width = "650px"),
                          tabsetPanel(id = "net_vis_panels",
                                      tabPanel(title='KEGG', value = 1, useShinyjs(),
                                               imageOutput(
                                                 "pathway_image",
                                                 dblclick = "image_click",
                                                 brush = brushOpts(id = "image_brush", delayType = "debounce", resetOnNew = TRUE),
                                                 hover = hoverOpts(id = "image_hover", delay = 0)
                                               )        
                                      ),
                                      tabPanel(title='VisNet',
                                               visNetworkOutput("visnet", height = "800px")
                                      )
                          )
                        )
                      )
             ),
             tabPanel("Help",
                      tabsetPanel(
                        tabPanel("Welcome", includeMarkdown("help_pages/welcome.md")),
                        # tabPanel("Getting Started", includeMarkdown("help_pages/getting_started.md")),
                        tabPanel("User Interface Guide", includeMarkdown("help_pages/ui_guide.md")),
                        # tabPanel("Features", includeMarkdown("help_pages/features.md")),
                        tabPanel("Use Case", includeMarkdown("help_pages/use_case.md")),
                        tabPanel("FAQ & Troubleshooting", includeMarkdown("help_pages/faq.md")),
                        tabPanel("Resources", includeMarkdown("help_pages/resources.md")),
                        tabPanel("Contact", includeMarkdown("help_pages/contact.md"))
                      )
             )
  )
)
