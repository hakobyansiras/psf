library(shiny)
library(DT)
library(shinyWidgets)
library(shinyjs)
library(visNetwork)
library(plotly)

load("subtypes.RData")

shinyUI(
  
  fluidPage(
    
    tags$head(tags$script("
    // Enable navigation prompt
    window.onbeforeunload = function() {
        return 'Your changes will be lost!';
    };
    ")),
    
    tags$script('
      $(document).mousedown(function(ev){
        if(ev.which == 3)
        {
          Shiny.onInputChange("right_click", Math.random());
        }
      });
    '),
    
    tags$script('
      $(document).mouseup(function(ev){
        if(ev.which == 1)
        {
          Shiny.onInputChange("left_click", Math.random());
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
    
    #### multi selection script "Ctrl" ####
    tags$script('
    $(document).on("keydown", function (e) {
       if(e.which == 17) {
          Shiny.onInputChange("multiple_choice", "yes");
       }
    });
  '),
    tags$script('
    $(document).on("keyup", function (e) {
       if(e.which == 17) {
          Shiny.onInputChange("multiple_choice", "no");
       }
    });
  '),
    
    #### node delete hotkey "Delete" ####
    tags$head(tags$script('$(document).keyup(function(event){
                                      if(event.keyCode == 46){
                                      $("#node_delete").click();
                                      }});')),
    
    #### edge drawing hotkey "E" ####
    tags$head(tags$script('$(document).keyup(function(event){
                                      if(event.keyCode == 69){
                          $("#draw_edge").click();
                          }});')),
    
    # tags$head(tags$script('$(document).keyup(function(event){
    #                                   if(event.keyCode == 72){
    #                       $("#check_disconnected_nodes").click();
    #                       }});')),
    
    #### marks cleaner hotkey "Esc" ####
    tags$script('
      $(document).on("keydown", function (e) {
         if(e.which == 27) {
            Shiny.onInputChange("clear", Math.random());
         }
      });
    '),
    # tags$head(tags$script('document.getElementById("image").addEventListener("wheel", onWheel); 
    #                       function onWheel(e) {
    #                         Shiny.onInputChange("scroll_test", e.deltaY);
    # };'))
    
    
    #### right click dialog prevention script ####
    tags$head(tags$script('document.addEventListener("contextmenu", function(e) {

                           e.preventDefault();
                           }, false);'
    )),
    
    
    uiOutput("hover_text"),
    
    #### app right click dialog initiation triger ####
    tags$head(tags$style('
                       #dialog_content {
                         position: absolute;
                         width: 150px;
                         border-radius: 0px;
                         z-index: 100;
                         padding: 0px;
                         background-color: rgba(255,255,255, 0);
                         }
                         ')),
    tags$script('
                    $(document).ready(function(){
                    // id of the plot
                    $(document).mousedown(function(e){ 
                    if(e.which == 3)
                    {
                    // ID of uiOutput
                    $("#dialog_content").css({
                    top: (e.pageY + 5) + "px",
                    left: (e.pageX + 5) + "px"
                    });
                    }
                    });     
                    });
                 '),
    
    #### right click dialogs ####
      useShinyjs(),
      hidden(
        div(id = "dialog_content",
            wellPanel(
              style = "padding: 0px; position:absolute; z-index:100; background-color: rgba(255,255,255, 0); width: 150px; border-radius: 0px;",
              div(style="display: inline-block; vertical-align:top; background-color: #ffffff; padding: 4px; width: 148px; border-radius: 0px; position: sticky; max-height: 28px;", prettyCheckbox("ingoing_edge", "Edge in", value = FALSE, shape = "curve", status = "success")),
              div(id = "dir_switch", actionButton("edge_direction_switch", "Reverse edge direction", style='text-align: left; width: 148px; border-style: none; padding:4px; border-radius: 0px;')),
              div(id = "edge_search_dialog", actionButton("edge_search_dialog", "Search edge", style='text-align: left; width: 148px; border-style: none; padding:4px; border-radius: 0px;')),
              div(id = "add_edge_button", actionButton("add_edge", "Add edge", style = 'text-align: left; width: 148px; border-style: none; padding:4px; border-radius: 0px;')),
              div(id = "connect_to_button", actionButton("connect_to", "Connect to", style = 'text-align: left; width: 148px; border-style: none; padding:4px; border-radius: 0px;')),
              div(id = "event_association", actionButton("associate_with_event_node", "Associate with event", style = 'text-align: left; width: 148px; border-style: none; padding:4px; border-radius: 0px;')),
              div(id = "transition_maker_button", actionButton("make_transition", "Switch to transition node", style = 'text-align: left; width: 148px; border-style: none; padding:4px; border-radius: 0px;')),
              div(id = "edge_delete_button", actionButton("edge_delete", "Remove edge", style='text-align: left; width: 148px; border-style: none; padding:4px; border-radius: 0px;')),
              # div(id = "add_new_node", actionButton("new_node", "Add node", style='text-align: left; width: 148px; border-style: none; padding:4px; border-radius: 0px;')),
              div(id = "delete_nodes", actionButton("node_delete", "Delete node(s)", style='text-align: left; width: 148px; border-style: none; padding:4px; border-radius: 0px;')),
              div(id = "edit_edge", actionButton("edit_edge", "Edit edge attrs", style='text-align: left; width: 148px; border-style: none; padding:4px; border-radius: 0px;')),
              div(id = "change_exp_value", actionButton("change_exp", "Change FC value", style='text-align: left; width: 148px; border-style: none; padding:4px; border-radius: 0px;')),
              div(id = "add_drug_button", actionButton("add_drug", "Apply drug", style='text-align: left; width: 148px; border-style: none; padding:4px; border-radius: 0px;')),
              div(id = "change_edge_weight_button", actionButton("change_edge_weight", "Change edge weight", style='text-align: left; width: 148px; border-style: none; padding:4px; border-radius: 0px;')),
              div(id = "color_scale_fixer", style="display: inline-block; vertical-align:top; background-color: #ffffff; padding: 4px; width: 148px; border-radius: 0px; position: sticky; max-height: 28px;", prettyCheckbox("fix_scale", "Fix color scale", value = FALSE, shape = "curve", status = "success"))
            )
        )
      ),
    
    column(2, 
           fixedPanel(
             titlePanel(span("PSF toolkit", style = "color:#1f992f"), windowTitle = "PSF toolkit"),
             wellPanel(
             style="background-color: #f2f2f2; height: 900px; width: 250px;",
             
             div(style="display: inline-block;vertical-align:top; width: 210px; word-wrap: break-word;", htmlOutput('pathway_name_and_source')),
             div(style="display: inline-block;vertical-align:top; width: 210px;", tags$hr(style="border-color: #1f992f;")),
             
             tabsetPanel(id = "app_mode",
               tabPanel(title = "Curation", 
                        div(style="display: inline-block;vertical-align:top;",
                            br(),
                            useShinyjs(),
                            div(style="display: inline-block;vertical-align:top; width: 210px;", selectizeInput("pathway_source", label = "Pathway source", choices = c("Pathway collection", "DF", "Editor"), selected = "Pathway collection")),
                            div(id = "collection_loader", style="display: inline-block;vertical-align:top; width: 210px;", fileInput("load_collection", label = "Load collection", accept = ".zip")),
                            # htmlOutput('pathway_loading_error_text'),
                            div(id = "pathway_name_selector", style="display: inline-block;vertical-align:top; width: 210px;", selectizeInput("selected_pathway", label = "Add from kegg", choices = "", selected = "")),
                            div(id = "pathway_source_radio", radioButtons("pathway_source_selector", label = NULL, choices = list("KEGG" = 1, "Local" = 2), selected = 2, inline = T)),
                            actionButton("load_pathway", label = "Load")
                            
                        ),
                        
                        br(),
                        tags$hr(style="border-color: #1f992f;"),
                        
                        div(style="display: inline-block;vertical-align:top; width: 210px;",
                            dropdown(label = "Network editing",
                                     # div(style="display: inline-block;vertical-align:top;", actionButton("save", "Save work")),
                                     div(style="display: inline-block;vertical-align:top;", actionButton("draw_edge", "Draw_edge")),
                                     div(style="display: inline-block;vertical-align:top;", actionButton("check_disconnected_nodes", "Highlight disconnected node(s)")),
                                     div(style="display: inline-block;vertical-align:top;", actionButton("check_duplicated_nodes", "Highlight duplicated nodes"))
                            )
                        ),
                        
                        div(downloadButton("download_collection", label = "Download collection", style='display: block; margin: auto;')),
                        # div(style="display: inline-block;vertical-align:top; width: 210px;",
                        #     dropdown(label = "Download",
                        #              div(downloadButton("download_table", label = "Download pathway", style='display: block; margin: auto;')),
                        #              div(downloadButton("download_pathway", label = "Save pathway", style='display: block; margin: auto;'))
                        #     )
                        # ),
                        
                        div(style="display: inline-block;vertical-align:top; width: 210px;",
                            dropdown(label = "Presenting modes",
                                     div(style="display: inline-block;vertical-align:top;width: 130px;", checkboxInput("show_changes", label = "Show changes", value = FALSE)),
                                     div(style="display: inline-block;vertical-align:top;width: 130px;", checkboxInput("show_modified_nodes", label = "Show node changes", value = FALSE))
                            )
                        ),
                        div(style="display: inline-block;vertical-align:top; width: 210px;",
                            actionButton("undo", label = NULL, icon("undo", verify_fa = FALSE)),
                            actionButton("redo", label = NULL, icon("redo", verify_fa = FALSE))
                        )
                        # div(style="display: inline-block;vertical-align:top;", checkboxInput("ingoing_edge", "Edge in", width = "148px")),
                        # div(style="display: inline-block;vertical-align:top;", uiOutput("edge_search")) 
                        
                        ),
               tabPanel("Vis", 
                        br(),
                        div(style="display: inline-block;vertical-align:top; width: 210px;",
                            # dropdown(label = "Upload exp",
                            div(style="display: inline-block;vertical-align:top;", fileInput("file", label = NULL, placeholder = "Exp matrix"))
                            # )
                        ),
                        useShinyjs(),
                        hidden(
                          div(id = "sample_selection_type_input", style="display: inline-block;vertical-align:top; width: 210px;",
                              radioButtons("sample_selection_type", label = NULL, choices = list("All" = 1, "Sample" = 2, "Group" = 3), selected = 1, inline = T)
                          ),
                          div(id = "sample_selection_input", style="display: inline-block;vertical-align:top; width: 210px;",
                              selectizeInput("selected_sample", label = "Select sample", choices = "")
                          )
                        ),
                        hidden(
                          div(id = "vis_buttons",
                              div(style="display: inline-block;vertical-align:top; width: 210px;",
                                  div(style="display: inline-block; vertical-align:top; width: 60px; top: 20px", checkboxInput("log_norm_checkbox", label = "Log", value = TRUE)),
                                  div(style="display: inline-block; vertical-align:top;", radioButtons("node_coloring_type", label = "Select mapping value", choices = list("PSF" = 1, "FC" = 2), selected = 1, inline = T))
                              ),
                              # div(style="display: inline-block;vertical-align:top; width: 210px;",
                              #     actionButton("map_fc_values", label = "Map FC values")
                              # ),
                              # div(style="display: inline-block;vertical-align:top; width: 210px;",
                              #     actionButton("psf_run", label = "Calculate psf")
                              # ),
                              # div(style="display: inline-block;vertical-align:top; width: 210px;",
                              #     selectizeInput("node_influence", label = "Node partial influence", choices = "", selected = NULL)
                              # ),
                              #### will work on this part ####
                              # div(style="display: inline-block;vertical-align:top; width: 210px;",
                              #     actionButton("reset_changed_values", label = "Reset modifications")
                              # ),
                              div(style="display: inline-block;vertical-align:top; width: 230px;",
                                  dropdown(label = "Partial influence analysis",
                                           div(style="display: inline-block;vertical-align:top;", selectizeInput("influence_type", label = "Influence type", choices = c("None", "any", "+", "-"), selected = "None")),
                                           div(style="display: inline-block;vertical-align:top;", selectizeInput("influence_node_num", label = "Number of nodes", choices = c(1:5, "all"), selected = "1")),
                                           div(style="display: inline-block;vertical-align:top;", actionButton("get_partial_influence", label = "Highlight part inf. nodes"),
                                               htmlOutput('partial_inf_err')
                                           )
                                  )
                              )
                              # div(style="display: inline-block;vertical-align:top; width: 210px;",
                              #     downloadButton("download_signal_table", label = "Download signal table", style='display: block; margin: auto;')
                              # )
                          )
                        ),
                        
                        htmlOutput('fc_table_load_error')
                        )
             ),
             
             
             # div(style="display: inline-block;vertical-align:top; width: 210px;",
             #     dropdown(label = "Load from file",
             #              div(style="display: inline-block;vertical-align:top;", fileInput("kegg_data", label = h3("Load RData file")))
             #     )
             # ),
             
             # div(style="display: inline-block;vertical-align:top; width: 210px;", 
             #     selectizeInput("app_mode", label = "App mode", choices = c("Curation", "Vis"), selected = "Curation")
             #     ),
             
             conditionalPanel(
               condition = "input.app_mode == 'Curation'",
               
             ),
             conditionalPanel(
               condition = "input.app_mode == 'Vis'",
               
             )
             
           )
           )
    ),
    
    column(10, 
           #### pathway image ####
           fluidRow(
             
             #### edge adding panel ####
             useShinyjs(),
             hidden(
               div(id = "edge_panel",
                   absolutePanel(
                     top = 20, right = 20,
                     width = 600,
                     draggable = TRUE,
                     fixed = TRUE,
                     cursor = "auto",
                     wellPanel(
                       style = "background-color: rgb(217, 217, 217, 0.8); z-index:2000;",
                       actionButton("close_table", "X", style='position: absolute; right: 8px; top: 8px; text-align: center; width: 18px; height 15px; border-style: solid; border-radius: 10px; background-color: rgb(255, 0, 0, 0.9); font-size: 80%; padding: 0px; align :'),
                       div(style = 'position: absolute; left: 10px; top: 12px;', prettySwitch(inputId = "direction_side", label = "In", value = TRUE)),
                       div(style="display: inline-block; vertical-align:top; left: 80px; top: 0px; padding: 1px; width: 148px; position: absolute; max-height: 28px;", checkboxInput("direction_state", label = "Direction", value = FALSE)),
                       br(),
                       DT::dataTableOutput("edge_table"),
                       div(id = "add_selected_edge", actionButton("add_selected_edge", "Add edge", style='position: absolute; left: 10px; top: 82%;'))
                     )
                   )
               )
             ),
             # tags$head(tags$style(
             #   HTML('
             #        #edge_creating_panel {position: absolute; background-color: rgb(217, 217, 217, 0.8); height 300px;}')
             # )),
             hidden(
               div(id = "edge_creating_panel",
                   absolutePanel(
                     style="top: 77px; right: 20px; width: 400px; height: 240px; position: fixed; cursor: move; left: 459px; background-color: rgba(217, 217, 217, 0.8);",
                     top = 20, right = 20,
                     width = 400,
                     draggable = TRUE,
                     fixed = TRUE,
                     cursor = "auto",
                     class = "panel panel-default",
                     actionButton("close_panel", "X", style='position: absolute; right: 8px; top: 8px; text-align: center; width: 18px; height 15px; border-style: solid; border-radius: 10px; background-color: rgb(255, 0, 0, 0.9); font-size: 80%; padding: 0px; align :'),
                     hidden(div(id = "node_selecion_message", style = 'padding-left: 8px; padding-top: 8px;', HTML('<p style="font-size: 16px; color: #669900;"><b>Select nodes to connect with.</b></p>'))),
                     div(style = 'padding-left: 8px; padding-top: 8px;', htmlOutput('edge_adding_error')),
                     hidden(div(id = "image_edge_direction_switch", style = 'position: absolute; left: 255px; top: 8px; width: 50px;', prettySwitch(inputId = "image_edge_direction", label = "In", value = TRUE))),
                     div(style = "padding-right: 10px; padding-left: 10px; padding-top: 10px;",
                         fluidRow(
                           column(6, div(selectizeInput("subtype1", label = "Edge subtype1*", choices = edge_subtype1, selected = ""))),
                           column(6, div(style="position:absolute; top: 32px; left: 30px; outline: 1px solid; outline-offset: 5px;", htmlOutput("edge_info")))
                         )
                     ),
                     div(style = "padding-right: 10px; padding-left: 10px;",
                         fluidRow(
                           column(6, div(selectizeInput("subtype2", label = "Edge subtype2", choices = c("activation", "binding/association", "compound", "dephosphorylation", "dissociation", "expression", "glycosylation", "indirect effect", "inhibition", "phosphorylation", "reaction", "repression", "state change", "ubiquitination", ""), selected = ""))),
                           column(6, div(selectizeInput("interaction_source", "Edge source*", choices = c("image", "database", "literature", ""), selected = "", multiple = TRUE))),
                           column(6, div(numericInput("", "Edge weight*", value = 1, min = 0, max = 1)))
                         )
                     ),
                     fluidRow(
                       div(id = "add_edge", actionButton("create_edge", "Create edge", style='display: block; margin: auto;'))
                     )
                   )
               )
             ),
             hidden(
               div(id = "event_node_adding_panel",
                   absolutePanel(
                     style="top: 77px; right: 20px; width: 400px; position: fixed; cursor: move; left: 459px; background-color: rgba(217, 217, 217, 0.8);",
                     top = 20, right = 20,
                     width = 400,
                     draggable = TRUE,
                     fixed = TRUE,
                     cursor = "auto",
                     class = "panel panel-default",
                     actionButton("close_node_panel", "X", style='position: absolute; right: 8px; top: 8px; text-align: center; width: 18px; height 15px; border-style: solid; border-radius: 10px; background-color: rgb(255, 0, 0, 0.9); font-size: 80%; padding: 0px; align :'),
                     div(style = 'padding-left: 8px; padding-top: 8px;', htmlOutput('node_adding_error')),
                     div(style = "padding-right: 10px; padding-left: 10px; padding-top: 10px;",
                         fluidRow(
                           column(6, selectizeInput('event_node_name', label = 'Node name*', choices = NULL, selected = "", options = list(create = TRUE) )),
                           column(6, selectizeInput("node_source", "Node source*", choices = c("image", "literature", ""), selected = "", multiple = TRUE))
                         )
                     ),
                     fluidRow(
                       column(6, div(style = 'padding-left: 11px; padding-top: 8px;', prettyCheckbox("event_node_mode", "Add event node", value = FALSE, shape = "curve", status = "success"))),
                       column(6, div(actionButton("add_node", "Create node", style='display: block; margin: auto;')))
                     )
                   )
               )
             ),
             useShinyjs(),
             hidden(
               div(id = "exp_change_panel", 
                   absolutePanel(
                     # top = 20, right = 20,
                     width = 250,
                     draggable = TRUE,
                     fixed = T,
                     cursor = "auto",
                     wellPanel(
                       style = "background-color: rgb(217, 217, 217, 0.8); z-index:1;",
                       numericInput("exp_slider", label = "Exp value", value = 1, step = 0.001),
                       selectizeInput("node_function_selector", label = "Node function", choices = c("mean", "max", "min", "gm_mean", "product"), selected = "mean"),
                       actionButton("change_exp_submit", label = "Change attrs and calc PSF"),
                       actionButton("close_slider", "X", style='position: absolute; right: 8px; top: 8px; text-align: center; width: 18px; height 15px; border-style: solid; border-radius: 10px; background-color: rgb(255, 0, 0, 0.9); font-size: 80%; padding: 0px; align :')
                       # div(style = 'position: absolute; left: 10px; top: 12px;', prettySwitch(inputId = "direction_side", label = "In", value = TRUE)),
                       # div(style="display: inline-block; vertical-align:top; left: 80px; top: 0px; padding: 1px; width: 148px; position: absolute; max-height: 28px;", checkboxInput("direction_state", label = "Direction", value = FALSE)),
                       # br(),
                       # DT::dataTableOutput("edge_table"),
                       # div(id = "add_selected_edge", actionButton("add_selected_edge", "Add edge", style='position: absolute; left: 10px; top: 82%;'))
                     )
                   )
               )
             ),
             useShinyjs(),
             hidden(
               div(id = "edge_attr_editing_panel_vis",
                   absolutePanel(
                     # style="top: 77px; right: 20px; width: 270px; height: 200px; position: fixed; cursor: move; left: 459px; background-color: rgba(217, 217, 217, 0.8);",
                     # top = 20, right = 20,
                     width = 250,
                     draggable = TRUE,
                     fixed = TRUE,
                     cursor = "auto",
                     wellPanel(
                       style = "background-color: rgb(217, 217, 217, 0.8); z-index:2000;",
	                     actionButton("close_edge_attr_vis_panel", "X", style='position: absolute; right: 8px; top: 8px; text-align: center; width: 18px; height 15px; border-style: solid; border-radius: 10px; background-color: rgb(255, 0, 0, 0.9); font-size: 80%; padding: 0px; align :'),
	                     selectizeInput("vis_subtype", label = "Edge type", choices = c("activation", "inhibition"), selected = ""),
	                     numericInput('vis_edge_weight', label = "Edge weight", min = 0, max = 1, step = 0.001, value = 1),
	                     actionButton("edit_vis_edge", "Edit edge", style='display: block; margin: auto;')
                     )
                   )
               )
             ),
             
             tabsetPanel(id = "net_vis_panels",
               tabPanel(title='KEGG', value = 1, useShinyjs(),
                        div(id = "image",
                            # style="height:580px; overflow-y: scroll; overflow-x: scroll;",
                            htmlOutput('pathway_load_error')),
                        imageOutput("pathway_image",
                                    dblclick = "image_click",
                                    brush = brushOpts(id = "image_brush", delayType = "debounce",resetOnNew = TRUE),
                                    hover = hoverOpts(id = "image_hover", delay = 0) 
                        )
                        # plotlyOutput("sink_plot", inline = TRUE)
                        
               ),
               tabPanel(title='VisNet', value = 2,
                        div(style="display: inline-block;vertical-align:top; width: 150px; margin-top: 1px;", textInput("network_search", label = NULL, placeholder = "Gene, chemical ...")),
                        div(style="display: inline-block;vertical-align:top; width: 150px; margin-top: 1px;", actionButton("search_go", label = "Search", icon = icon("search", lib = "glyphicon"))),
                        
                        # includeHTML("static.html"),
                        
                        visNetworkOutput("visnet", height = "800px"),
                        
                        # div(style="background-image: url('http://rest.kegg.jp/get/hsa04062/image')", 
                        #     visNetworkOutput("visnet", height = "800px")
                        # ),
                        useShinyjs(),
                        hidden(
                          div(id = "node_attrs",
                              absolutePanel(
                                top = 500, right = 500,
                                width = 300,
                                draggable = TRUE,
                                fixed = TRUE,
                                cursor = "auto",
                                wellPanel(
                                  style = "background-color: rgb(217, 217, 217, 0.8); z-index:2000;",
                                  actionButton("close_node_attr_panel", "X", style='position: absolute; right: 8px; top: 8px; text-align: center; width: 18px; height 15px; border-style: solid; border-radius: 10px; background-color: rgb(255, 0, 0, 0.9); font-size: 80%; padding: 0px; align :'),
                                  textInput("node_label", label = "Node name*"),
                                  textInput("node_genes", label = "Node ids*"),
                                  selectizeInput("node_type", label = "Node type*", choices = c("gene", "compound", "map", "linker", "event", ""), selected = ""),
                                  selectizeInput("node_function", label = "Node psf function", choices = c("mean", "min", "max", "sum"), selected = "mean"),
                                  numericInput("node_width", label = "Node width", value = 46),
                                  numericInput("node_height", label = "Node height", value = 17),
                                  numericInput("node_fc", label = "Node fold change", value = 1),
                                  div(style = 'padding-left: 8px; padding-top: 8px;', htmlOutput('node_data_adding_error')),
                                  actionButton("submit_node_attrs", label = "Save")
                                )
                              )
                          )
                        ),
                        
                        useShinyjs(),
                        hidden(
                          div(id = "edge_attrs",
                              absolutePanel(
                                top = 500, right = 500,
                                width = 300,
                                draggable = TRUE,
                                fixed = TRUE,
                                cursor = "auto",
                                wellPanel(
                                  style = "background-color: rgb(217, 217, 217, 0.8); z-index:2000;",
                                  actionButton("close_edge_attr_panel", "X", style='position: absolute; right: 8px; top: 8px; text-align: center; width: 18px; height 15px; border-style: solid; border-radius: 10px; background-color: rgb(255, 0, 0, 0.9); font-size: 80%; padding: 0px; align :'),
                                  selectizeInput("edge_type", label = "Edge type", choices = c("ECrel", "GErel", "PCrel", "PPrel"), selected = "PPrel"),
                                  selectizeInput("subtype1_visnet", label = "Subtype 1*", choices = c("activation", "binding/association", "compound", "dephosphorylation", "dissociation", "expression", "glycosylation", "indirect effect", "inhibition", "phosphorylation", "reaction", "repression", "state change", "ubiquitination", ""), selected = ""),
                                  selectizeInput("subtype2_visnet", label = "Subtype 2", choices = edge_subtype1, selected = ""),
                                  selectizeInput("interaction_source_visnetwork", "Edge source*", choices = c("image", "database", "literature", ""), selected = "image", multiple = TRUE),
                                  numericInput("edge_wight", label = "Edge weight", value = 1),
                                  div(style = 'padding-left: 8px; padding-top: 8px;', htmlOutput('edge_data_adding_error')),
                                  actionButton("submit_edge_attrs", label = "Save")
                                )
                              )
                          )
                        ),
                        
                        tabsetPanel(id = "undernetwork_tables", 
                                    tabPanel(title = "Node table",
                                             dataTableOutput("network_node_table", height = "200px", width = "100%")),
                                    tabPanel(title = "Edge table",
                                             dataTableOutput("network_edge_table", height = "200px", width = "100%"))
                        )
                        )
             )
           )
           )
  )
)
