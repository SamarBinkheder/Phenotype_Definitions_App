#This is the main one

library(shiny)
library(shinydashboard)
library(DBI)
library(pool)
library(DT)
library(bslib)
library(tidyverse)
library(glue)
#library(shinyWidgets)

# Connect to the MySQL database (replace with your database details)
db <- config::get("azuredb")
pool <- dbPool(
  drv = RMariaDB::MariaDB(), 
  dbname = db$database,
  host = db$server,
  port = db$port,
  username = db$uid,
  password = db$pwd,
  ssl.ca = db$ssl_ca,
  sslmode = "require")

tbl <- db$table
onStop(function() {poolClose(pool)})

data <- read_csv("phenotypes.csv",col_names = F)

# Define UI
ui <- page_sidebar(
  #Select theme 
  theme = bs_theme(bootswatch = "journal"),
  title = "Clinical Phenotype Knowledgebase (CliPheKB)",

  sidebar = sidebar(
    selectizeInput('query', choices = NULL, label = h6("Choose a phenotype:")),
    checkboxInput("term_boundaries", "Apply Exact Match", value = FALSE),  # Add the checkbox for term boundaries
    tags$hr(),  # Add a horizontal line between buttons
    checkboxGroupInput("icons", h6("Filter by Phenotype Keywords:"),choiceNames = list(" ICD", " CPT", " SNOMED", "Laboratory","Medication"), choiceValues = list("ICD", "CPT", "SNOMED", "Laboratory","Medication")),
    selectInput("keyword_filter", h6("Keywords Boolean Operator:"), choices = c("OR", "AND"), selected = "OR"),
    card(textOutput("txt")),
    tags$hr(),  # Add a horizontal line between buttons
    actionButton("submit", "Search / Apply Filters", class = "btn-info"),
    actionButton("clear_results", "Reset Results", class = "btn-info"),
    downloadButton("download_csv", "Download CSV", class = "btn-info")
  
  ),
  
  card( card_header(class ="bg-dark" , "Definition"), 
        fluidRow(column(width = 12, h4(""), 
        div(style = "font-size: 16px; font-weight: bold;",
        textOutput("result_count"))), column(width = 12, DTOutput("datatable"))),
        )
  
)


# Define server logic
server <- function(input, output, session) {
  updateSelectizeInput(session, 'query', choices = data$X1, server = TRUE)
  
  observeEvent(input$submit, {

    # Get the selected keywords
    keywords <- input$icons
    
    query <- glue_sql("SELECT * FROM {`tbl`}  WHERE {`tbl`}.term = {input$query}", .con = pool)
    
    query_result <- dbGetQuery(pool, query)
    
    query_result$original_sentence <- query_result$original_sentence %>% map(~ .x[!.x == 00]) %>% map_chr(rawToChar)
    
    query_result$original_sentence <- iconv(query_result$original_sentence,"WINDOWS-1252","UTF-8")
    
    # Apply the space filter to query_result based on term boundaries checkbox
    if (input$term_boundaries) {
      term_pattern <- paste0("\\b", input$query, "\\b")
      query_result <- query_result %>%
        filter(str_detect(query_result$original_sentence, regex(term_pattern, ignore_case = TRUE)))
    }
    
    # Filter the query result based on the selected keywords
    if (length(keywords) > 0) {
      if (input$keyword_filter == "OR") {
        query_result <- query_result %>% 
          filter(reduce(lapply(keywords, function(x) str_detect(str_to_lower(query_result$original_sentence), str_to_lower(x))), `|`))
      } else {
        query_result <- query_result %>% 
          filter(reduce(lapply(keywords, function(x) str_detect(str_to_lower(query_result$original_sentence), str_to_lower(x))), `&`))
      }
    }
    
    #### Save file for download
    query_result_save <- reactive({ query_result %>% select(PMID = pmid, Phenotype =term, term_id,Definition = original_sentence)})
    ####
    
    
    # Results and Links for Pubmed and replaced PMID by nothing in the PMID column
    df_reactive <- reactive({ query_result %>% 
        select(PMID = pmid, Phenotype =term, term_id,start=term_start_position, end =term_end_position, Definition = original_sentence) %>%
        mutate(PMID=paste0("<a href='https://pubmed.ncbi.nlm.nih.gov/", gsub("PMID", "", PMID), "/' target='_blank'>",  gsub("PMID", "", PMID),  "</a>")) %>%
        mutate(Definition = map_chr(Definition,~gsub(str_to_lower(str_sub(Definition,start=start,end=end)),
                                                     paste0("<span style='background-color: #FFFF00'>",
                                                     str_sub(Definition,start=start,end=end),"</span>"),.x,ignore.case = TRUE, useBytes = T)))
    })
    
    # Display the table with highlighted words
    output$datatable <- renderDT({datatable(df_reactive() %>% 
                                              as_tibble() %>% 
                                              select(PMID, Phenotype, Definition),escape = F, options = list(searching = TRUE, paging = TRUE, lengthMenu = c(10, 50, 100, 200, 500)))})
    

    # Display the number of results and unique PMIDs
    output$result_count <- renderText({
      n_results <- nrow(df_reactive())
      n_pmids <- length(unique(df_reactive()$PMID))
      if (n_results == 1) {paste("1 result found in", n_pmids, "article")} 
      else {paste(n_results, "results found in", n_pmids, "articles")} 
      }) 
    
  
    
    # Display the selected keywords
    output$txt <- renderText({
      if (length(keywords) > 0) {
        if (input$keyword_filter == "OR") {
          icons <- paste(keywords, collapse = " OR ")
        } else {
          icons <- paste(keywords, collapse = " AND ")
        }
        paste("The results are filtered by:",input$query, " AND (", icons, ")")
      } else {
        ""
      }
    })
  
    ### Start Download part #######
    output$download_csv <- downloadHandler(
      filename = function() {
        paste("Phenotype_pheno_data_", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        query_result_save() %>% 
          as_tibble() %>% 
          select(PMID, Phenotype, Definition) %>% 
          mutate(PMID = gsub("PMID", "", PMID)) %>% # remove HTML tags from PMID column
          write.csv(file = file, row.names = FALSE)
      }
    )
    ### End download part######
    
    })
  
  
  
  #Reset Page
  observeEvent(input$clear_results, {
    # reset Page to the first option
    updateSelectizeInput(session, 'query', choices = data$X1, server = TRUE)
    # Clear filters here
    updateCheckboxGroupInput(session, "icons", selected = character(0))
    updateSelectInput(session, "keyword_filter", selected = "OR")
    
    # Update the table and result count
    df_reactive <- reactive({ data.frame() })
    output$datatable <- renderDT({datatable(df_reactive() %>% 
                                              as_tibble() %>% 
                                              select(PMID, Phenotype, Definition),escape = F, options = list(searching = FALSE, paging = TRUE,lengthMenu = c(10, 50, 100, 200)))})
    output$txt <- renderText("")
    output$result_count <- renderText("")
    
    # Trigger the search button click to update the results
    session$sendInputMessage("submit", list())
    session$sendInputMessage("apply_filter", list())

    
  })
  
  
}


# Run the Shiny app
shinyApp(ui, server)


