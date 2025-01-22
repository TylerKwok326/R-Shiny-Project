library(shiny)
library(shinydashboard)
library(DT)
library(colourpicker)
library(shinythemes)
library(ggplot2)
library(pheatmap)
library(plotly)
library(gridExtra)
library(dplyr)
library(readr)
library(RColorBrewer)
library(stringr)
library(fgsea)
library(msigdbr)
library(GEOquery)
library(fgsea)
library(org.Hs.eg.db)

options(shiny.maxRequestSize = 30*1024^2) 

# UI Definition
ui <- dashboardPage(
  dashboardHeader(title = "Bioinformatics Analysis App"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("About", tabName = "about", icon = icon("info-circle")),
      menuItem("Sample Information", tabName = "samples", icon = icon("table")),
      menuItem("Counts Analysis", tabName = "counts", icon = icon("chart-bar")),
      menuItem("Differential Expression", tabName = "DE", icon = icon("dna")),
      menuItem("Gene Set Enrichment", tabName = "GSEA", icon = icon("project-diagram"))
    )
  ),
  dashboardBody(
    tabItems(
      # About tab
      tabItem(tabName = "about",
              fluidRow(
                box(width = 12,
                    h1("mRNA-Seq Analysis of Huntington's Disease Brain Tissue"),
                    h2("Study Overview"),
                    p("This application analyzes RNA sequencing data from a study examining transcriptional changes in Huntington's Disease (HD). The study focused on post-mortem BA9 brain tissue, comparing samples from HD patients with neurologically normal controls."),
                    h2("Dataset Details"),
                    tags$ul(
                      tags$li(strong("Tissue Source: "), "Human prefrontal cortex (BA9) region"),
                      tags$li(strong("Sample Size: "), "69 total samples",
                              tags$ul(
                                tags$li("20 Huntington's Disease patients"),
                                tags$li("49 neurologically normal controls")
                              )
                      ),
                      tags$li(strong("Technology: "), "High-throughput RNA sequencing"),
                      tags$li(strong("Total Genes Detected: "), "28,087"),
                      tags$li(strong("Differentially Expressed Genes: "), "5,480 (19% of detected genes, FDR<0.05)")
                    )
                )
              )
      ),
      
      # Samples tab
      tabItem(tabName = "samples",
              fluidRow(
                box(width = 3,
                    fileInput("samples_csvfile", "Load Sample CSV file",
                              accept = c(".csv"))
                ),
                box(width = 9,
                    tabsetPanel(
                      tabPanel("Summary", tableOutput("samples_summary")),
                      tabPanel("Data Table", dataTableOutput("samples_table")),
                      tabPanel("Visualizations", plotOutput("samples_plots"))
                    )
                )
              )
      ),
      
      # Counts tab
      tabItem(tabName = "counts",
              fluidRow(
                box(width = 3,
                    fileInput("counts_csvfile", "Load Counts CSV file",
                              accept = c(".csv")),
                    sliderInput("counts_slider_var", "Minimum variance percentile:",
                                min = 0, max = 100, value = 50),
                    sliderInput("counts_slider_num", "Minimum non-zero samples:",
                                min = 0, max = 69, value = 35)
                ),
                box(width = 9,
                    tabsetPanel(
                      tabPanel("Summary", tableOutput("counts_table")),
                      tabPanel("Filter Plots", plotOutput("counts_filter_plots")),
                      tabPanel("Heatmap", plotOutput("counts_heatmap")),
                      tabPanel("PCA", sidebarLayout(
                        sidebarPanel(
                          selectInput("PC1_choice", "Select PC1",
                                      choices = paste0("PC", 1:69)),
                          selectInput("PC2_choice", "Select PC2",
                                      choices = paste0("PC", 1:69)),
                          selectInput("pca_group", "Color points by:",
                                      choices = c("None")),
                          colourInput("pca_color1", "Group 1 color", "red"),
                          colourInput("pca_color2", "Group 2 color", "blue")
                        ),
                        mainPanel(
                          plotOutput("counts_PCA")
                        )
                      ))
                    )
                )
              )
      ),
      
      # DE tab
      tabItem(tabName = "DE",
              fluidRow(
                box(width = 3,
                    fileInput("DE_csvfile", "Load DE Results CSV",
                              accept = c(".csv")),
                    selectInput("DE_xchoice", "X-axis variable", 
                                choices = c("baseMean", "log2FoldChange", "pvalue", "padj")),
                    selectInput("DE_ychoice", "Y-axis variable",
                                choices = c("baseMean", "log2FoldChange", "pvalue", "padj")),
                    colourInput("DE_upcolor", "Upregulated genes color", "red"),
                    colourInput("DE_downcolor", "Downregulated genes color", "blue"),
                    colourInput("DE_nonsigcolor", "Non-significant genes color", "grey"),
                    sliderInput("DE_slider", "Significance threshold (-log10)",
                                min = -35, max = 0, value = -15, step = 1)
                ),
                box(width = 9,
                    tabsetPanel(
                      tabPanel("Results Table", tableOutput("DE_table")),
                      tabPanel("Volcano Plot", plotOutput("DE_volcano"))
                    )
                )
              )
      ),
      
      # GSEA tab
      tabItem(tabName = "GSEA",
              tabsetPanel(
                tabPanel("Pathway Visualization",
                         sidebarLayout(
                           sidebarPanel(
                             fileInput("GSEA_csvfile", "Load GSEA Results CSV",
                                       accept = c(".csv")),
                             sliderInput("GSEA_Bar", "Number of pathways to show:",
                                         min = 5, max = 50, value = 20)
                           ),
                           mainPanel(
                             plotOutput("GSEA_barplot", click = "barplot_click"),
                             conditionalPanel(
                               condition = "typeof input.barplot_click !== 'undefined'",
                               wellPanel(
                                 verbatimTextOutput("selected_pathway_info")
                               )
                             )
                           )
                         )
                ),
                tabPanel("Results Table",
                         sidebarLayout(
                           sidebarPanel(
                             sliderInput("GSEA_table_pvalue", "P-value threshold (-log10):",
                                         min = -10, max = 0, value = -5),
                             radioButtons("GSEA_pathways_choice", "Pathway direction:",
                                          choices = c("All", "Positive", "Negative")),
                             downloadButton("download_NES_table", "Download Results")
                           ),
                           mainPanel(
                             dataTableOutput("GSEA_results_table")
                           )
                         )
                ),
                tabPanel("Scatter Plot",
                         sidebarLayout(
                           sidebarPanel(
                             sliderInput("GSEA_scatter_pvalue", "P-value threshold (-log10):",
                                         min = -10, max = 0, value = -5)
                           ),
                           mainPanel(
                             plotOutput("GSEA_scatterplot")
                           )
                         )
                )
              )
      )
    )
  )
)
server <- function(input, output, session) {
  #---------------------------------------------------------Samples Tab------------------------------------------------------------
  samples_load_data <- reactive({
    req(input$samples_csvfile)
    file <- input$samples_csvfile
    if (is.null(file)) {
      return(NULL)
    } else {
      datafile <- read_csv(file$datapath) %>% 
        return()
    }
  })
  
  observe({
    req(samples_load_data())
    updateSelectInput(session, "pca_group",
                      choices = c(names(samples_load_data())))
  })
  
  samples_table_summary <- function(dataf) {
    # Create summary dataframe with column name and type
    summary_df <- data.frame(
      "Column Name" = names(dataf),
      "Type" = sapply(dataf, function(x) {
        if (is.numeric(x)) {
          if (all(x == floor(x), na.rm = TRUE)) "integer" 
          else "double"
        }
        else if (is.factor(x) || is.character(x)) "factor"
        else class(x)[1]
      }),
      "Mean (sd) or Distinct Values" = sapply(dataf, function(x) {
        if (is.numeric(x)) {
          # Handle cases with all NA values
          if (all(is.na(x))) {
            return("NA")
          }
          mean_val <- mean(x, na.rm = TRUE)
          sd_val <- sd(x, na.rm = TRUE)
          # Format based on the magnitude of the values
          if (abs(mean_val) >= 1000 || abs(sd_val) >= 1000) {
            sprintf("%.0f (+/- %.0f)", mean_val, sd_val)
          } else if (abs(mean_val) >= 100 || abs(sd_val) >= 100) {
            sprintf("%.1f (+/- %.1f)", mean_val, sd_val)
          } else {
            sprintf("%.2f (+/- %.2f)", mean_val, sd_val)
          }
        } else {
          # For categorical variables, count frequencies
          freq_table <- table(x)
          paste(names(freq_table), " (", freq_table, ")", collapse = ", ", sep = "")
        }
      }),
      stringsAsFactors = FALSE
    )
    
    # Add row count and column count as attributes
    attr(summary_df, "n_rows") <- nrow(dataf)
    attr(summary_df, "n_cols") <- ncol(dataf)
    
    return(summary_df)
  }
  
  library(ggplot2)
  library(dplyr)
  library(reshape2)
  
  samples_plot <- function(dataf) {
    colnames(dataf) <- gsub('-', '_', colnames(dataf))
    
    continuous_vars <- c("Age_of_death", "PMI", "RIN")
    
    plot_data <- melt(dataf[, c(continuous_vars, "Diagnosis")],
                      id.vars = "Diagnosis",
                      variable.name = "variable",
                      value.name = "value")
    
    ggplot(plot_data, aes(x = value, fill = Diagnosis)) +
      geom_density(alpha = 0.5) +
      facet_wrap(~variable, scales = "free", ncol = 2) +
      labs(title = "Distribution of Continuous Variables",
           x = "Value",
           y = "Density") +
      theme_minimal() +
      theme(plot.title = element_text(size = 16),
            strip.text = element_text(size = 12),
            axis.title = element_text(size = 12),
            axis.text = element_text(size = 10),
            legend.position = "bottom")
  }
  
  output$samples_summary <- renderTable({
    req(samples_load_data())  # Make sure data is loaded
    summary_table <- samples_table_summary(samples_load_data())
    # Add a header row with dataset dimensions
    attr_text <- sprintf("Dataset dimensions: %d rows × %d columns", 
                         attr(summary_table, "n_rows"), 
                         attr(summary_table, "n_cols"))
    summary_table <- rbind(c(attr_text, "", ""), summary_table)
    summary_table
  }, rownames = FALSE, width = "100%")  
  
  output$samples_table <- renderDataTable({
    samples_load_data()
  })
  
  output$samples_plots <- renderPlot({
    samples_plot(samples_load_data())
  })
  
  
  #---------------------------------------------------------Counts Tab------------------------------------------------------------
  counts_load_data <- reactive({
    req(input$counts_csvfile)
    file=input$counts_csvfile
    if (is.null(file))
    {return(NULL)} 
    else
    {datafile=read_csv(file$datapath)} %>% 
      return()
  })
  
  counts_filter_summary <- function(dataf, slider_var, slider_num){
    stats_df <- dataf
    stats_df$variance <- apply(dataf[,-1], 1, var)
    stats_df <- arrange(stats_df, by=variance)
    stats_df$rank <- rank(stats_df$variance)
    stats_df$nonzero <- rowSums(stats_df!=0)
    filtered <- filter(stats_df, rank >= nrow(stats_df)*(slider_var/100) & nonzero >= slider_num)
    num_samples = ncol(stats_df)-4
    num_genes = nrow(stats_df)
    num_genes_filtered = nrow(filtered)
    perc_genes_filtered = (num_genes_filtered/num_genes)*100
    num_genes_out = num_genes - num_genes_filtered
    perc_genes_out = (num_genes_out/num_genes)*100
    final_stats <- tibble(num_samples=num_samples,
                          num_genes=num_genes,
                          num_genes_filtered=num_genes_filtered,
                          percent_filtered=perc_genes_filtered,
                          num_genes_filtered_out=num_genes_out,
                          percent_filtered_out=perc_genes_out)
    return(final_stats)
  }
  
  counts_filter_plots <- function(dataf, slider_var, slider_num) {
    if (!is.null(dataf)) {
      # First plot - Median vs Variance
      # Convert to regular data frame and handle column selection
      plot_tib1 <- as.data.frame(dataf) %>%
        mutate(
          # Use numeric column indices instead of select
          Median = apply(.[, -1], MARGIN = 1, FUN = median),
          Variance = apply(.[, -1], MARGIN = 1, FUN = var)
        )
      perc_val <- quantile(plot_tib1$Variance, probs = slider_var/100)
      plot_tib1$thresh <- ifelse(plot_tib1$Variance >= perc_val, "TRUE", "FALSE")
      
      p1 <- ggplot(plot_tib1, aes(Median, Variance)) +
        geom_point(aes(color=thresh), alpha=0.75) +
        scale_color_manual(values = c("FALSE" = "black", "TRUE" = "blue")) +
        labs(title = 'Plot of Median vs Variance', 
             subtitle = "Genes filtered out are in red. X and Y axes are log-scaled.") +
        scale_y_log10() +
        scale_x_log10() +
        theme_bw() +
        theme(legend.position = 'bottom')
      
      # Second plot - Median vs Non-zero samples
      plot_tib2 <- as.data.frame(dataf)
      plot_tib2$Median <- apply(plot_tib2[, -1], MARGIN = 1, FUN = median)
      plot_tib2[plot_tib2 == 0] <- NA
      plot_tib2$no_zeros <- rowSums(is.na(plot_tib2))
      plot_tib2$thresh <- ifelse(plot_tib2$no_zeros <= slider_num, "TRUE", "FALSE")
      
      p2 <- ggplot(plot_tib2, aes(Median, no_zeros)) +
        geom_point(aes(color=thresh), alpha=0.75) +
        scale_color_manual(values = c("FALSE" = "black", "TRUE" = "blue")) +
        scale_x_log10() +
        labs(title = 'Plot of Median vs Number of Non-Zero genes',
             subtitle = "Genes filtered out are in red. X-axis is log scaled.") +
        theme_bw() +
        ylab('Number of samples with zero count') +
        theme(legend.position = 'bottom')
      
      grid.arrange(p1, p2, nrow = 1)
    }
  }
  
  counts_heatmap <- function(dataf, slider_var, slider_num) {
    stats_df <- dataf
    stats_df$variance <- apply(dataf[,-1], 1, var)
    stats_df <- arrange(stats_df, by=variance)
    stats_df$rank <- rank(stats_df$variance)
    stats_df$nonzero <- rowSums(stats_df!=0)
    filtered <- filter(stats_df, rank >= nrow(stats_df)*(slider_var/100) & nonzero >= slider_num)
    mat <- as.matrix(filtered[,2:70])
    rownames(mat) <- filtered$gene
    p <- pheatmap(log10(mat+1),
                  scale="row",
                  color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(11),
                  cluster_rows=TRUE,
                  cluster_cols=TRUE,
                  show_rownames=FALSE,
                  show_colnames=TRUE,
                  legend=TRUE
    )
  }
  
  counts_PCA <- function(counts_data, meta_data, PC1_choice, PC2_choice, group_var = NULL, color1, color2) {
    # Remove the first column if it's an index or gene names
    expr_mat <- as.matrix(counts_data[,-1])
    rownames(expr_mat) <- counts_data[[1]]  # save row names if needed
    expr_mat <- t(expr_mat)  # transpose for PCA
    
    # Perform PCA
    pca <- prcomp(expr_mat, center = TRUE, scale = FALSE)
    
    # Calculate variance explained
    var_explained <- (pca$sdev^2 / sum(pca$sdev^2)) * 100
    
    # Create plot data
    plot_data <- as.data.frame(pca$x)
    
    # Extract PC numbers from choices
    pc1_num <- as.numeric(gsub("PC", "", PC1_choice))
    pc2_num <- as.numeric(gsub("PC", "", PC2_choice))
    
    # Add grouping information if provided
    if (!is.null(group_var) && group_var != "None" && !is.null(meta_data)) {
      # Make sure rownames match between PCA and metadata
      if (all(rownames(plot_data) %in% meta_data$Sample)) {
        plot_data$group <- meta_data$Diagnosis[match(rownames(plot_data), meta_data$Sample)]
        
        p <- ggplot(plot_data, aes(x = plot_data[,pc1_num], 
                                   y = plot_data[,pc2_num],
                                   color = group)) +
          geom_point(size = 3) +
          scale_color_manual(values = c(color1, color2)) +
          labs(x = sprintf("PC%d (%.1f%%)", pc1_num, var_explained[pc1_num]),
               y = sprintf("PC%d (%.1f%%)", pc2_num, var_explained[pc2_num]),
               color = group_var) +
          theme_bw()
      } else {
        warning("Sample names don't match between counts and metadata")
        p <- ggplot(plot_data, aes(x = plot_data[,pc1_num], 
                                   y = plot_data[,pc2_num])) +
          geom_point(size = 3) +
          labs(x = sprintf("PC%d (%.1f%%)", pc1_num, var_explained[pc1_num]),
               y = sprintf("PC%d (%.1f%%)", pc2_num, var_explained[pc2_num])) +
          theme_bw()
      }
    } else {
      p <- ggplot(plot_data, aes(x = plot_data[,pc1_num], 
                                 y = plot_data[,pc2_num])) +
        geom_point(size = 3) +
        labs(x = sprintf("PC%d (%.1f%%)", pc1_num, var_explained[pc1_num]),
             y = sprintf("PC%d (%.1f%%)", pc2_num, var_explained[pc2_num])) +
        theme_bw()
    }
    
    return(p)
  }
  
  
  output$counts_table <- renderTable({counts_filter_summary(dataf = counts_load_data(), 
                                                            slider_var = input$counts_slider_var, 
                                                            slider_num=input$counts_slider_num)})
  
  output$counts_filter_plots <- renderPlot({counts_filter_plots(dataf = counts_load_data(),
                                                                slider_var = input$counts_slider_var, 
                                                                slider_num = input$counts_slider_num)})
  
  output$counts_heatmap <- renderPlot({counts_heatmap(dataf = counts_load_data(), 
                                                      slider_var = input$counts_slider_var,
                                                      slider_num = input$counts_slider_num)})
  
  output$counts_PCA <- renderPlot({
    req(counts_load_data())
    counts_PCA(
      counts_data = counts_load_data(),
      meta_data = samples_load_data(),
      PC1_choice = input$PC1_choice,
      PC2_choice = input$PC2_choice,
      group_var = input$pca_group,
      color1 = input$pca_color1,
      color2 = input$pca_color2
    )
  })
  #---------------------------------------------------------DE Tab------------------------------------------------------------
  DE_load_data <- reactive({
    req(input$DE_csvfile)
    file <- input$DE_csvfile
    if (is.null(file)) {
      return(NULL)
    } else {
      # Read the CSV file
      datafile <- read.csv(file$datapath, header = TRUE, sep = ",")
      
      # Check if the first column needs to be renamed
      # Only rename if it's not already named appropriately
      if (names(datafile)[1] %in% c("", "X", "x", "row.names")) {
        names(datafile)[1] <- "Gene"
      }
      
      return(datafile)
    }
  })
  
  #' Volcano plot
  
  # Make sure the volcano plot function can handle these values:
  DE_volcano_plot <- function(dataf, x_name, y_name, slider, up_color, down_color, nonsig_color, fc_threshold = 1) {
    # Create a factor for coloring points
    dataf <- dataf %>%
      mutate(
        regulation = case_when(
          !!sym(x_name) >= fc_threshold & !!sym(y_name) < 1*10^(slider) ~ "Upregulated",
          !!sym(x_name) <= -fc_threshold & !!sym(y_name) < 1*10^(slider) ~ "Downregulated",
          TRUE ~ "Not Significant"
        )
      )
    
    # Create the plot
    p <- ggplot(data = dataf, 
                aes(x = !!sym(x_name), 
                    y = -log10(!!sym(y_name)),
                    color = regulation)) + 
      geom_point() +
      theme_minimal() +
      theme(legend.position = "bottom") +
      scale_color_manual(
        values = c("Upregulated" = up_color, 
                   "Downregulated" = down_color, 
                   "Not Significant" = nonsig_color)
      ) +
      labs(x = x_name, 
           y = str_glue("-log10({y_name})"),
           color = "Regulation Status")
    
    return(p)
  }
  
  
  # Update the table filtering function accordingly:
  DE_draw_table <- function(dataf, slider) {
    # Make sure padj column exists
    if(!"padj" %in% names(dataf)) {
      return(dataf)  # Return unfiltered if no padj column
    }
    
    dataf %>%
      arrange(pvalue) %>% 
      filter(padj < 10^slider) %>%
      mutate(across(matches("^p(value|adj)$"), 
                    ~formatC(., digits = 2, format = "e"))) %>%
      return()
  }
  
  # return volcano output 
  output$DE_volcano <- renderPlot({
    DE_volcano_plot(
      dataf = DE_load_data(), 
      slider = input$DE_slider, 
      x_name = input$DE_xchoice, 
      y_name = input$DE_ychoice,
      up_color = input$DE_upcolor,
      down_color = input$DE_downcolor,
      nonsig_color = input$DE_nonsigcolor,
      fc_threshold = 1  
    )
  })
  
  #return table output
  output$DE_table <- renderTable({DE_draw_table(dataf = DE_load_data(), 
                                                slider = input$DE_slider)})
  
  #---------------------------------------------------------GSEA Tab------------------------------------------------------------
  # Data loading function aligned with UI component
  GSEA_load_data <- reactive({
    req(input$GSEA_csvfile)
    file <- input$GSEA_csvfile
    if (is.null(file)) {
      return(NULL)
    }
    data <- try(read_csv(file$datapath, col_names = TRUE))
    if (inherits(data, "try-error")) {
      stop("Error reading file. Please ensure it's a valid CSV.")
    }
    return(data)
  })
  
  # Visualization functions aligned with UI inputs
  GSEA_bar <- function(dataf, pathways_slider, wrap_width = 100) {
    filtered <- dataf %>% 
      arrange(padj) %>% 
      slice_head(n = pathways_slider) %>% 
      mutate(
        status = ifelse(NES > 0, 'Upregulated', 'Downregulated'),
        pathway = gsub("\\_", " ", pathway),
        pathway = str_wrap(pathway, width = wrap_width)
      ) %>%
      arrange(status) %>%
      mutate(pathway = factor(pathway, levels = unique(pathway)))
    
    ggplot(filtered, aes(x = reorder(pathway, NES), 
                         y = NES, 
                         fill = status)) +
      geom_col() +
      scale_fill_manual(
        values = c("Downregulated" = "#4575b4", "Upregulated" = "#d73027"),
        labels = c("Negative", "Positive")
      ) +
      coord_flip() +
      theme_bw() +
      theme(
        axis.text = element_text(size = 10),
        legend.position = "bottom"
      ) +
      labs(
        x = "Pathway", 
        y = "Normalized Enrichment Score",
        fill = "Regulation",
        title = "GSEA Pathway Analysis"
      )
  }
  
  GSEA_table <- function(dataf, pvalue_choice, type_pathways_choice) {
    data <- dataf %>% 
      mutate(
        status = ifelse(NES > 0, 'Positive', 'Negative'),
        pathway = str_replace_all(pathway, '_', ' ')
      ) %>% 
      filter(padj <= 10^(pvalue_choice)) %>%  
      dplyr::select(pathway, NES, pval, padj, size, status) %>%  # Added dplyr:: prefix
      mutate(across(where(is.numeric), round, digits = 3))
    
    if (type_pathways_choice != 'All') {
      data <- filter(data, status == type_pathways_choice)
    }
    
    return(data)
  }
  
  GSEA_scatter <- function(dataf, pvalue_choice) {
    dataf %>%
      mutate(
        significant = padj < 10^(pvalue_choice),  # This determines color
        neglog10padj = -log10(padj)  # This is correct for y-axis
      ) %>%
      ggplot(aes(x = NES,  # This is correct for x-axis
                 y = neglog10padj, 
                 color = significant)) +
      geom_point(alpha = 0.6) +
      scale_color_manual(
        values = c("TRUE" = "#d73027", "FALSE" = "grey50"),
        labels = c(
          paste('>', 10^(pvalue_choice)),
          paste('<', 10^(pvalue_choice))
        )
      ) +
      theme_bw() +
      labs(
        x = "Normalized Enrichment Score",
        y = "-log10(Adjusted P-value)",
        color = "Significant",
        title = "GSEA Results Scatter Plot"
      ) +
      theme(legend.position = "bottom")
  }
  
  # Outputs aligned with UI components
  output$GSEA_barplot <- renderPlot({
    req(GSEA_load_data())
    GSEA_bar(
      dataf = GSEA_load_data(), 
      pathways_slider = input$GSEA_Bar
    )
  }, height = function() {
    min(2500, max(900, input$GSEA_Bar * 30))
  })
  
  output$GSEA_results_table <- renderDataTable({
    req(GSEA_load_data())
    datatable(
      GSEA_table(
        dataf = GSEA_load_data(),
        pvalue_choice = input$GSEA_table_pvalue,
        type_pathways_choice = input$GSEA_pathways_choice
      ),
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        striped = TRUE
      )
    )
  })
  
  output$GSEA_scatterplot <- renderPlot({
    req(GSEA_load_data())
    GSEA_scatter(
      dataf = GSEA_load_data(), 
      pvalue_choice = input$GSEA_scatter_pvalue
    )
  })
  
  output$download_NES_table <- downloadHandler(
    filename = function() {
      paste0("GSEA_results_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write_csv(
        GSEA_table(
          GSEA_load_data(), 
          input$GSEA_table_pvalue, 
          input$GSEA_pathways_choice
        ), 
        file
      )
    }
  )
  output$selected_pathway_info <- renderPrint({
    req(input$barplot_click, GSEA_load_data())
    
    # Get the clicked y-position (which corresponds to the pathway)
    click <- input$barplot_click
    
    # Filter data to match what's shown in the plot
    plot_data <- GSEA_load_data() %>%
      arrange(padj) %>%
      slice_head(n = input$GSEA_Bar) %>%
      mutate(
        status = ifelse(NES > 0, 'Upregulated', 'Downregulated'),
        pathway = gsub("\\_", " ", pathway)
      ) %>%
      arrange(status)
    
    # Find the clicked pathway
    pathway_index <- round(click$y)
    if(pathway_index > 0 && pathway_index <= nrow(plot_data)) {
      selected_pathway <- plot_data[pathway_index, ]
      # Format the output
      cat("Pathway Details:\n")
      cat("Name:", selected_pathway$pathway, "\n")
      cat("NES:", round(selected_pathway$NES, 3), "\n")
      cat("P-value:", format.pval(selected_pathway$pval, digits = 3), "\n")
      cat("Adjusted p-value:", format.pval(selected_pathway$padj, digits = 3), "\n")
      cat("Size:", selected_pathway$size, "\n")
    }
  })
}
shinyApp(ui = ui, server = server) 