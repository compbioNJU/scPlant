
sample_data <- reactive({
  withProgress(value=1,
               message="Loading data",
               detail="please wait for a while ...",
               {
                 shinyjs::hide(id = "nw_tohide1")
                 shinyjs::hide(id = "nw_tohide3")
                 meanRAS <- sapply(split(names(Idents(SeuratObj)), Idents(SeuratObj)),
                                   function(cells) colMeans(rasMat[cells, ]))
                 meanRAS <- t(scale(t(meanRAS), center = T, scale=T))

                 meanExp <- suppressMessages(Seurat::AverageExpression(SeuratObj, assays = DefaultAssay(SeuratObj),
                                                                       slot = 'data')[[DefaultAssay(SeuratObj)]])
                 meanExp <- meanExp[colnames(rasMat), ] %>% as.matrix()
                 meanExp <- t(scale(t(meanExp), center = T, scale=T))

                 rssMat_long <- reshape2::melt(rssMat)
                 colnames(rssMat_long) <- c('regulon', 'celltype', 'specificScore')
                 rssMat_long$regulon <- as.character(rssMat_long$regulon)
                 rssMat_long$celltype <- as.character(rssMat_long$celltype)
                 rssMat_long <- rssMat_long[, c(2,1,3)]

                 markers <- list()
                 markers <- SeuratObj@misc$Allmarkers

                 groupColor <- setNames(CellFunTopic::scPalette2(length(unique(Idents(SeuratObj)))), unique(Idents(SeuratObj)))
                 colors <- list(groupColor=groupColor)

                 ldaOut <- SeuratObj@misc$ldaOut
                 pws <- ID2Description(SeuratObj, by = "GO")
                 betaDF <- tidytext::tidy(ldaOut, matrix = "beta") %>% dplyr::mutate(descrip=unname(pws[term]))
                 gammaDF <- tidytext::tidy(ldaOut, matrix = "gamma")

               }
  )
  list(rssMat=rssMat, rssMat_long=rssMat_long,
       tf_target=tf_target, meanRAS=meanRAS, meanExp=meanExp,
       SeuratObj=SeuratObj, markers=markers, colors=colors, ldaOut=ldaOut,
       betaDF=betaDF, gammaDF=gammaDF, pws = pws)
})



GSEAresulttable <- reactive({
  nn <- paste0("GSEAresult_", input$switchTerm)
  if (nn %in% names(sample_data()$SeuratObj@misc)) {
    dt2 <- sample_data()$SeuratObj@misc[[nn]]
    if (nn == "GSEAresult_GO") {
      dt2$QuickGO <- sprintf('<a href="https://www.ebi.ac.uk/QuickGO/term/%s"><i class=\"fa fa-external-link-alt\"></i> %s</a>', dt2$ID, dt2$ID)
      # dt2$QuickGO <- sprintf('<a href="https://www.ebi.ac.uk/QuickGO/term/%s" class=\"fa fa-external-link-alt\"> %s</a>', dt2$ID, dt2$ID)
    }
    dt2 <- dt2 %>% dplyr::mutate(`-log10(FDR)`=-log10(p.adjust)) %>% dplyr::select("cluster","ID", "Description", "-log10(FDR)", everything())
  } else {
    return(NULL)
    }

  DT::datatable(
    data = dt2,
    class = 'display',
    escape = FALSE,
    rownames = FALSE,
    filter = "top",
    plugins = "natural",
    extensions = c("RowGroup","Buttons"),
    callback = htmlwidgets::JS(paste0(
      "table.rowGroup().",
      ifelse(input$tableShowSetting, "enable()", "disable()"),
      ".draw();"
    )),
    options = list(
      scrollCollapse = T,
      autoWidth = T,
      orderClasses = T,
      paging = T,
      pageLength = 10,
      rowGroup = list(dataSrc = 0),
      fixedHeader = T,
      scrollX = T,
      scrollY = T,
      searchHighlight = TRUE,
      dom = "Bfrtip",
      fnDrawCallback = htmlwidgets::JS("
            function() {
              HTMLWidgets.staticRender();
            }
          "),
      buttons = list(
        list(
          extend = "colvis",
          columns = c(3:14),
          text = "Columns to show"
        ),
        list(
          extend = "collection",
          text = "Download",
          buttons = list(
            list(
              extend = "csv",
              filename = "enriched_pathways",
              title = "Enriched pathways"
            ),
            list(
              extend = "excel",
              filename = "enriched_pathways",
              title = "Enriched pathways"
            )
          )
        )
      ),
      columnDefs = list(
        list(
          visible = F,
          targets = c(5,12,13)
        ),
        list(
          className = "dt-center",
          targets = c(0,1,4:7,11,14)
        ),
        list(
          className = "dt-left",
          targets = c(3,8:10)
        ),
        list(
          width = "220px",
          className = "dt-left",
          targets = 2
        ),
        list(width = "20px", targets = c(4,5)),
        list(width = "15px", targets = c(11)),
        list(width = "100px", targets = c(14))
      )
    )
  ) %>%
    formatStyle(
    columns = "cluster",
    color = styleEqual(
      unique(dt2$cluster),
      sample_data()$colors$groupColor[unique(dt2$cluster)]
    ),
    fontWeight = "bold"
  ) %>%
    formatStyle(
      columns = "ONTOLOGY",
      backgroundColor = styleEqual(
        unique(dt2$ONTOLOGY),
        grDevices::adjustcolor(scPalette(length(unique(dt2$ONTOLOGY))),0.3)
      ),
      fontWeight = "bold"
    ) %>%
    formatStyle(
    columns = 'p.adjust',
    background = styleColorBar(dt2$p.adjust,
                               '#FC4E07', angle = -90),
    backgroundSize = '100% 50%',
    backgroundRepeat = 'no-repeat',
    backgroundPosition = 'center',
    fontWeight = "bold"
  ) %>%
    formatStyle(
    columns = 'pvalue',
    background = styleColorBar(dt2$pvalue,
                               '#00AFBB', angle = -90),
    backgroundSize = '100% 50%',
    backgroundRepeat = 'no-repeat',
    backgroundPosition = 'center',
    fontWeight = "bold"
  ) %>%
    formatStyle(
      columns = 'qvalues',
      background = styleColorBar(dt2$qvalues,
                                 "#F39C11", angle = -90),
      backgroundSize = '100% 50%',
      backgroundRepeat = 'no-repeat',
      backgroundPosition = 'center',
      fontWeight = "bold"
    ) %>%
    # formatStyle(
    #   columns = "qvalues",
    #   color = styleInterval(seq(min(dt2$qvalues), max(dt2$qvalues), length.out = 49),
    #                         colorRampPalette(c("#F8BF76", "#DB8B0A"))(50)),
    #   background = styleColorBar(dt2$qvalues,
    #                              "#F39C11", angle = -90),
    #   backgroundSize = "98% 18%",
    #   backgroundRepeat = "no-repeat",
    #   backgroundPosition = "bottom",
    #   fontWeight = "bold"
    # ) %>%
    formatStyle(
      columns = "-log10(FDR)",
      color = styleInterval(seq(min(dt2$`-log10(FDR)`), max(dt2$`-log10(FDR)`), length.out = 49),
                            colorRampPalette(c("#807DBA", "#3F007D"))(50)),
      background = styleColorBar(dt2$`-log10(FDR)`,
                                 "#6A51A3", angle = -90),
      backgroundSize = "98% 18%",
      backgroundRepeat = "no-repeat",
      backgroundPosition = "bottom",
      fontWeight = "bold"
    ) %>%
    formatStyle(
      columns = "NES",
      backgroundColor = styleInterval(seq(min(dt2$NES), max(dt2$NES), length.out = 49),
                                      colorRampPalette(c("#FFFFFF", "#B03C2D"))(50)),
      fontWeight = "bold"
      ) %>%
    formatStyle(
      columns = "enrichmentScore",
      backgroundColor = styleInterval(seq(min(dt2$enrichmentScore), max(dt2$enrichmentScore), length.out = 49),
                                      colorRampPalette(c("#FFFFFF", "#DB8B0A"))(50)),
      fontWeight = "bold"
    ) %>%
    formatSignif(
      columns = c("-log10(FDR)","enrichmentScore", "NES", "pvalue","p.adjust","qvalues"),
      digits = 3
    )
})

output$GSEAtable <- DT::renderDataTable(
  GSEAresulttable(),
  server = TRUE
)






output$grouping <- renderUI({
  nn <- paste0("GSEAresult_", input$switchTerm)
  if (nn %in% names(sample_data()$SeuratObj@misc)) {
    tags$span(
      dropdownButton(
        tags$h4(icon("eye"), "Options"),
        tags$hr(),
        materialSwitch(
          inputId = "tableShowSetting",
          label = tagList(icon("object-group"), "Grouping"),
          status = "danger",
          value = F
        ),
        circle = F,
        right = T,
        inline = T,
        status = "danger",
        icon = icon("gear"),
        size = "sm",
        width = "300px",
        tooltip = tooltipOptions(title = "Options", placement = "top")
      ),
      style = "float:right;"
    )
  } else {
    return(NULL)
  }
})



output$pathwayList <- renderUI({
  if (paste0("GSEAresult_", input$switchTerm) %in% names(sample_data()$SeuratObj@misc)) {
    rids <- input$GSEAtable_rows_selected
    GSEAresult <- sample_data()$SeuratObj@misc[[paste0("GSEAresult_", input$switchTerm)]]
    nn <- paste(GSEAresult[rids, "ID"], GSEAresult[rids, "Description"], sep = " ")
    dat <- GSEAresult[rids, "ID"]
    if (!is.null(dat)) {
      dat <- setNames(object = dat, nn)
    }
    tags$span(
      pickerInput(
        inputId = "selectedpw",
        label = "Selected:  ",
        choices = dat,
        selected = GSEAresult[rids, "ID"],
        options = list(
          `actions-box` = TRUE,
          size = 10,
          `selected-text-format`= "count",
          width = "300px"
        ),
        multiple = TRUE,
        inline = TRUE
      ),
      style = "float:left;"
    )
  } else {
    return(NULL)
  }
})






output$select_plotui2 <- renderUI({
  tagList(tags$b("No GSEA result of "), shinydashboardPlus::dashboardLabel(input$switchTerm, status = "info", style = "square"), tags$b("were found."))
})

observeEvent(input$switchTerm, {
  if (!(paste0("GSEAresult_", input$switchTerm) %in% names(sample_data()$SeuratObj@misc))) {
    shinyjs::hide(id = "pa")
    shinyjs::hide(id = "tohide1")
    shinyjs::show(id = "tohide2")
  } else {
    shinyjs::show(id = "pa")
    shinyjs::show(id = "tohide1")
    shinyjs::hide(id = "tohide2")
  }
})



observeEvent(input$button2, {
  dataTableProxy('GSEAtable') %>% selectRows(NULL)
})





pws <- eventReactive(input$button3, input$selectedpw)

inter <- reactiveValues(pwss=NULL)


observeEvent({
  input$button3
}, {
  if (length(input$selectedpw) > 1) {
    inter$pwss <- pws()
  }else{
    show_alert(title="Warnings",
               type="warning",
               text="Please select >=2 pathways to generate plots.")
  }
})


gseaHeatmap1 <- reactive({
  gseaHeatmap(sample_data()$SeuratObj, by = input$switchTerm, pathwayIDs = inter$pwss, toshow = input$toshow, topPath = input$topPath,
              colour = input$hmpcol, scale = input$scale, fontsize_row = input$fontsize_row)
})
# output$gseaheatmap <- renderPlot(gseaHeatmap1())
output$gseaheatmap <- renderPlot({
  nn <- paste0("GSEAresult_", input$switchTerm)
  if (nn %in% names(sample_data()$SeuratObj@misc)) {
    gseaHeatmap1()
  } else {
    return(NULL)
  }
})



# clustercorplot1 <- reactive({
#   if (input$similarity == "Jaccard index") {
#     clustercorplot_jaccard(sample_data()$SeuratObj, by = input$switchTerm, pathwayIDs = inter$pwss, vertex.size.cex=input$vertex.size.cex,
#                            vertex.label.cex=input$vertex.label.cex, edge.max.width=input$edge.max.width, color.use = sample_data()$colors$groupColor,
#                            vertex.label.dist=input$vertex.label.dist, link_threshold=input$link_threshold_clustercorplot)
#   } else {
#     clustercorplot(sample_data()$SeuratObj, by = input$switchTerm, pathwayIDs = inter$pwss, vertex.size.cex=input$vertex.size.cex,
#                    vertex.label.cex=input$vertex.label.cex, edge.max.width=input$edge.max.width, color.use = sample_data()$colors$groupColor,
#                    vertex.label.dist=input$vertex.label.dist, link_threshold=input$link_threshold_clustercorplot)
#   }
# })
#
# # output$clustercorplot <- renderPlot(clustercorplot1())
# output$clustercorplot <- renderPlot({
#   nn <- paste0("GSEAresult_", input$switchTerm)
#   if (nn %in% names(sample_data()$SeuratObj@misc)) {
#     clustercorplot1()
#   } else {
#     return(NULL)
#   }
# })





hierarchyplot_tree1 <- reactive({
  hierarchyplot_tree(sample_data()$SeuratObj, by = input$switchTerm, pathwayIDs = inter$pwss, topaths = input$topath1,
                     cluster_cutree_k = input$cluster_cutree_k, pathway_cutree_k = input$pathway_cutree_k,
                     vertex.size.cex = input$vertex.size.cex2, color.use.cluster = sample_data()$colors$groupColor,
                     vertex.label.cex = input$vertex.label.cex2,
                     edge.max.width = input$edge.max.width2, alpha.edge = input$alpha.edge)
})

# output$hierarchyplot_tree <- renderPlot(hierarchyplot_tree1())
output$hierarchyplot_tree <- renderPlot({
  if (paste0("GSEAresult_", input$switchTerm) %in% names(sample_data()$SeuratObj@misc)) {
    hierarchyplot_tree1()
  } else {
    return(NULL)
  }
})

# output$cosineByGSEA <- renderPlot({
#   Cosine_networkByGSEA2(sample_data()$SeuratObj, layout=igraph::layout_with_fr, cos_sim_thresh=input$cosineByGSEA_thr,
#                        p.adjust_thresh=0.05, node.cex=6, width_range=c(0.8, 4), vertex.label.dist=0.8,
#                        text_size = input$cosineByGSEA_text, cluster.color = sample_data()$colors$groupColor2)
# })

output$edgeBundling <- renderPlot({
  edge_bundling_GSEA(sample_data()$SeuratObj, link_threshold=input$edgeBundling_thr,
                      p.adjust_thresh=0.05, link_width=0.9, method=input$edgeBundling_method,
                      node.by='cellType', group.by="cell_type")
})


showOnePw <- function() {
  output$pathwayScatterplot <- renderPlot({
    nn <- paste0("GSEAresult_", input$switchTerm)
    if (nn %in% names(sample_data()$SeuratObj@misc)) {
      pathwayScatterplot(sample_data()$SeuratObj, by = input$switchTerm, pathwayID = interr$pwidd,
                         reduction = input$reduction1, colour = input$sctcol, pointsize = input$pointsize,
                         label = input$labelornot, label.size = input$label.size)
    } else {
      return(NULL)
    }
  })
}

pwid <- eventReactive(input$button4, input$selectedpw)

interr <- reactiveValues(pwidd=NULL)

observeEvent({
  input$button4
}, {
  if (length(input$selectedpw) == 1) {
    interr$pwidd <- pwid()
    showOnePw()
  }else{
    show_alert(title="Warnings",
               type="warning",
               text="Please select only one pathway.")
  }
})

observe({
  if (length(input$selectedpw) == 0) {
    interr$pwidd <- sample_data()$SeuratObj@misc$GSEAresult_GO %>%
      dplyr::slice_min(order_by = p.adjust, n=1, with_ties = F) %>% dplyr::pull(ID)
    showOnePw()
  }
})


output$pathwayScatterplotdl <- downloadHandler(
  filename = "pathwayScatterplot.pdf",
  content = function(file) {
    pdf(file, width = 13, height = 10)
    p <- pathwayScatterplot(sample_data()$SeuratObj, by = input$switchTerm, pathwayID = interr$pwidd,
                            reduction = input$reduction1, colour = input$sctcol, pointsize = input$pointsize,
                            label = input$labelornot, label.size = input$label.size)
    print(p)
    dev.off()
  },
  contentType = "image/pdf"
)


output$gseaHeatmapdl <- downloadHandler(
  filename = "gseaHeatmapplot.pdf",
  content = function(file) {
    pdf(file, width = 13, height = 10)
    p <- gseaHeatmap1()
    print(p)
    dev.off()
  },
  contentType = "image/pdf"
)



output$hierarchyplot_treedl <- downloadHandler(
  filename = "hierarchyplot.pdf",
  content = function(file) {
    pdf(file, width = 13, height = 10)
    hierarchyplot_tree(sample_data()$SeuratObj, by = input$switchTerm, pathwayIDs = inter$pwss, topaths = input$topath1,
                       vertex.size.cex = input$vertex.size.cex2,
                       vertex.label.cex = input$vertex.label.cex2,
                       edge.max.width = input$edge.max.width2, alpha.edge = input$alpha.edge)
    dev.off()
  },
  contentType = "image/pdf"
)


########################################################################
######################### Topic modeling #########################

output$tm_echt <- renderEcharts4r({
  topicNW3(sample_data()$betaDF, topn=input$tm_topn1, pws=sample_data()$pws)
})


# output$tm_topicProb <- renderPlot({
#   tt <- gsub("Topic ", "", input$tm_echt_clicked_data$source)
#   if (length(tt)<1) tt <- unique(sample_data()$betaDF$topic)[1]
#   topicProb(sample_data()$SeuratObj, reduction="umap", topic=tt, pointSize=0.1) +
#     Seurat::NoAxes() + theme(plot.title = element_text(colour = "black", size = 20, hjust = 0.5, face = "bold"))
# })
#
# output$tm_Topterms <- renderPlot({
#   tt <- gsub("Topic ", "", input$tm_echt_clicked_data$source)
#   if (length(tt)<1) tt <- unique(sample_data()$betaDF$topic)[1]
#   Topterms_Topic(sample_data()$betaDF, Topic = tt, topn = 20, text.size = 6) +
#     theme(plot.title = element_text(colour = "black", size = 20, hjust = 0.5, face = "bold"))
# })

showtopic <- function() {
  output$tm_topicProb <- renderPlot({
    topicProb(sample_data()$SeuratObj, reduction="umap", topic=ttt$topicid, pointSize=0.1) +
      Seurat::NoAxes() + theme(plot.title = element_text(colour = "black", size = 20, hjust = 0.5, face = "bold"))
  })

  output$tm_Topterms <- renderPlot({
    Topterms_Topic(sample_data()$betaDF, Topic = ttt$topicid, topn = 20, text.size = 6) +
      theme(plot.title = element_text(colour = "black", size = 20, hjust = 0.5, face = "bold"))
  })
}

ttt <- reactiveValues(topicid=NULL)

observe({
  if (is.null(ttt$topicid)) {
    ttt$topicid <- unique(sample_data()$betaDF$topic)[1]
    showtopic()
  }
})

observeEvent({
  input$submit
}, {
  ttt$topicid <- input$tm_topics
  showtopic()
})

observeEvent({
  input$tm_echt_clicked_data
}, {
  ttt$topicid <- gsub("Topic ", "", input$tm_echt_clicked_data$source)
  showtopic()
})


output$tm_sankey <- renderSankeyNetwork(
  plot_sankey(sample_data()$gammaDF, topn=1)
)


output$tm_heatmap <- renderPlot(
  # cluster_topic_hmp(sample_data()$ldaOut)
  cluster_topic_hmp2(sample_data()$ldaOut)
)


output$tm_cosine_network_cluster <- renderPlot(
  cosine_network_cluster(sample_data()$ldaOut,
                         layout="fr",
                         cos_sim_thresh=input$cos_sim_thresh,
                         radius=input$radiuscosnet)
)


output$tm_cosine_network_term <- renderPlot(
  cosine_network_term(sample_data()$SeuratObj,
                      cosine_cal_by = "Topic modeling",
                      topn = input$tm_cosnet_topn,
                      layout = input$tm_layout2,
                      cos_sim_thresh = input$cos_sim_thresh2,
                      radius = input$radiuscosnet2,
                      text_size = input$tm_text_size)
)

























