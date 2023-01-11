


output$markerGeneTable <- renderDataTable({
  dt <- sample_data()$markers
  # commoncols <- c('cluster','gene','name','pctPie','pctBar','avgBar','avg_log2FC','pct.diff','p_val','p_val_adj')
  # commoncols <- c(commoncols[1:3],setdiff(colnames(dt),commoncols),commoncols[-c(1:3)])
  dt %<>% dplyr::select('cluster','gene', everything())

  DT::datatable(
    data = dt,
    escape = F,
    rownames = F,
    filter = "top",
    plugins = "natural",
    extensions = c("Buttons"),
    options = list(
      scrollCollapse = T,
      autoWidth = T,
      paging = T,
      pageLength = 10,
      rowGroup = list(dataSrc = 0),
      fixedHeader = T,
      ## dom = "t",
      ## scrollY = "540px",
      dom = "Bfrtip",
      scrollX = T,
      scrollY = T,
      searchHighlight = TRUE,
      buttons = list(
        list(
          extend = "colvis",
          # columns = c(0:8),
          text = "Columns to show"
        ),
        list(
          extend = "collection",
          text = "Download",
          buttons = list(
            list(
              extend = "csv",
              filename = "Marker_genes",
              title = "Marker_genes"
            ),
            list(
              extend = "excel",
              filename = "Marker_genes",
              title = "Marker_genes"
            )
          )
        )
      ),
      columnDefs = list(
        list(
          className = "dt-center",
          targets = c(0,1,2)
        ),
        list(
          className = "dt-left",
          targets = c(3:8)
        )
      ),
      fnDrawCallback = htmlwidgets::JS("
            function() {
              HTMLWidgets.staticRender();
            }
          ")
    )
  ) %>%
    spk_add_deps() %>%
    formatStyle(
      columns = "cluster",
      color = styleEqual(
        unique(dt$cluster),
        CellFunTopic::scPalette2(length(unique(dt$cluster)))
      ),
      fontWeight = "bold"
    ) %>%
    formatStyle(
      columns = 'p_val_adj',
      background = styleColorBar(dt$p_val_adj,
                                 '#FC4E07', angle = -90),
      backgroundSize = '100% 50%',
      backgroundRepeat = 'no-repeat',
      backgroundPosition = 'center',
      fontWeight = "bold"
    ) %>%
    formatStyle(
      columns = 'p_val',
      background = styleColorBar(dt$p_val,
                                 '#00AFBB', angle = -90),
      backgroundSize = '100% 50%',
      backgroundRepeat = 'no-repeat',
      backgroundPosition = 'center',
      fontWeight = "bold"
    ) %>%
    formatStyle(
      columns = "avg_log2FC",
      color = styleInterval(seq(min(dt$avg_log2FC), max(dt$avg_log2FC), length.out = 49),
                            colorRampPalette(c("#F8BF76", "#DB8B0A"))(50)),
      background = styleColorBar(dt$avg_log2FC,
                                 "#F39C11", angle = -90),
      backgroundSize = "98% 18%",
      backgroundRepeat = "no-repeat",
      backgroundPosition = "bottom",
      fontWeight = "bold"
    ) %>%
    formatStyle(
      columns = "pct.diff",
      color = styleInterval(seq(min(dt$pct.diff), max(dt$pct.diff), length.out = 49),
                            colorRampPalette(c("#807DBA", "#3F007D"))(50)),
      background = styleColorBar(dt$pct.diff,
                                 "#6A51A3", angle = -90),
      backgroundSize = "98% 18%",
      backgroundRepeat = "no-repeat",
      backgroundPosition = "bottom",
      fontWeight = "bold"
    ) %>%
    formatStyle(
      columns = "pct.1",
      backgroundColor = styleInterval(seq(min(dt$pct.1), max(dt$pct.1), length.out = 49),
                                      colorRampPalette(c("#FFFFFF", "#B03C2D"))(50)),
      fontWeight = "bold"
    ) %>%
    formatStyle(
      columns = "pct.2",
      backgroundColor = styleInterval(seq(min(dt$pct.2), max(dt$pct.2), length.out = 49),
                                      colorRampPalette(c("#FFFFFF", "#DB8B0A"))(50)),
      fontWeight = "bold"
    ) %>%
    formatRound(
      columns = "avg_log2FC"
    ) %>%
    formatRound(
      columns = "pct.diff"
    ) %>%
    formatSignif(
      columns = c("p_val_adj", "p_val"),
      digits = 3
    )
})


getSelectedGenes <- function(){
  rids <- input$markerGeneTable_rows_selected
  dat <- sample_data()$markers

  if(length(rids) > 0){
    genes <- dat[rids, 'gene']
  }else{
    genes <- NULL
  }
  genes
}

output$geneList <- renderUI({
  genes <- getSelectedGenes()
  tags$span(
    pickerInput(
      inputId = "selectedgenes",
      label = "Selected:  ",
      choices = genes,
      selected = genes,
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
})




generatePlot <- function(genes){
  withProgress(message="Generating ", {
    # incProgress(1/3, detail="FeaturePlot ...")
    # plst <- lapply(genes, function(x){
    #   pp <- FeaturePlot(sample_data()$SeuratObj, reduction = 'umap', order=TRUE,
    #                     cols=c('lightgrey', pals::brewer.blues(10)),
    #                     coord.fixed = F,# min.cutoff='q1',
    #                     pt.size=1, combine=T, features = x) + Seurat::NoAxes()
    #   pp
    # })
    # p1 <- ggarrange(plotlist = plst, common.legend = T)

    incProgress(1/2, detail="Violin Plot ...")
    p2 <- Seurat::VlnPlot(sample_data()$SeuratObj, features=genes, stack =T, flip=T,
                          fill.by = "ident", cols=sample_data()$colors$groupColor) +
      NoLegend() +
      theme(axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.text.y = element_blank())

    incProgress(1/2, detail="Bubble Plot ...")
    p3 <- DotPlot(sample_data()$SeuratObj, features = genes) +
      scale_color_gradientn(colours=rev(pals::brewer.rdylbu(20)), guide = "colourbar") +
      theme(axis.text.y=element_text(angle=0, hjust=1, size=8),
            axis.text.x=element_text(angle=45, hjust=1, size=12),
            axis.title = element_blank()) + ggpubr::rotate()

  })
  # return(list(p1=p1,p2=p2,p3=p3))
  return(list(p2=p2,p3=p3))
}


updatePlot <- function(genes){
  pout <- generatePlot(genes)
  # output$feature_plot <- renderPlot({
  #   pout$p1
  # })
  output$vln_plot <- renderPlot({
    pout$p2
  })
  output$bubble_plot <- renderPlot({
    pout$p3
  })
}


output$geneTable <- renderUI({
  dataTableOutput("markerGeneTable")
})


observeEvent(input$btn1, {
  dataTableProxy('markerGeneTable') %>% selectRows(NULL)
})


updateSinglePlot <- function(gene) {
  output$singleGene_featureplot <- renderPlot({
    pp1 <- FeaturePlot(sample_data()$SeuratObj, reduction = 'umap', order=TRUE,
                     cols=c('lightgrey', pals::brewer.blues(10)),
                     coord.fixed = F,# min.cutoff='q1',
                     pt.size=1, combine=T, features = gene) + NoAxes()
    pp1
  })
  output$singleGene_violinplot <- renderPlot({
    pp2 <- VlnPlot(sample_data()$SeuratObj, features=gene, combine = T, pt.size=0, cols=sample_data()$colors$groupColor) +
      ggpubr::theme_pubr() + NoLegend() + theme(axis.text.x=element_text(angle=45, hjust=1))
    pp2
  })
}


observe({
  if (length(input$selectedgenes) == 0) {
    gene <- sample_data()$markers %>% dplyr::ungroup() %>%
      dplyr::slice_max(order_by = avg_log2FC, n=1, with_ties = F) %>% dplyr::pull(gene)
    updateSinglePlot(gene)
  }
})


observeEvent({
  input$btn2
}, {
  if (length(input$selectedgenes) == 1) {
    gene <- input$selectedgenes

    updateTabsetPanel(
      session = session,
      inputId = "tabset",
      selected = "single"
    )
    updateSinglePlot(gene)

  }else{
    show_alert(title="Warnings",
               type="warning",
               text="Please select only one gene")
  }
})



observeEvent({
  input$btn3
}, {
  if (length(input$selectedgenes)>1) {
    genes <- input$selectedgenes
    updateTabsetPanel(
      session = session,
      inputId = "tabset",
      selected = "multiple"
    )
    updatePlot(unique(genes))
  }else{
    sendSweetAlert(session=session, title="Warnings",
                   type="warning",
                   text="Please select >=2 genes in the above table to generate plots")
  }
})






