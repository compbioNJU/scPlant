


RSStable <- reactive({
  RSSdf <- sample_data()$rssMat_long

  DT::datatable(
    data = RSSdf,
    class = 'display',
    rownames = FALSE,
    filter = "top",
    plugins = "natural",
    extensions = c("Buttons"),
    options = list(
      scrollCollapse = T,
      autoWidth = T,
      paging = T,
      pageLength = 20,
      fixedHeader = T,
      scrollX = T,
      scrollY = T,
      searchHighlight = TRUE,
      dom = "Bfrtip",
      buttons = list(
        list(
          extend = "collection",
          text = "Download",
          buttons = list(
            list(
              extend = "csv",
              filename = "regulons",
              title = "regulons"
            ),
            list(
              extend = "excel",
              filename = "regulons",
              title = "regulons"
            )
          )
        )
      ),
      fnDrawCallback = htmlwidgets::JS("
            function() {
              HTMLWidgets.staticRender();
            }
          "),
      columnDefs = list(
        list(
          className = "dt-center",
          targets = c(0,1)
        ),
        list(
          width = "100px",
          className = "dt-left",
          targets = 2
        )
      )
    )
  ) %>%
    formatStyle(
      columns = "celltype",
      color = styleEqual(
        unique(RSSdf$celltype),
        scPalette(length(unique(RSSdf$celltype)))
      ),
      fontWeight = "bold"
    ) %>%
    formatStyle(
      columns = "specificScore",
      color = styleInterval(seq(min(RSSdf$specificScore), max(RSSdf$specificScore), length.out = 49),
                            colorRampPalette(c("#807DBA", "#3F007D"))(50)),
      background = styleColorBar(RSSdf$specificScore,
                                 "#6A51A3", angle = -90),
      backgroundSize = "98% 18%",
      backgroundRepeat = "no-repeat",
      backgroundPosition = "bottom",
      fontWeight = "bold"
    )
})

output$RSStable <- DT::renderDataTable(
  RSStable()
)


output$RSS_topnPlot <- renderPlot({
  topRegulons(rssMat=sample_data()$rssMat, topn = input$nw_topn1)
})


observeEvent(input$nw_btn2, {
  dataTableProxy('RSStable') %>% selectRows(NULL)
})


observeEvent(input$nw_btn1, {
  shinyjs::show(id = "nw_tohide1")

  rids <- input$RSStable_rows_selected
  regulons <- unique(sample_data()$rssMat_long[rids, 'regulon'])
  tf_target <- sample_data()$tf_target[sample_data()$tf_target$TF %in% regulons, ]

  output$targetTable <- DT::renderDataTable({
    DT::datatable(
      data = tf_target,
      class = 'display',
      rownames = FALSE,
      filter = "top",
      plugins = "natural",
      extensions = c("Buttons"),
      options = list(
        scrollCollapse = T,
        autoWidth = T,
        paging = T,
        pageLength = 15,
        fixedHeader = T,
        scrollX = T,
        scrollY = T,
        searchHighlight = TRUE,
        dom = "Bfrtip",
        buttons = list(
          list(
            extend = "collection",
            text = "Download",
            buttons = list(
              list(
                extend = "csv",
                filename = "targets",
                title = "targets"
              ),
              list(
                extend = "excel",
                filename = "targets",
                title = "targets"
              )
            )
          )
        ),
        fnDrawCallback = htmlwidgets::JS("
            function() {
              HTMLWidgets.staticRender();
            }
          "),
        columnDefs = list(
          list(
            className = "dt-center",
            targets = c(0,1)
          ),
          list(
            width = "100px",
            className = "dt-left",
            targets = 2
          )
        )
      )
    ) %>%
      formatStyle(
        columns = "TF",
        color = styleEqual(
          unique(tf_target$TF),
          scPalette(length(unique(tf_target$TF)))
        ),
        fontWeight = "bold"
      ) %>%
      formatStyle(
        columns = "importance_score",
        color = styleInterval(seq(min(tf_target$importance_score), max(tf_target$importance_score), length.out = 49),
                              colorRampPalette(c("#807DBA", "#3F007D"))(50)),
        background = styleColorBar(tf_target$importance_score,
                                   "#6A51A3", angle = -90),
        backgroundSize = "98% 18%",
        backgroundRepeat = "no-repeat",
        backgroundPosition = "bottom",
        fontWeight = "bold"
      )
  })

  output$target_topnPlot <- renderPlot({
    toptargets(tf_target, topn = input$nw_topn2, regulons = regulons)
  })

  if (length(regulons) > 1) {
    shinyjs::show(id = "nw_tohide3")
    output$ras_exp_hmpPlot <- renderPlot({
      ras_exp_hmp2(meanRAS = sample_data()$meanRAS, meanExp = sample_data()$meanExp, regulons)
    })
  } else {shinyjs::hide(id = "nw_tohide3")}

})


# observeEvent({
#   input$nw_btn3
# }, {
#   rids <- input$RSStable_rows_selected
#   if (length(rids) == 1) {
#     shinyjs::show(id = "nw_tohide2")
#     regulon <- sample_data()$rssMat_long[rids, 'regulon']
#     output$ras_exp_plot <- renderPlot({
#       ras_exp_scatter(sample_data()$SeuratObj, rasMat=sample_data()$rasMat, gene=regulon, reduction = 'umap')
#     })
#   }else{
#     show_alert(title="Warnings",
#                type="warning",
#                text="Please select only one regulon.")
#   }
# })



