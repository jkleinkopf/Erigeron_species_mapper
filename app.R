# Erigeron (sensu lato) species mapper
# - Species table -> selectize dropdown (with color swatches)
# - Phylogeny (ggtree) below the dropdown; tips use the same colored circles as the map

library(shiny)
library(shinyWidgets)
library(shinylive)
library(httpuv)
library(tidyverse)
library(leaflet)
library(RColorBrewer)
library(htmltools)
library(bslib)
library(geojsonsf)
library(rmapshaper)
library(sf)
library(ape)
library(ggtree)

### Load Erigeron phylogeny
tree <- read.tree("erigeron_ingroup_all_phylo_2026-10-07.tree") # Close to done, but needs updating probably
tree$tip.label <- str_trim(gsub(pattern = "_", replacement = " ", tree$tip.label))

### Load species names data table
erigeron_sp <- read.csv("erigeron_sampling.csv")
erigeron_sp_subset <- erigeron_sp %>%
  arrange(scientificName) %>%
  select(scientificName, sectName, country) %>%
  mutate(scientificName = str_trim(scientificName)) %>%
  distinct() %>% 
  filter(scientificName %in% tree$tip.label)

erigeron_sp_subset[is.na(erigeron_sp_subset)] <- ""

### Load specimen data and clean up for mapping
gbif_parq <- nanoparquet::read_parquet("gbif_erigeron_acceptedNames.parquet") %>%
  filter(!is.na(species), species != "") %>%
  mutate(species = str_trim(species),
         scientificName2 = str_trim(scientificName2))

obs_tbl <- gbif_parq %>%
  rename(lat = decimalLatitude, lng = decimalLongitude) %>%
  filter(species %in% erigeron_sp_subset$scientificName)

### ONE shared, reproducible color mapping (map, dropdown and tree all use this)
sp_levels <- sort(unique(obs_tbl$species))
set.seed(42)  # fixed seed -> same colors every launch
sp_cols <- setNames(sample(scales::hue_pal(l = 60, c = 100)(length(sp_levels))), sp_levels)
pal <- colorFactor(palette = unname(sp_cols), domain = sp_levels, levels = sp_levels)

### Species table: one row per species (for the dropdown)
species_tbl <- erigeron_sp_subset %>%
  filter(scientificName %in% obs_tbl$species) %>%
  rename(species = scientificName, section = sectName) %>%
  group_by(species) %>%
  summarise(
    section = paste(unique(section[section != ""]), collapse = "; "),
    country = paste(unique(country[country != ""]), collapse = "; "),
    .groups = "drop"
  ) %>%
  mutate(color = unname(sp_cols[species]))

### HTML for each dropdown entry (color dot + species + section + country)
dot <- function(col) sprintf(
  "<span style='display:inline-block;width:12px;height:12px;border-radius:50%%;border:1px solid #000;margin-right:6px;vertical-align:middle;background:%s'></span>", col)
picker_content <- with(species_tbl, paste0(
  dot(color),
  "<i>", htmlEscape(species), "</i>",
  ifelse(section != "", paste0(" <span style='color:#555'>| ", htmlEscape(section), "</span>"), "")
))
picker_tokens <- with(species_tbl, paste(species, section, country))

### Build UI
ui <- fluidPage(
  fluidRow(
    column(10, titlePanel("Erigeron (sensu lato) species mapper")),
    column(2, div(style = "text-align:right; margin-top:20px;",
                  actionButton("aboutBtn", "About", icon = icon("circle-info"))))
  ),
  tags$head(
    tags$link(rel = "shortcut icon", href = "./mountain_and_flower.png")
  ),
  fluidRow(
    column(
      5,
      h4("Select one or more species"),
      pickerInput(
        "species_sel", label = NULL, width = "100%", multiple = TRUE,
        choices = species_tbl$species,
        choicesOpt = list(
          content = picker_content,   # HTML: color dot + species + section + country
          tokens  = picker_tokens     # extra text the search box matches on
        ),
        options = pickerOptions(
          liveSearch = TRUE,
          actionsBox = TRUE,                  # "Select All / Deselect All" buttons
          selectedTextFormat = "count > 2",   # "3 selected" instead of a long string
          countSelectedText = "{0} species selected",
          noneSelectedText = "Species | Section ",
          size = 12,
          container = "body"
        )
      ),
      actionButton("clearSelectionBtn", "Clear Selection"),
      plotOutput("tree_plot", height = "75vh")
    ),
    column(
      7,
      fluidRow(
        column(
          3,
          sliderInput("range", "Collection year range:", min = 1900, max = 2025, value = c(1900, 2025))
        ),
        column(
          2,
          switchInput(
            inputId = "activeFilter",
            label = "Show only imaged specimens",
            value = TRUE,
            onStatus = "success",
            offStatus = "danger"
          )
        )
      ),
      leafletOutput("map", height = "78vh")
    )
  )
)


# Server
server <- function(input, output, session) {
  
  ### "About" pop-up
  observeEvent(input$aboutBtn, {
    showModal(modalDialog(
      title = "About this app",
      HTML("<strong>How to use:</strong> Select one or more species from the drop-down menu (upper left). Species names and sections are shown and searchable (e.g., you can select all Erigeron subsection Tridactylia). 
            In the future, I hope to have an interactive and zoom-able phylogeny. Please let me know if you want to suggest any other improvements or ideas at jak566@cornell.edu!
            <br><br><strong>About:</strong> Data downloaded from GBIF shows preserved specimens collected between 1900-2025. Species not in the phylogeny are omitted here, and not all species are represented (i.e., rare taxa)
            <br><br><strong>GBIF data:</strong> GBIF.org (15 October 2025) GBIF Occurrence Download <a href='https://doi.org/10.15468/dl.n983dr' target='_blank'>https://doi.org/10.15468/dl.n983dr</a>
           <br><br><strong>Phylogeny:</strong> Phylogenetic relationships inferred from a coalescent species-tree analysis of 1030 Erigeron-specific single-copy loci (ms in preparation)
           <br><br><strong>Author:</strong> Joseph is a PhD candidate advised by Dr. Hannah Marx. <strong>Email</strong>: jak566@cornell.edu. <strong>Website</strong>: <a href='https://jkleinkopf.github.io/' target='_blank'>jkleinkopf.github.io</a>"),
      easyClose = TRUE,          # click outside or press Esc to close
      footer = modalButton("Close"),
      size = "m"
    ))
  })
  
  observeEvent(input$clearSelectionBtn, {
    updatePickerInput(session, "species_sel", selected = character(0))
  })
  
  # Reactive filter: only runs when the selection changes
  filtered_obs <- reactive({
    selected_species <- input$species_sel
    if (length(selected_species) == 0) return(NULL)
    
    dat <- obs_tbl %>%
      filter(species %in% selected_species, between(year, input$range[1], input$range[2]))
    
    if (input$activeFilter) dat %>% filter(mediaType == "StillImage") else dat
  })
  
  ### Phylogeny: selected species get the same fill color as their map circles;
  ### everything else is light grey.
  n_tips <- Ntip(tree)
  max_x  <- max(node.depth.edgelength(tree))
  
  output$tree_plot <- renderPlot({
    sel <- input$species_sel
    
    tipdat <- tibble(label = tree$tip.label) %>%
      mutate(
        is_sel   = label %in% sel,
        tip_fill = ifelse(is_sel, unname(sp_cols[label]), "transparent"),
        tip_line = ifelse(is_sel, "black", "transparent"),
        txt_col  = ifelse(is_sel, unname(sp_cols[label]), "black")
      )
    
    w <- session$clientData$output_tree_plot_width
    req(w)
    
    lab_size <- w / 250
    
    ggtree(tree, ladderize = T, layout = "circular", linewidth = 0.2) %<+% tipdat +
      geom_tiplab(aes(color = txt_col), fontface = "italic", size = lab_size, offset = max_x * 0.02) +
      geom_tippoint(aes(fill = tip_fill, color = tip_line), shape = 21, size = lab_size * 1.5, stroke = 0.2) +
      scale_fill_identity() +
      scale_color_identity()
  })
  
  ### Initialize leaflet map
  output$map <- renderLeaflet({
    static_options <- leafletOptions(
      zoomControl = TRUE, dragging = TRUE, minZoom = 2, maxZoom = 18,
      preferCanvas = TRUE
    )
    leaflet(options = static_options) %>%
      setView(lng = -105, lat = 45, zoom = 3) %>%
      addProviderTiles(providers$Esri.WorldTopoMap,    group = "Esri Topo (default)") %>%
      addProviderTiles(providers$OpenStreetMap.Mapnik, group = "OpenStreetMap") %>%
      addProviderTiles(providers$OpenTopoMap,          group = "OpenTopoMap (contours)") %>%
      addProviderTiles(providers$Esri.WorldImagery,    group = "Satellite (Esri)") %>%
      addLayersControl(
        baseGroups = c("Esri Topo (default)", "OpenStreetMap", "OpenTopoMap (contours)", "Satellite (Esri)"),
        options = layersControlOptions(collapsed = TRUE),
        position = "topright"
      )
  })
  
  ### Update points without reloading the base map
  observe({
    dat <- filtered_obs()
    
    leafletProxy("map", data = dat) %>%
      clearMarkers() %>%
      clearControls()
    
    if (!is.null(dat) && nrow(dat) > 0) {
      leafletProxy("map", data = dat) %>%
        addCircleMarkers(
          lng = ~lng, lat = ~lat, radius = 6, weight = 1, color = "black",
          fillColor = ~pal(species), fillOpacity = 0.65,
          popup = ~paste0("<i>", scientificName2, "</i>", "<br>", "<strong>Institution:</strong> ", institutionCode,
                          "<br>", "<strong>Year collected:</strong> ", year,
                          "<br><strong>GBIF Record:</strong> <a href='https://www.gbif.org/occurrence/", gbifID,
                          "' target = '_blank'> Link</a>") %>% lapply(htmltools::HTML),
          label = ~paste0("<i>", scientificName2, "</i>", "<br><strong>Click to view more info</strong>") %>%
            lapply(htmltools::HTML)
        ) %>%
        addLegend(position = "bottomright", pal = pal, values = ~species,
                  title = "Selected Species", opacity = 1)
    }
  })
}

shinyApp(ui, server)
