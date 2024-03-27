
## ---- Design_diagram function
Design_diagram <- function(form, data, Filter = NA, Colour = NULL, direction = "vertical") {
    require(ggraph)
    require(igraph)
  if(is.null(Colour)) {
    Colour <- 'Blank'
    data <- data %>% mutate(Blank = 1)
  }

    form <- lme4::lFormula(form,
      data = data,
      control = lme4::lmerControl(
        check.nobs.vs.nlev = "ignore",
        check.nobs.vs.nRE = "ignore",
        check.nlev.gtr.1 = "ignore"
      )
    )
  ## add nested variables
  dat <- as.data.frame(data) #form$fr
  terms <- form$fr %>% terms()
  resp <- terms %>% attr("response")
  fixed_terms <- terms %>% attr("varnames.fixed") %>% `[`(-1)
  fixed_terms <- ifelse(length(fixed_terms) == 0, "", fixed_terms)
  
  random_terms <- terms %>% attr("predvars.random") %>% all.vars() %>% `[`(-1) 
  factors <- terms %>% attr("factors") 
  if (fixed_terms == "") {
      nested_terms <- rev(terms %>% attr("term.labels"))
  } else {
      nested_terms <- rev(terms %>% attr("term.labels") %>% str_subset(fixed_terms, negate = TRUE))
  }
    levels <- length(nested_terms)
    if (levels > 1) {
        for (i in 1:(levels-1)) {
            dat <- dat %>%
                bind_cols(
                    !!sym(nested_terms[i]) := dimnames(form$reTrms$Ztlist[[i]])[[1]][form$reTrms$Ztlist[[i]]@i + 1]
                )
        }
    }
  edges <- data.frame(
    from = 'Root',
    to = unique(dat[,nested_terms[levels]]),
    Height = levels + 1,
    Level = random_terms[levels]) %>%
    mutate(Name = to,
           Label = Name) %>%
    unique() #%>%
    if (levels > 1) {
        for (i in levels:2) {
            edges <- edges %>%
                bind_rows(
                    dat %>%
                    ## filter(Nest == first(Nest)) %>% 
                    {if (length(Filter)>0 & !is.na(Filter)) filter(., !!Filter[[1]])
                     else .
                    } %>% 
                    dplyr::select(from = !!nested_terms[[i]],
                                  to = !!nested_terms[[i-1]],
                                  Label = !!random_terms[[i-1]],
                                  Name = nested_terms[[i-1]],
                                  Colour = !!sym(Colour)) %>%
                    mutate(Height = i,
                           Level = random_terms[i-1]) %>% 
                    distinct() %>%
                    droplevels()
                )
        }
    }
    if (nrow(edges) != nrow(dat)) {
        edges <- edges %>%
            bind_rows(
                dat %>%
                ## filter(Nest == first(Nest)) %>%
                {if (length(Filter)>0 & !is.na(Filter)) filter(., !!Filter[[1]])
                 else .
                } %>% 
                mutate(N = 1:n(),
                       to = paste0('Rep', !!nested_terms[[1]], N),
                       Label = NA,
                       Name = to,
                       Height = 0,
                       Level = 'Reps',
                       Colour = !!sym(Colour)
                       ) %>%
                dplyr::select(from = !!nested_terms[[1]],
                              to,
                              Label,
                              Name,
                              Height,
                              Level,
                              Colour
                              )
            )
    }
  
  vertices <- edges %>%
    dplyr::select(Name, Height, Level, Label) %>%
    distinct() %>%
    rbind(data.frame(Name = 'Root', Height = 6, Level = 'Root', Label = NA)) 
  heights <- edges %>% dplyr::select(Level, Height) %>% distinct()

  graph <- graph_from_data_frame(edges, vertices = vertices)                                          
  library(igraph)
  library(ggraph)
  g <- ggraph(graph, layout = "dendrogram", height = Height) +
    ## ggraph(graph, layout = "igraph", algorithm = 'tree', circular = FALSE, height = edges$Height) +
    ## geom_edge_diagonal() +
    ## geom_edge_diagonal(aes(alpha = after_stat(index), colour = Colour), show.legend = c(alpha = FALSE, colour = TRUE)) +
    geom_edge_diagonal(aes(colour = Colour)) +
    ## geom_edge_fan(aes(alpha = after_stat(index))) +
    ## geom_edge_elbow(aes(alpha = after_stat(index))) +
    geom_node_label(aes(label = Label)) +
    ## geom_node_label(aes(label = Label, size = Height)) +
    theme_classic() 
  
  if (direction == 'vertical') {
    g <- g +
    scale_y_continuous('', breaks = heights$Height, labels = heights$Level) +
    theme(
      axis.text.y = element_text(size = rel(2)),
      panel.grid.major.y = element_line(),
      axis.line = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      ## axis.title.y = element_text(margin = margin(l = 5, unit = 'cm')),
      legend.position = 'top',
      ## legend.justification = c(1,1),
      legend.direction = 'horizontal') 
    } else {
      g <- g + scale_y_reverse('', breaks = heights$Height, labels = as.character(heights$Level)) +
        theme(
          axis.text.x = element_text(size = rel(2)),
          panel.grid.major.x = element_line(),
          axis.line = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          ## axis.title.x = element_text(margin = margin(l = 5, unit = 'cm')),
          legend.position = 'top',
          ## legend.justification = c(1,1),
          legend.direction = 'horizontal') +
        coord_flip()
    }
  
  g + guides(color = guide_none())
}
## ----end



## start by stripping out some levels
data1 <- data
data <- data1

data2 <- data |>
  dplyr::select(REEF_NAME, SITE_NO, TRANSECT_NO, FRAME, REPORT_YEAR) |>
  distinct() |>
  mutate(Y = 1) |>
  mutate(Photo = str_extract(FRAME, "\\d+$")) |>
  mutate(Photo = ifelse(REEF_NAME == first(REEF_NAME) &
    SITE_NO == first(SITE_NO) &
    TRANSECT_NO == first(TRANSECT_NO),
  Photo, NA
  )) |>
  mutate(REPORT_YEAR = factor(REPORT_YEAR)) #|> 
  #filter(!is.na(Photo))
data <- data2

Design_diagram(
  form = Y ~ REPORT_YEAR + (1 | REEF_NAME/SITE_NO/TRANSECT_NO/Photo),
  data = data2,
  Colour = "REPORT_YEAR"
)

data <- data |>
  mutate(
  filter(!(REEF_NAME != first(REEF_NAME)


  edges <- data |>
    mutate(Origin = 'Root', Height = 6, Name = SECTOR, Level = 'Region',
      forground_colour = NA, background_colour = "white") |>
    select(from = Origin, to = SECTOR, Height, Name, Level, forground_colour, background_colour) |>
    distinct() |>
    rbind(
      data |>
        mutate(Height = 5,
          Name = REEF_NAME,
          Level = 'Reef',
          forground_colour = NA,
          background_colour = "white") |>
        select(from = SECTOR,
          to = REEF_NAME,
          Height,
          Name,
          Level,
          forground_colour,
          background_colour) |>
        distinct()
    ) |> 
    rbind(
      data |>
        mutate(Height = 3,
          Name = SITE_NO,
          Site = paste(REEF_NAME, SITE_NO),
          Level = 'Site',
          forground_colour = NA,
          background_colour = "white") |>
        select(from = REEF_NAME,
          to = Site,
          Height,
          Name = Site,
          Level,
          forground_colour,
          background_colour) |>
        distinct()
    ) |> 
    rbind(
      data |>
        mutate(
          Height = 2,
          Name = TRANSECT_NO,
          Site = paste(REEF_NAME, SITE_NO),
          Transect = paste(Site, TRANSECT_NO),
          Level = "Transect",
          forground_colour = NA,
          background_colour = "white") |> 
        select(from = Site,
          to = Transect,
          Height,
          Name = Transect,
          Level,
          forground_colour,
          background_colour) |>
        distinct()
    ) ## |> 
    ## rbind(
    ##   data |>
    ##     mutate(Height = 1,
    ##       REPORT_YEAR = factor(REPORT_YEAR),
    ##       Name = REPORT_YEAR,
    ##       Level = 'Year',
    ##       forground_colour = REPORT_YEAR,
    ##       background_colour = "white") |>
    ##     select(from = TRANSECT_NO,
    ##       to = REPORT_YEAR,
    ##       Height,
    ##       Name,
    ##       Level,
    ##       forground_colour,
    ##       background_colour) |>
    ##     distinct()
    ##  ) |> 
    ## rbind(
    ##   data |>
    ##     mutate(Height = 0,
    ##       REPORT_YEAR = factor(REPORT_YEAR),
    ##       PHOTO = paste(REPORT_YEAR, str_extract(FRAME, "\\d+$")), 
    ##       Name = PHOTO,
    ##       Level = 'Photo',
    ##       ## forground_colour = REPORT_YEAR,
    ##       forground_colour = REPORT_YEAR,
    ##       background_colour = NA) |>
    ##     select(from = REPORT_YEAR,
    ##       to = PHOTO,
    ##       Height,
    ##       Name,
    ##       Level,
    ##       forground_colour,
    ##       background_colour) |>
    ##     distinct()
    ## )  

  reef_names <- data |>
    pull(REEF_NAME) |>
    unique()

  ## edges <- edges |>
  ##   mutate(Name = case_when(
  ##     Level == "Reef" ~ ifelse(Name == reef_names[1], Name, "")
  ##     ## Level == "Site" ~ ifelse(from == reef_names[1], Name, "")
  ##   ))
  ## ## edges <- edges |> filter(Name != "")
  
  vertices <- edges |>
    select(Name, Height, Level, forground_colour, background_colour) |>
    distinct() |>
    rbind(data.frame(Name = 'Root', Height = 7, Level = 'Root',
      forground_colour = NA, background_colour = "white")) |>
    mutate(
      Name2 = case_when(
        Level == 'Photo' ~ str_replace_all(Name, '.*',''),
        Level == 'Reef' ~ ifelse(Name == reef_names[1], Name, ""), #str_replace_all(Name, '.*',''),
        Level == 'Site' ~ str_replace(Name, '.*(S.*)', '\\1'),
        Level == 'Transect' ~ str_replace(Name, '.*(T.*)', '\\1'),
        Level == 'Zone' ~ str_replace(Name, '[^_]*_[^_]*_', ''),
        !Level %in% c('Site') ~ str_replace(Name, "(.*)", "\\1")),
      size = Height
    )
 
  vertices <- vertices |>
    mutate(Name2 = ifelse(Level == "Colony", NA, Name2)) |>
    mutate(size = ifelse(str_detect(Name, ".*D$") & Level == "Zone", size - 0.2, size)) |>
    mutate(size = ifelse(str_detect(Name, ".*S$") & Level == "Zone", size + 0.2, size)) |>
    mutate(background_colour = ifelse(Name2 == "", NA, background_colour))


  levels <- edges |> dplyr::select(Level, Height) |> distinct()
  graph <- igraph::graph_from_data_frame(edges, vertices = vertices)


  ggraph(graph, layout = "dendrogram", circular = FALSE, height = size, repel = TRUE) +
    geom_edge_diagonal(show.legend = FALSE) +
    geom_node_point() +
    geom_node_label(aes(label = Name2, size = Level, fill = background_colour),
      label.size = NA, repel = FALSE, show.legend = FALSE) +
    ## geom_edge_diagonal(aes(colour = edges$forground_colour), show.legend = FALSE) +
    ## geom_node_point() +
    ## geom_node_label(aes(label = Name2, size = Level, colour = forground_colour,
      ## fill = background_colour),
      ## label.size = NA, repel = FALSE, show.legend = FALSE) +
    theme_classic() +
    scale_y_continuous('', breaks = levels$Height, labels = levels$Level) +
    scale_size_manual(breaks = c('Photo','Reef', 'Region', 'Root', 'Site', 'Transect', 'Year'),
      values = c(3,5,4,4,4,3,3)) +
    scale_fill_manual(values = c("white"), na.value = "#00000000") +
    theme(
      axis.text.y = element_text(size = rel(2)),
      panel.grid.major.y = element_line(),
      axis.line = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = 'top',
      legend.direction = 'horizontal')
  


