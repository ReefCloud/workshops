## ---- packages
library(knitr)
library(tidyverse)
library(easystats)
library(sf)
## ----end

## Read in data

### Fixed

## ---- read_data_fixed
data_fixed <- read_csv("../data/reef_data_synthetic_fixed.csv", trim_ws = TRUE)
## ----end

## ---- glimpse_data_fixed
data_fixed |> glimpse()
## ----end

## ---- head_data_fixed
data_fixed |> head()
## ----end

## ---- str_data_fixed
data_fixed |> str()
## ----end

## ---- codebook_data_fixed
data_fixed |>
  datawizard::data_codebook() |>
  knitr::kable()
## ----end

### Random

## ---- read_data_random
data_random <- read_csv("../data/reef_data_synthetic_random.csv", trim_ws = TRUE)
## ----end

## ---- glimpse_data_random
data_random |> glimpse()
## ----end

## ---- head_data_random
data_random |> head()
## ----end

## ---- str_data_random
data_random |> str()
## ----end

## ---- codebook_data_random
data_random |>
  datawizard::data_codebook() |>
  knitr::kable()
## ----end


## Excluding extraneous fields 

## ---- focus_data_fixed
data_fixed <- data_fixed |>
  dplyr::select(site_id,
    site_name,
    site_latitude,
    site_longitude,
    survey_start_date,
    survey_depth,
    survey_transect_number,
    image_id,
    image_quality,
    point_id,
    point_no,
    point_machine_classification
    )
data_fixed |> head()
## ----end

## ---- focus_data_random
data_random <- data_random |>
  dplyr::select(site_id,
    site_name,
    site_latitude,
    site_longitude,
    survey_start_date,
    survey_depth,
    survey_transect_number,
    image_id,
    image_quality,
    point_id,
    point_no,
    point_machine_classification
    )
data_random |> head()
## ----end

## Excluding images of poor quality

## ---- poor_images_data_fixed
data_fixed <-
  data_fixed |>
  dplyr::filter(is.na(image_quality) | image_quality != 0)  
## ----end

## ---- poor_images_data_random
data_random <-
  data_random |>
  dplyr::filter(is.na(image_quality) | image_quality != 0)  
## ----end

## Lengthen the data

## ---- lengthen_data_fixed
data_fixed <-
  data_fixed |>
  pivot_longer(cols = matches("point_.*_classification"),
    names_to = "type",
    values_to = "classification"
    ) 
## ----end

## ---- lengthen_data_random
data_random <-
  data_random |>
  pivot_longer(cols = matches("point_.*_classification"),
    names_to = "type",
    values_to = "classification"
    ) 
## ----end

## Joining to the label set lookup

## ---- make_labelset
labelset <- tribble(
  ~CODE, ~DESCRIPTION, ~"FUNCTIONAL GROUP", ~"KEYBOARD SHORTCUT CODE",
  "HCC", "Hard coral", "Hard coral", NA,
  "SC", "Soft coral", "Soft coral", NA,
  "MA", "Macroalgae", "Macroalgae", NA
)
write_csv(labelset, file = "../data/labelset.csv")
## ----end

## ---- read_data_labelset
labelset <- read_csv("../data/labelset.csv", trim_ws = TRUE)
## ----end

## ---- glimpse_data_labelset
labelset |> glimpse()
## ----end

## ---- head_data_labelset
labelset |> head()
## ----end

## ---- str_data_labelset
labelset |> str()
## ----end

## ---- codebook_data_labelset
labelset |>
  datawizard::data_codebook() |>
  knitr::kable()
## ----end


## ---- join_labelset_fixed
data_fixed <-
  data_fixed |>
  left_join(labelset |>
              dplyr::select(CODE, GROUP = `FUNCTIONAL GROUP`),
              by = c("classification" = "CODE")
    )
data_fixed |> as.data.frame() |> head() 
## ----end

## ---- join_labelset_random
data_random <-
  data_random |>
  left_join(labelset |>
              dplyr::select(CODE, GROUP = `FUNCTIONAL GROUP`),
              by = c("classification" = "CODE")
    )
data_random |> as.data.frame() |> head() 
## ----end



## Tally up the points

## ---- tally_data_fixed
data_fixed <- 
  data_fixed |> 
  group_by(across(c(starts_with("site"),
    starts_with("survey"),
    type,
    image_id,
    GROUP))
  ) |>
  summarise(COUNT = n(), .groups = "keep") |> 
  ungroup(GROUP) |>
  mutate(TOTAL = sum(COUNT)) |>
  ungroup() 
## ----end

## ---- tally_data_random
data_random <- 
  data_random |> 
  group_by(across(c(starts_with("site"),
    starts_with("survey"),
    type,
    image_id,
    GROUP))
  ) |>
  summarise(COUNT = n(), .groups = "keep") |> 
  ungroup(GROUP) |>
  mutate(TOTAL = sum(COUNT)) |>
  ungroup() 
## ----end
## Recode transects

## ---- recode_transects_data_fixed
data_fixed <- 
  data_fixed |>
  mutate(transect_id = paste0(site_id, survey_depth, survey_transect_number)) 
## ----end

## ---- recode_transects_data_random
data_random <- 
  data_random |>
  mutate(transect_id = paste0(site_id, survey_depth, survey_transect_number)) 
## ----end

## Fill gaps

## ---- fill_gaps_data_fixed
GROUPS <- data_fixed |> pull(GROUP) |> unique()
data.filler <- data_fixed |> 
  dplyr::select(
    starts_with("site"),
    survey_start_date,
    survey_depth,
    transect_id,
    image_id,
    type,
    TOTAL) |> 
  distinct() |> 
 tidyr::crossing(GROUP = GROUPS) 

data_fixed <-
  data_fixed |> 
  full_join(data.filler) |>
  group_by(
    across(c(starts_with("site"),
      survey_start_date,
      survey_depth,
      transect_id,
      image_id,
      type,
      GROUP
    ))) |> 
  mutate(COUNT = ifelse(is.na(COUNT), 0, COUNT),
    TOTAL = max(TOTAL, na.rm = TRUE)
  ) 
## ----end

## ---- fill_gaps_data_random
GROUPS <- data_random |> pull(GROUP) |> unique()
data.filler <- data_random |> 
  dplyr::select(
    starts_with("site"),
    survey_start_date,
    survey_depth,
    transect_id,
    image_id,
    type,
    TOTAL) |> 
  distinct() |> 
 tidyr::crossing(GROUP = GROUPS) 

data_random <-
  data_random |> 
  full_join(data.filler) |>
  group_by(
    across(c(starts_with("site"),
      survey_start_date,
      survey_depth,
      transect_id,
      image_id,
      type,
      GROUP
    ))) |> 
  mutate(COUNT = ifelse(is.na(COUNT), 0, COUNT),
    TOTAL = max(TOTAL, na.rm = TRUE)
  )
## ----end

## Sum to transect level

## ---- sum_transect_data_fixed
data_fixed <- 
  data_fixed |>
  ungroup(image_id) |>
  summarise(COUNT = sum(COUNT),
    TOTAL = sum(TOTAL)
  ) |> 
  ungroup() |> 
  droplevels()
## ----end

## ---- sum_transect_data_random
data_random <- 
  data_random |>
  ungroup(image_id) |>
  summarise(COUNT = sum(COUNT),
    TOTAL = sum(TOTAL)
  ) |> 
  ungroup() |> 
  droplevels()
## ----end


## Generate year field

## ---- generate_year_data_fixed
data_fixed <-
  data_fixed |>
  mutate(Year = lubridate::year(survey_start_date),
    TropYear = lubridate::year(survey_start_date + months(3))
  ) 
## ----end

## ---- generate_year_data_random
data_random <-
  data_random |>
  mutate(Year = lubridate::year(survey_start_date),
    TropYear = lubridate::year(survey_start_date + months(3))
  ) 
## ----end

## Generate reef field

## ---- generate_reef_data_fixed
data_fixed <-
  data_fixed |>
  mutate(Reef_id = str_replace(site_name, "(.*) Site.*", "\\1"))
## ----end

## ---- generate_reef_data_random
data_random <-
  data_random |>
  mutate(Reef_id = str_replace(site_name, "(.*) Site.*", "\\1"))
## ----end

## Visualisations


## ---- visualisation_data_fixed
data_fixed |>
  filter(type == "point_machine_classification", GROUP == "Hard coral") |> 
  ggplot(aes(y =  COUNT/TOTAL, x = survey_start_date, colour = factor(survey_depth))) +
  geom_point() +
  geom_line(aes(group = transect_id)) + 
  facet_wrap(~Reef_id + site_name)
## ----end

## ---- visualisation_data_random
data_random |>
  filter(type == "point_machine_classification", GROUP == "Hard coral") |> 
  ggplot(aes(y =  COUNT/TOTAL, x = survey_start_date, colour = factor(survey_depth))) +
  geom_point() 
  ## geom_line(aes(group = transect_id)) + 
  ## facet_wrap(~Reef_id + site_name)
## ----end
