library(here)
library(tidyverse)
library(obistools)

topp_data <- read_delim("https://docs.google.com/spreadsheets/d/e/2PACX-1vQBlYK7E8hINx9KH-0omPQXfLyrt97GMTDt7TtEw-QXH1tAl1EoSLFmiX4Ix3SRJQ/pub?gid=1045607690&single=true&output=tsv", delim = "\t")
sp_list <- read_delim("https://docs.google.com/spreadsheets/d/e/2PACX-1vQBlYK7E8hINx9KH-0omPQXfLyrt97GMTDt7TtEw-QXH1tAl1EoSLFmiX4Ix3SRJQ/pub?gid=247461962&single=true&output=tsv", delim = "\t")

# process data
occ <- topp_data %>%
  mutate(
    occurrenceStatus = "present",
    samplingProtocol = case_when(
      year == 2023 ~ "Top Predators census (page 88) https://doi.org/10.5281/zenodo.8013721",
      year == 2024 ~ "Top predators (page 80) https://doi.org/10.5281/zenodo.11653689"
    ),
    # format time into hh:mm
    time_start = format(strptime(time_start_utc_text, "%R"), "%R"),
    time_stop = format(strptime(time_stop_utc_text, "%R"), "%R"),
    # format time into hh:mmZ/hh:mmZ or only to hh:mm if time stop is NA
    eventTime = case_when(
      !is.na(time_start) & is.na(time_stop) ~ paste0(time_start, "Z"),
      !is.na(time_start) & !is.na(time_stop) ~ paste0(time_start, "Z/", time_stop, "Z"),
      TRUE ~ NA
    ),
    # ISO-8601 date yyyy-mm-dd
    eventDate = paste0(
      year, "-",
      str_pad(month, width = 2, pad = "0"), "-",
      str_pad(day, width = 2, pad = "0")
    ),
    # presume negatives for longitude and latitude
    decimalLatitude_start = round(latitude_start_degree * -1.0 + latitude_start_decimal_minute/60, 4),
    decimalLongitude_start = round(longitude_start_degree * -1.0 + longitude_start_decimal_minute/60, 4),
    decimalLatitude_end = round(latitude_stop_degree * -1.0 + latitude_stop_decimal_minute/60, 4),
    decimalLongitude_end = round(longitude_stop_degree * -1.0 + longitude_stop_decimal_minute/60, 4),
    footprintWKT = 
      case_when(!is.na(decimalLongitude_end) & !is.na(decimalLatitude_end) & decimalLongitude_start != decimalLongitude_end ~ 
                  str_c(
                    "LINESTRING (", 
                    format(decimalLongitude_start, nsmall = 4), " ", 
                    format(decimalLatitude_start, nsmall = 4), ", ", 
                    format(decimalLongitude_end, nsmall = 4), " ", 
                    format(decimalLatitude_end, nsmall = 4), ")"
                  ),
                TRUE ~ 
                  str_c(
                    "POINT (", 
                    format(decimalLongitude_start, nsmall = 4), " ", 
                    format(decimalLatitude_start, nsmall = 4), ")"
                  )
    )
  ) %>%
  mutate(centroid = calculate_centroid(footprintWKT)) %>%
  unnest_wider(centroid) %>%
  mutate(decimalLatitude = sprintf("%.4f", decimalLatitude),
         decimalLongitude = sprintf("%.4f", decimalLongitude),
         coordinateUncertaintyInMeters = case_when(
           as.integer(coordinateUncertaintyInMeters) == 0 ~ 30,
           TRUE ~ as.integer(coordinateUncertaintyInMeters)
           ))


# events

events <- occ %>% distinct(
  eventID, eventDate, year, month, day, decimalLatitude, decimalLongitude, 
  coordinateUncertaintyInMeters, footprintWKT, verbatimLocality, locality, waterBody,
  eventTime, eventRemarks, samplingProtocol, reportedWeather, time_start, time_stop,
  `Sea State (BF)`, `Cloud cover (%)`, `Vessel speed (nm)`)


events %>% distinct(year, samplingProtocol) # check protocol per expedition

# humboldt

# repeat species list for each event
humboldt <- expand_grid(events, sp_list) %>%
  select(eventID, year, targetTaxonomicScope, reportedWeather, time_start, time_stop) %>%
  mutate(
    verbatimTargetScope = "birds and marine mammals", 
    eventDurationValue = as.numeric(hm(time_stop) - hm(time_start), units = "mins"),
    eventDurationUnit = "minutes",
    taxonCompletenessReported = "notReported",
    isTaxonomicScopeFullyReported = "true",
    isAbsenceReported = "false",
    inventoryTypes = "restrictedSearch",
    protocolNames = "Top predators census",
    protocolDescriptions = case_when(
      year == 2023 ~ "Continuous monitoring of birds and marine mammals (species identification and headcount) is performed from the bridge or a spot offering the best visibility on deck. Bird/mammal standard counts are 30 min non-stop observation with binoculars for identification (if required) and age/sex determination when possible. A 600 mm tele objective camera is used for documentation and identification of species that pose identification issues in the field (e.g. Catharacta spp., Pachyptila spp.). GPS ship position and climatic conditions are recorded at each start and end position of counts. Counts are performed during daylight (from dawn to dusk), while visibility permitting (counts must be stopped when visibility is poor due to heavy fog or precipitation) to avoid bias in animal detection and subsequent false population estimates. Equipment used for the survey: Binoculars Leica Ultravid 10*32, 600 mm Long lense SONY camera (XR10iv), Garmin Oregon 600 GPS",
      year == 2024 ~ "Continuous monitoring of birds and marine mammals (species identification and headcount) is performed from the bridge or a spot offering the best visibility on deck. Bird/mammal standard counts are 30 min non-stop observation with binoculars for identification (if required) and age/sex determination when possible. A 300 mm tele objective camera is used for documentation and identification of species that pose identification issues in the field (e.g. Catharacta spp., Pachyptila spp.). GPS ship position and climatic conditions are recorded at each start and end position of counts. Counts are performed during daylight (from dawn to dusk), while visibility permitting (counts are interrupted when visibility is poor due to heavy fog or precipitation) to avoid bias in animal detection and subsequent false population estimates. Equipment used for the survey: Binoculars Leica Ultravid 10*32, 600 mm Long lense SONY camera (XR10iv), Garmin Oregon 600 GPS"
    ),
    protocolReferences = case_when(
      year == 2023 ~ "Danis, Bruno, Maria Amenabar, Annette Bombosch, Axelle Brusselman, Marius Buydens, Bruno Delille, Martin Dogniez, et al. “Report of the TANGO 1 Expedition to the West Antarctic Peninsula”. Zenodo, June 7, 2023. https://doi.org/10.5281/zenodo.8013722.",
      year == 2024 ~ "Danis, Bruno. “Report of the TANGO 2 Expedition to the West Antarctic Peninsula”. Zenodo, June 14, 2024. https://doi.org/10.5281/zenodo.11653690."
    ),
    isAbundanceReported = "true",
    isLeastSpecificTargetCategoryQuantityInclusive = "false",
    hasVouchers = "false",
    hasMaterialSamples = "false",
    samplingPerformedBy = "Henri Robert | Bruno Danis",
    isSamplingEffortReported = "true",
    samplingEffortProtocol = "Continuous monitoring of birds and marine mammals (species identification and headcount) is performed from the bridge or a spot offering the best visibility on deck. Bird/mammal standard counts are 30 min non-stop observation with binoculars for identification (if required) and age/sex determination when possible.",
    samplingEffortValue = 30,
    samplingEffortUnit = "minutes continuous monitoring"
         ) %>%
  select(eventID, targetTaxonomicScope, reportedWeather, verbatimTargetScope, 
         eventDurationValue, eventDurationUnit, taxonCompletenessReported, 
         isTaxonomicScopeFullyReported, isAbsenceReported, inventoryTypes, 
         protocolNames, protocolDescriptions, protocolReferences, isAbundanceReported,
         isLeastSpecificTargetCategoryQuantityInclusive, hasVouchers, 
         hasMaterialSamples, samplingPerformedBy, isSamplingEffortReported,
         samplingEffortProtocol, samplingEffortValue, samplingEffortUnit)

# emof

emof <- events %>%
  # cast these columns to character so that they can be pivot_longer (same type)
  mutate(across(c(`Sea State (BF)`, `Cloud cover (%)`, `Vessel speed (nm)`), as.character)) %>%
  pivot_longer(
    cols = c("Sea State (BF)", "Cloud cover (%)", "Vessel speed (nm)"),
    names_to = "measurementType",
    values_to = "measurementValue",
    values_drop_na = TRUE
  ) %>%
  mutate(
    measurementTypeID = case_when(
      measurementType == "Sea State (BF)" ~ "http://vocab.ices.dk/?ref=1705",
      measurementType == "Vessel speed (nm)" ~ "http://vocab.nerc.ac.uk/collection/P01/current/APSAGP01/"
    ),
    measurementUnit = case_when(
      measurementType == "Cloud cover (%)" ~ "Percent",
      measurementType == "Vessel speed (nm)" ~ "Knots"
    ),
    measurementUnitID = case_when(
      measurementType == "Cloud cover (%)" ~ "http://vocab.nerc.ac.uk/collection/P06/current/UPCT/",
      measurementType == "Vessel speed (nm)" ~ "http://vocab.nerc.ac.uk/collection/P06/current/UKNT/"
    )
  ) %>%
  select(
    eventID, 
    measurementType,
    measurementTypeID,
    measurementValue,
    measurementUnit,
    measurementUnitID
  )

# write files

write_tsv(events, here("data", "output", "event.txt"), na = "")
write_tsv(occ, here("data", "output", "occurrence.txt"), na = "")
write_tsv(humboldt, here("data", "output", "humboldt.txt"), na = "")
write_tsv(emof, here("data", "output", "emof.txt"), na = "")



# # plot wkt on map to check locations. Seems to only corresponds to the location if both lat lon are negatives
# library(tmap)
# occ$wkt = st_as_sfc(occ$footprintWKT)
# occ = st_sf(occ, crs = 4326)
# class(occ)
# tm_shape(occ) + tm_lines()
