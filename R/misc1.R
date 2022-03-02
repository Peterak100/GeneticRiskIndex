setwd("/home/peter/GeneticRiskIndex/R")
galah_config(email = "peterkriesner@yahoo.com")

# Retrieve state counts from ALA
### use galah_identify() since v1.4 instead of select_taxa()
### use atlas_counts() instead of ala_counts()
# taxa or taxa$ala_search_term

## retrieve state counts from ALA
## this yields count for all items in 'taxa' (25212). Need to group them...
get_state_counts <- function(taxa) {
  galah_call() |>
    galah_identify(taxa$ala_search_term) |>
    galah_filter(year >= 1960, basisOfRecord == BASIS2,
                 stateProvince == STATE) |>
    atlas_counts(type = "record", limit = NULL)
}

## retrieve all counts from ALA
get_all_counts <- function(taxa) {
  galah_call() |>
    galah_identify(taxa$ala_search_term) |>
    galah_filter(year >= 1960, basisOfRecord == BASIS2) |>
    atlas_counts(type = "record", limit = NULL)
}
 

## original version of function
get_state_counts <- function(taxa) {
  ala_counts(
    taxa = select_taxa(taxa$ala_search_term), 
    filters = ALA_FILTERS,
    group_by = "species",
    type = "record",
    limit = NULL
  )
}

# "Varanus rosenbergi"
# "Tiliqua"
# cl1048 is IBRA7 Region
# cl10930 is NRM Region 2017
# cl22 is State or Territory
stumpy1 <- galah_call() |> galah_identify("Tiliqua") |>
  galah_filter(ALA_FILTERS)|>
  galah_select(basisOfRecord, cl1048, cl10930, year, group = "basic") |>
  atlas_occurrences()
write_csv(stumpy1, "/media/sf_Ubuntu_102/temp1/stumpy2.csv")

check1a <- galah_call() |>
  galah_identify("Tiliqua") |>
  galah_filter(year >= 1960, basisOfRecord == BASIS2,
               stateProvince == STATE) |>
  galah_group_by(scientificName) |>
  atlas_counts(type = "record", limit = NULL)

# need to include a call to galah_select()
check2 <- galah_call() |>
  galah_identify("Varanus rosenbergi") |>
  galah_filter(year >= 1960, basisOfRecord == BASIS2,
               stateProvince == STATE) |>
  galah_select(basisOfRecord, species, cl1048, group = "basic") |>
  atlas_occurrences()

### this works:
check3 <- galah_call() |>
  galah_identify("Varanus rosenbergi") |>
  galah_filter(stateProvince == STATE) |>
  galah_select(basisOfRecord, species, group = "basic") |>
  atlas_occurrences()
write_csv(check3, "/media/sf_Ubuntu_102/temp1/check3.csv")

# 89 possibilities for IBRA7 regions
check4 <- search_field_values("cl1048", limit = 100)


