setwd("/home/peter/GeneticRiskIndex/R")
galah_config(email = "peterkriesner@yahoo.com")

# Retrieve state counts from ALA
### use galah_identify() since v1.4 instead of select_taxa()
### use atlas_counts() instead of ala_counts()
get_state_counts <- function(taxa) {
  atlas_counts(
    taxa = galah_identify(taxa$ala_search_term), 
    galah_filter(year >= 1960, basisOfRecord == BASIS),
    type = "record",
    limit = NULL
  )
}

# "Varanus rosenbergi"
# "Tiliqua"
# cl1048 is IBRA7 Region
# cl10930 is NRM Region 2017
# cl22 is State or Territory
stumpy1 <- galah_call() |> galah_identify("Varanus rosenbergi") |>
  galah_filter(year >= 1960, year <= 2021, basisOfRecord == "HUMAN_OBSERVATION")|>
  galah_select(basisOfRecord, cl1048, cl10930, year, group = "basic") |>
  atlas_occurrences()
write_csv(stumpy1, "/media/sf_Ubuntu_102/temp1/stumpy2.csv")

