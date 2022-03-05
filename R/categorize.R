### Prefiltering ###########################################################

# Main function called from scripts
# Categorizes risk where distance/resistance are not needed
# This has to be done in chunks as ALA database dies somewhere 
# around 1000 rows
precategorize_risk <- function(taxa) {
    n <- nrow(taxa)
    r <- rep(1:ceiling(n/GALAH_MAXROWS), each=GALAH_MAXROWS)[1:n]
    s <- lapply(split(taxa, r), precategorize_chunk)
    # Split apply combine chunks
    do.call(rbind, s)
}

# Precategorise risk for a chunk of the dataframe.
# This filters out high and low count values, regionally
# irrelevant species and data deficient species.
precategorize_chunk <- function(taxa) {
   taxa %>% 
     add_count_cols() %>%
     add_risk_col() %>%
     label_many_observations() %>%
     label_few_observations() %>%
     label_not_assessed() %>%
     label_isolation_by_resistance() %>%
     label_isolation_by_distance() %>%
     label_low_regional_relevance() %>%
     identity
}

#### ALA queries #############################################

# Add columns for each taxon of statewide and national observation counts
# creates two temporary dataframes and renames column headings
# don't need to rename 'count' column for all_counts
add_count_cols <- function(taxa) {
  cat("\nGetting taxon counts for initial filtering...\n")
  count_taxa <- filter(taxa, assess == "ALA") 
  state_counts <- get_state_counts(count_taxa) |>
    rename(state_count = count, ala_search_term = species)
  all_counts <- get_all_counts(count_taxa) |> rename(ala_search_term = species)
  taxa <- left_join(taxa, state_counts, by = "ala_search_term", all.x=TRUE)
  taxa <- left_join(taxa, all_counts, by = "ala_search_term", all.x=TRUE)
  cat("Counts retrieved successfully.\n\n")
  return(taxa)
}

# Retrieve state counts from ALA
# ala_search_term field needs to be equivalent to 'species' on ALA
# to test for a single taxon (say #6): taxa6 <- taxa[6,]
# then: x <- get_state_counts(taxa6)
# use BASIS2 because CONFIG.TOML file apparently cannot contain a list
### since v1.4
### use galah_identify() instead of select_taxa()
### use atlas_counts() instead of ala_counts()
### use galah_group_by() instead of group_by()
get_state_counts <- function(taxa) {
  galah_call() |>
    galah_identify(taxa$ala_search_term) |>
    galah_filter(year >= as.character(TIME_START),
                 year <= as.character(TIME_END),
                 basisOfRecord == BASIS2,
                 stateProvince == STATE) |>
    galah_group_by("species") |>
    atlas_counts(type = "record", limit = NULL)
}

# Retrieve national counts from ALA
# don't need to rename
get_all_counts <- function(taxa) {
  galah_call() |>
    galah_identify(taxa$ala_search_term) |>
    galah_filter(year >= as.character(TIME_START),
                 year <= as.character(TIME_END),
                 basisOfRecord == BASIS2) |>
    galah_group_by("species") |>
    atlas_counts(type = "record", limit = NULL)
}


### Filtering ########################################################

# Add a column that classifies risk. Initially "unknown".
add_risk_col <- function(taxa) {
  taxa$risk <- rep(NA, length(taxa$ala_search_term))
  return(taxa)
}

# Label very common species as "abundant"
## current thresholds: > 100K records for Vic, or > 300K records for ALA total
label_many_observations <- function(taxa) {
  ids <- taxa$state_count > MAX_OBSERVATIONS | taxa$count > MAX_OBS_TOTAL
  taxa$risk[ids] <- "abundant"
  taxa$filter_category[ids] <- "many_observations"
  return(taxa)
}

# Label very rare species as "rare"
### was previously based on taxa$state_count
label_few_observations <- function(taxa) {
  ids <- taxa$count < MIN_OBSERVATIONS
  taxa$risk[ids] <- "rare"
  taxa$filter_category[ids] <- "few_observations"
  return(taxa)
}

# Label species not relevant to STATE e.g. Victoria
label_low_regional_relevance <- function(taxa) {
  # TODO: is the proportion enough?
  ids <- (taxa$state_count / taxa$count) < MIN_PROP_IN_STATE
  taxa$risk[ids] <- "widespread"
  taxa$filter_category[ids] <- "low_proportion_in_state"
  return(taxa)
}

label_not_assessed <- function(taxa) {
  ids <- taxa$assess != "ALA"
  taxa$risk[ids] <- "not_assessed"
  taxa$filter_category[ids] <- "not_ALA_taxon"
  return(taxa)
}

## disperse_model is column 40
label_isolation_by_distance <- function(taxa) {
  taxa$filter_category[is.na(taxa$filter_category) &
                taxa$disperse_model == "Distance"] <- "isolation_by_distance"
  return(taxa)
}

label_isolation_by_resistance <- function(taxa) {
  taxa$filter_category[is.na(taxa$filter_category) &
                taxa$disperse_model == "Habitat"] <- "isolation_by_resistance"
  return(taxa)
}

# Label species for which there is not enough data.
label_data_deficient <- function(taxa) {
  # TODO: add something here
  # not sure what data deficient means in practice
  return(taxa)
}
