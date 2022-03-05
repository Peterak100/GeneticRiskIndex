source("common.R")
source("categorize.R")
source("observations.R")
source("resistance.R")
source("distance.R")
source("fire_severity.R")

# Load main taxa dataframe from csv ##########################################
taxa <- read.csv(BATCH_TAXA_CSV_PATH, header = TRUE)
# head(taxa)

##############################################################################
# Precategorize based on counts
# Categorize risk using queries to ALA: slow.
## adds 4 columns to BATCH_TAXA_CSV:
## "state_count", "count", "risk", "filter_category"
precategorized_taxa <- precategorize_risk(taxa)
# head(precategorized_taxa[,c(1,2,99:102)],14)
# precategorized_taxa$filter_category

# Taxa to access based on distance metrics
# this is saved to the /data/groupings folder
write_csv(precategorized_taxa,
          file.path(groupingspath, "precategorized_taxa.csv"))


#############################################################################
## Manual single taxon observations and preclustering for testing:
taxon <- precategorized_taxa[1, ]
taxon$ala_search_term
## these 3 functions are from observations.R
## filter_observations(taxon) has 3 (or 4) sub functions...
## the last of these is from fire_severity.R
## precluster_observations(taxon) has 2 sub functions...

obs <- load_or_download_obs(taxon, taxapath, force_download=FALSE) |>
  filter_observations(taxon) |>
  precluster_observations(taxon)
# check preclustering...
head(obs[,7:9],43) |> as.data.frame()

taxonpath <- taxon_path(taxon, taxapath)
# combines "x" and "y" columns into single "geometry" column
shapes <- sf::st_as_sf(obs, coords = c("x", "y"), crs = METRIC_EPSG)
### c.f. observations.R line 234,
scaled_eps <- taxon$epsilon * 1000 * EPSILON_SENSITIVITY_SCALAR

## Create a full-sized raster for preclusters
preclustered <- buffer_preclustered(shapes, scaled_eps)
cat("Preclusters:", nrow(preclustered), "\n")
## this can now be plotted to check result...
plot(preclustered)

precluster_rast <- shapes_to_raster(preclustered, taxon, mask_layer, taxonpath)
## this can also be plotted
### COULD NOT FIND FUNCTION shape_to_raster (TYPO?)


# pixel_freq <- freq(precluster_rast)
# pixel_freq
# colnames(shapes)
# left_join(shapes, pixel_freq, copy=TRUE, by=c("precluster" = "value")) %>%
#   write_csv(file.path(taxapath, "preclusters.csv"))
#
# process_observations(taxon, mask_layer, taxapath, error=TRUE)


##############################################################################
# Clustering for Isolation by distance and resistance taxa

isolation_taxa <- filter(precategorized_taxa, filter_category %in%
            c("isolation_by_distance", "isolation_by_resistance"))
nrow(isolation_taxa)
head(isolation_taxa)

## 4th March 5pm checked above here!!!!!!!

# load/download, filter and precluster observations for all taxa
preclustered_isolation_taxa <- process_observations(isolation_taxa,
            mask_layer, taxapath, throw_errors=TRUE)
# FIXME: there are lots of errors here, something to do with SF coercion. 
# Use throw_errors=TRUE to use traceback()
preclustered_isolation_taxa$error
head(preclustered_isolation_taxa)
nrow(preclustered_isolation_taxa)


##############################################################################
# Main CSV output

# Write csv for all categorized taxa

# FIXME: this duplicates columns giving them x/y suffixes.
# The idea is to fill in precategorized_taxa with
# any new columns from preclustered_isolation_taxa, and use values from
# preclustered_isolation_taxa for other columns, where they are different.

# get the added columns form preclustered_isolation_taxa (*order matters here*)
diff_cols <- setdiff(names(preclustered_isolation_taxa),
              names(precategorized_taxa))
# first update the precategorized_taxa with matching columns from
# preclustered_isolation_taxa then left_join with the remaining columns
# in preclustered_isolation_taxa
categorized_taxa <- rows_update(precategorized_taxa, 
            preclustered_isolation_taxa[, names(precategorized_taxa)],
            by = "ala_search_term") %>% 
  left_join(preclustered_isolation_taxa[, c("ala_search_term", diff_cols)],
            by = "ala_search_term")

write_csv(categorized_taxa, file.path(groupingspath, "catagorized_taxa.csv"))
head(categorized_taxa)
categorized_taxa$error


##############################################################################
# Circuitscape/isolation by resistance output

# Write csv for taxa that we need to process with Circuitscape
isolation_by_resistance_taxa <- filter(preclustered_isolation_taxa,
            is.na(risk), filter_category == "isolation_by_resistance")
head(isolation_by_resistance_taxa)
nrow(isolation_by_resistance_taxa)
write_csv(isolation_by_resistance_taxa, file.path(groupingspath,
            "isolation_by_resistance_taxa.csv"))

# Write as a single column job-list
job_file <- file(file.path(datapath, "batch_jobs.txt"))
underscored <- gsub(" ", "_", isolation_by_resistance_taxa$ala_search_term)
writeLines(underscored, job_file)
close(job_file)

# Download and write raster files for Circuitscape resistance models
prepare_resistance_files(isolation_by_resistance_taxa, taxapath)
