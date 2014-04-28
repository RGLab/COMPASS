## interface matches that of COMPASS
check_meta <- function(data, counts, meta,
                       individual_id, sample_id) {

  if (!is.data.frame(meta))
    stopf("'meta' must be a 'data.frame'")

  ## ensure that individual_id, sample_id are names in the metadata
  nm <- names(meta)
  if (!(individual_id %in% nm)) {
    stopf("Expected name '%s' in the metadata (individual_id)", individual_id)
  }

  if (!(sample_id %in% nm)) {
    stopf("Expected name '%s' in the metadata (sample_id)", sample_id)
  }

  ## ensure that the names in y_s, y_u are also in the names of meta
  all_names <- names(data)
  missing_names <- all_names[ !(all_names %in% meta[[sample_id]]) ]

  if (length(missing_names)) {
    warning("There are sample ids in 'data' that are not available in ",
            "the metadata passed. These sample ids are: ",
            paste(missing_names, collapse=", "))
  }



  ## convert the sample id, individual id to character if necessary
  meta[[individual_id]] <- as.character( meta[[individual_id]] )
  meta[[sample_id]] <- as.character( meta[[sample_id]] )
  meta

}
