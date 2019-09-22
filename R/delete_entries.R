# Private functions to subset the data
delete_entries = function(x, ids)
{
  x$mutations = x$mutations %>%
    filter(!id %in% ids)

  x$mutations_locations = x$mutations_locations %>%
    filter(!id %in% ids)

  x$mutations_annotations = x$mutations_annotations %>%
    filter(!id %in% ids)

  if(has_mobster_fits(x)) {
    warning("There are MOBSTER clusters, deleting mutations requires to re-compute the clusters")
  }

  if(has_viber_fits(x)) {
    warning("There are VIBER clusters, deleting mutations requires to re-compute the clusters")
  }

  x
}
