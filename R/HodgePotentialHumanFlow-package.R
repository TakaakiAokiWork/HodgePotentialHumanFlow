#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' Scalar potential of origin-destionation matrix by Hodge-Kodaira decompostion
#'
#' This function calculates scalar potential of human flow for a given origin-destionation matrix.
#' @param od_table A data.frame with three columns: origin, dest, and trips
#' @return a data.frame with two columns: geozone and negative potential
#' @export
scalar_potential = function(od_table){
  unique_geozones = unique(union(od_table$origin, od_table$dest))
  s = numeric(length(unique_geozones))
  names(s) = unique_geozones
  for(i in 1:nrow(od_table)){
    o = od_table[i,1]
    d = od_table[i,2]
    trips = od_table[i,3]
    s[o] = s[o] - trips
    s[d] = s[d] + trips
  }
  for(i in 1:length(s)){
    s[i] = s[i] / length(s)
  }
  df = stack(s)
  names(df) = c("NegativePotential", "geozone")
  return(df) 
}
## usethis namespace: end
NULL
