#' Scalar potential of origin-destionation matrix by Hodge-Kodaira decompostion
#' @useDynLib HodgePotentialHumanFlow, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @keywords internal
"_PACKAGE"

.onUnload = function(libpath){
  library.dynam.unload("HodgePotentialHumanFlow", libpath)
}

## usethis namespace: start
#' This function calculates scalar potential of human flow for a given origin-destionation matrix.
#'
#' @param od_table A data.frame with three columns: origin, dest, and trips
#' @return a data.frame with two columns: geozone and negative potential
#' @export
scalar_potential = function(od_table){
  unique_geozones = unique(union(od_table$origin, od_table$dest))
  res = pontential_in_complete_graph_case(od_table$origin, od_table$dest, od_table$trips, unique_geozones)
  
#  s = numeric(length(unique_geozones))
#  names(s) = unique_geozones
#  for(i in 1:nrow(od_table)){
#    o = od_table[i,1]
#    d = od_table[i,2]
#    trips = od_table[i,3]
#    s[o] = s[o] - trips
#    s[d] = s[d] + trips
#  }
#  for(i in 1:length(s)){
#    s[i] = s[i] / length(s)
#  }
#  df = stack(s)
  #names(df) = c("NegativePotential", "geozone")
  return(res$value) 
}
## usethis namespace: end
NULL
