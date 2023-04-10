#' Scalar potential of origin-destionation matrix by Hodge-Kodaira decompostion
#' @useDynLib HodgePotentialHumanFlow, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @keywords internal
"_PACKAGE"

.onUnload = function(libpath){
  library.dynam.unload("HodgePotentialHumanFlow", libpath)
}

## usethis namespace: start

#' This function calculates scalar potential of human flow for a given origin-destination matrix (complete (graph version).
#'
#' @param od_table A data.frame with three columns: origin, dest, and trips
#' @return a data.frame with two columns: zone and potential
#' @export
scalar_potential_on_complete_graph = function(od_table){
  unique_geozones = unique(union(od_table$origin, od_table$dest))
  res = pontential_in_complete_graph_case(od_table$origin, od_table$dest, od_table$trips, unique_geozones)
  message( sprintf("Percentage of gradient flow = %.2f",res$R2) )
  return(res$value) 
}

#' This function calculates scalar potential of human flow for a given origin-destination matrix with p-value (complete graph version).
#'
#' @param od_table A data.frame with three columns: origin, dest, and trips
#' @return This function returns a data.frame with three columns: zone, potential, p-value
#' @export
scalar_potential_on_complete_graph_with_pvalue = function(od_table){
  unique_geozones = unique(union(od_table$origin, od_table$dest))
  res = pontential_in_complete_graph_case(od_table$origin, od_table$dest, od_table$trips, unique_geozones)

  # Null model by Central limit theorem (see ref #2).
  #  P(s) \sim Normal(0, res$sd_null) 
  # Pesudo p-value, P(s > s_i) 
  res$value$pvalue = pnorm(res$value$potential, mean=0, sd = res$sd_null, lower.tail=F)
  return(res$value) 
}

#' This function calculates scalar potential of human flow for a given origin-destination matrix with p-value (unweighted graph version)
#'
#' @param flow_on_edges A data.frame with three columns: vertex1, vertex2, and netflow.
#' - num_samples: Monte calro samples (default = 1e5)
#' - seed: Monte calro samples (default = -1, then generate it by std::random_device)
#' @return This function returns a data.frame with three columns: zone, potential, p-value
#' @export
scalar_potential_on_graph = function(flow_on_edges, num_samples = 1e5, seed =-1){
  unique_geozones = unique(union(flow_on_edges$vertex1, flow_on_edges$vertex2))
  weight = rep(1, length(flow_on_edges$vertex1) )
  res = pontential_on_graph(flow_on_edges$vertex1, flow_on_edges$vertex2, weight, flow_on_edges$netflow,unique_geozones, num_samples, seed)
  return(res$value) 
}


#' This function calculates scalar potential of human flow for a given origin-destination matrix with p-value (weighted graph version)
#'
#' @param flow_on_edges A data.frame with three columns: vertex1, vertex2, weight, and netflow.
#' - num_samples: Monte calro samples (default = 1e5)
#' - seed: Monte calro samples (default = -1, then generate it by std::random_device)
#' @return This function returns a data.frame with three columns: zone, potential, p-value
#' @export
scalar_potential_on_weighted_graph = function(flow_on_edges, num_samples = 1e5, seed =-1){
  unique_geozones = unique(union(flow_on_edges$vertex1, flow_on_edges$vertex2))
  res = pontential_on_graph(flow_on_edges$vertex1, flow_on_edges$vertex2, flow_on_edges$weight,  flow_on_edges$netflow, unique_geozones, num_samples, seed)
  return(res$value) 
}


## usethis namespace: end
NULL
