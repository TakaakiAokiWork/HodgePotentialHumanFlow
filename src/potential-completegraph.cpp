#include <Rcpp.h>
using namespace Rcpp;
#include <string>
#include <unordered_map>


//' First example
//' @param args string vector
//' @export
// [[Rcpp::export]]
Rcpp::List pontential_in_complete_graph_case(
    Rcpp::StringVector &origin, 
    Rcpp::StringVector &dest, 
    Rcpp::NumericVector &trips, 
    Rcpp::StringVector &unique_geozomes){

  // mapping zone string to index
  std::unordered_map<Rcpp::String, size_t> zones_to_index;
  for(size_t i = 0; i< size_t(unique_geozomes.size());++i){ zones_to_index[ unique_geozomes[i] ] = i; }

  Rcpp::NumericVector potential( unique_geozomes.size() );

  size_t nrow = origin.size();
  for(size_t i = 0; i < nrow;++i){
    potential[ zones_to_index[ origin[i] ] ] -= trips[i];
    potential[ zones_to_index[ dest[i]   ] ] += trips[i];
  }
  for(auto &x : potential){
    x /= unique_geozomes.size();
  }
 

  Rcpp::DataFrame df = Rcpp::DataFrame::create(Named("zone") = unique_geozomes, Named("potential") = potential);
  Rcpp::List res = Rcpp::List::create(Named("value") = df);
  return(res);
}
