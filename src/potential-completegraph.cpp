#include <Rcpp.h>
using namespace Rcpp;
#include <string>
#include <unordered_map>


//' First example
//' @param args string vector
//' @export
// [[Rcpp::export]]
List pontential_in_complete_graph_case(
    StringVector &origin, 
    StringVector &dest, 
    NumericVector &trips, 
    StringVector &unique_geozomes){

  // mapping zone string to index
  std::unordered_map<Rcpp::String, size_t> zones_to_index;
  for(size_t i = 0; i< size_t(unique_geozomes.size());++i){ zones_to_index[ unique_geozomes[i] ] = i; }

  // make a N x N netflow matrix,  Y_ij = M_ij - M_ji
  size_t N = unique_geozomes.size();
  NumericMatrix Y(N);
  Y.fill(0);
  for(size_t row = 0; row < size_t(origin.size());++row){
    size_t i = zones_to_index[ origin[row] ];
    size_t j = zones_to_index[ dest[row] ];
    Y(i,j) += trips[row];
    Y(j,i) -= trips[row];
  }

  // calc. potential
  NumericVector potential(N);
  for(size_t i = 0; i < N; ++i){
    for(size_t j = 0; j < N; ++j){
      potential[i] += Y(i,j);
    }
    potential[i] = - potential[i] / N;
  } 
  Rcpp::DataFrame df = Rcpp::DataFrame::create(Named("zone") = unique_geozomes, Named("potential") = potential);

  // calc. R2
  double R2_numerator = 0;
  double R2_denominator = 0;
  for(size_t i = 0; i < N; ++i){ 
    for(size_t j = i+1; j < N; ++j){ 
      R2_numerator += std::pow(potential[i] - potential[j], 2);
      R2_denominator += std::pow(Y(i, j) ,2);
    }
  }
  double R2= R2_numerator/ R2_denominator * 100;  // percentage
  double sum_Yij_squared = R2_denominator * 2; // upper triable

  // calc. standar deviation of normal distibution of null model
  double sd_null = sqrt(sum_Yij_squared / (N*N*N) ); 
                                                  
  List res = Rcpp::List::create(Named("value") = df, Named("R2") = R2, Named("sd_null") = sd_null);
  return(res);
}
