#include <RcppGSL.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_linalg.h>

using namespace Rcpp;
#include <string>
#include <unordered_map>
#include <limits>


bool is_equal(double a,double b){
  if(abs(a-b) < std::numeric_limits<double>::epsilon()){ return true; }
  return false;
}

//' First example
//' @param args string vector
//' @export
// [[Rcpp::export]]
List pontential_on_graph(
    StringVector &vertex1, 
    StringVector &vertex2, 
    NumericVector &netflow_R, 
    StringVector &unique_geozomes){

  // mapping zone string to index
  std::unordered_map<Rcpp::String, size_t> zones_to_index;
  for(size_t i = 0; i< size_t(unique_geozomes.size());++i){ zones_to_index[ unique_geozomes[i] ] = i; }

  // store given edges and netflow
  size_t N = unique_geozomes.size();
  std::vector<double> netflow = as< std::vector<double> >(netflow_R);
  size_t M = netflow.size();
  using Edge = std::pair<size_t, size_t>;
  std::vector<Edge> edgelist;
  for(size_t m = 0; m < M; ++m){
    int i = zones_to_index[ vertex1[m] ];
    int j = zones_to_index[ vertex2[m] ];
    edgelist.emplace_back(i,j);
  }

  // set up Laplacian matrix
  gsl_matrix *L = gsl_matrix_calloc(N,N);
  {
    std::vector<int> degree(N,0);
    for(size_t m = 0; m < edgelist.size(); ++m){
      int i = edgelist[m].first;
      int j = edgelist[m].second;
      gsl_matrix_set(L,i,j,-1);
      gsl_matrix_set(L,j,i,-1);
      degree[i] +=1;
      degree[j] +=1;
    }
    for(size_t i = 0; i < N; ++i){
      gsl_matrix_set(L, i, i, degree[i]);
    }
    // validation
    for(size_t i = 0; i < N; ++i){ 
      for(size_t j = 0; j < N; ++j){ 
        double Lij = gsl_matrix_get(L, i,j);
        double Lji = gsl_matrix_get(L, j,i);
        assert( is_equal( Lij, Lji ) );
      }
    }
  }

  // Singular Value Decomposition to solve Ax=b.
  gsl_matrix * mat_V = gsl_matrix_calloc(N,N);
  gsl_vector * singular = gsl_vector_calloc(N); 
  gsl_vector * work = gsl_vector_calloc(N); 
  gsl_linalg_SV_decomp(L, mat_V, singular, work); // L is replaced by U

  auto calculate_potential = [&L, &mat_V, &singular, &N, &edgelist](const std::vector<double> & Y, gsl_vector * s){
    gsl_vector * minus_divY = gsl_vector_calloc(N); // = - div Y
    for(size_t m = 0; m < edgelist.size(); ++m){
      int i = edgelist[m].first;
      int j = edgelist[m].second;
      double Yij = Y[m];
      minus_divY->data[i] += - Yij;
      minus_divY->data[j] += Yij;
    }
    gsl_linalg_SV_solve(L, mat_V, singular, minus_divY, s);
    // shift s to be mean <s> = 0
    double mean =  gsl_stats_mean(s->data, s->stride, s->size);
    gsl_vector_add_constant(s, - mean);
    gsl_vector_free(minus_divY);
  };

  // calc. potential
  gsl_vector * s = gsl_vector_calloc(N);
  calculate_potential(netflow, s); 

  NumericVector potential(N);
  for(size_t i = 0; i < N; ++i){ potential[i] = s->data[i]; } 
  Rcpp::DataFrame df = Rcpp::DataFrame::create(Named("zone") = unique_geozomes, Named("potential") = potential);

  // calc. R2
  double R2_numerator = 0;
  double R2_denominator = 0;
  for(size_t m = 0; m < edgelist.size(); ++m){
    int i = edgelist[m].first;
    int j = edgelist[m].second;
    double Yij = netflow[m];
    R2_numerator += std::pow(potential[i] - potential[j], 2);
    R2_denominator += std::pow(Yij ,2);
  }
  double R2= R2_numerator/ R2_denominator * 100;  // percentage

  List res = Rcpp::List::create(Named("value") = df, Named("R2") = R2);
  return(res);
}
