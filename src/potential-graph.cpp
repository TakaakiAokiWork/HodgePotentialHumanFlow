#include <RcppGSL.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_linalg.h>

using namespace Rcpp;
#include <string>
#include <unordered_map>
#include <random>
#include <algorithm>
#include <limits>

#include <omp.h>

bool is_equal(double a,double b){
  if(abs(a-b) < std::numeric_limits<double>::epsilon()){ return true; }
  return false;
}

//' First example
//' @param args string vector
//' @keywords internal
// [[Rcpp::export]]
List pontential_on_graph(
    StringVector &vertex1,  // (v1,v2,weight) => weighted edge list (default weight =1)
    StringVector &vertex2, 
    NumericVector &weight, 
    NumericVector &netflow_R, // Y_ij edge flow for v1,v2 
    StringVector &unique_geozomes,
    size_t num_samples,
    unsigned long int seed){

  // mapping zone (string) to index (size_t)
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
    std::vector<double> degree(N,0);
    for(size_t m = 0; m < edgelist.size(); ++m){
      int i = edgelist[m].first;
      int j = edgelist[m].second;
      double w = weight[m];
      gsl_matrix_set(L,i,j,-w);
      gsl_matrix_set(L,j,i,-w);
      degree[i] +=w;
      degree[j] +=w;
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
  gsl_vector * potential = gsl_vector_calloc(N);
  calculate_potential(netflow, potential); 


  /*
   *   Calculate potential distribution of the null model (Monte Calro) if num_samples > 0
   */
  std::vector<double> pvalue_above(N,0.0); // prob. when  s_i(null) > s_i(original)
  if (num_samples > 0){
    {
      Rprintf("Calculating p-values by monte calro (%ld steps).", num_samples);
      std::bernoulli_distribution dist_flip(0.5);

      // MC sampling
      std::vector<size_t> num_exceed_samples_in_nullmodel(N,0); // #samples if s_i(null) > s_i(original)

#pragma omp parallel
      {
#pragma omp single
        {
          Rprintf(" thread = %d", omp_get_num_threads() );
        }
        // shared local variables
        gsl_vector * tmp_potential = gsl_vector_calloc(N);
        std::vector<double> randomized_flow( netflow );
        std::vector<size_t> num_exceed_local(N,0);
        std::mt19937_64 rng;

#pragma omp critical
        {
          // random generator for each thread
          std::random_device rd;
          std::seed_seq seq{ rd(), rd(), rd(), rd()};
          rng.seed(seq);
        }

#pragma omp for schedule(static) nowait
        for(size_t n = 0; n < num_samples; ++n){
          // random order
          std::shuffle(randomized_flow.begin(), randomized_flow.end(), rng);
          // flip by prob. of 0.5
          for(size_t m = 0; m < randomized_flow.size(); ++m){
            if ( dist_flip(rng)){ randomized_flow[m] = - randomized_flow[m]; }
          }
          calculate_potential(randomized_flow, tmp_potential);
          for(size_t i =0; i < N; ++i){
            if (tmp_potential->data[i] > potential->data[i]){
              num_exceed_local[i] += 1;
            }
          }
        }
#pragma omp critical
        {
          for(size_t i =0; i < N; ++i){ num_exceed_samples_in_nullmodel[i] += num_exceed_local[i];} // reduction
        }

        gsl_vector_free(tmp_potential);
      } // end of parallel

      for(size_t i=0; i< N; ++i){
        pvalue_above[i] = double(num_exceed_samples_in_nullmodel[i]) / num_samples;
      }
    }
  }


  // calc. R2
  double R2_numerator = 0;
  double R2_denominator = 0;
  for(size_t m = 0; m < edgelist.size(); ++m){
    int i = edgelist[m].first;
    int j = edgelist[m].second;
    double Yij = netflow[m];
    R2_numerator += std::pow(potential->data[j] - potential->data[i], 2);
    R2_denominator += std::pow(Yij ,2);
  }
  double R2= R2_numerator/ R2_denominator * 100;  // percentage

  // make a return object
  NumericVector potential_R(N);
  for(size_t i = 0; i < N; ++i){ potential_R[i] = potential->data[i]; } 

  Rcpp::DataFrame df;
  if (num_samples >0){
    df = Rcpp::DataFrame::create(Named("zone") = unique_geozomes, Named("potential") = potential_R, Named("pvalue_above") = wrap(pvalue_above) );
  } else{
    df = Rcpp::DataFrame::create(Named("zone") = unique_geozomes, Named("potential") = potential_R );
  }

  List res = Rcpp::List::create(Named("value") = df, Named("R2") = R2, Named("seed_used") = seed);
  return(res);
}
