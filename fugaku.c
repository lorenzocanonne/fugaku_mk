#define _GNU_SOURCE
#include <bits/pthreadtypes-arch.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <stdbool.h>
#include <string.h>
#include <mpi.h>


#define MAX_COLOR_SMALL_OFFSET 64
#define MAX_COLOR_OFFSET 1024


//{ RANDOM GENERATORS

  //pick up Park_miller random numbers gp/solexa_3623/source gpquick_seed
  int gpquick_seed;
  extern int gpquick_seed;
  //#include "pch.h"
  // pch.h   Adapted from GPQUICK-2.1
  // W. Langdon cs.ucl.ac.uk 5 May 1994 (21 Feb 95 make it )
  // Use Park-Miller random numbers

  // QPQUICK
  // Standard header files for compilers supporting Pre-Compiled Headers

  //WBL 31 Mar 18  for gcc 4.8.5 NB C not C++
  //WBL 19 Oct 03  Port to Alpha gcc 3.3
  //WBL 21 Dec 02  Port to Microsoft Visual C++

  #include <limits.h>

  //////////  Random number functions cause portability problems.  Resolve here.
  // Define	rand_0to1()  (random int from 0 to 1)
  //			rnd(x) (random integer from 0 to x-1)
  //			rndomize() (initialize with a time based seed) 			
  // park-miller.cc   Park-Miller random number generator
  // W. Langdon cs.ucl.ac.uk 5 May 1994

  #ifndef LONG_MAX
  error LONG_MAX missing
  #endif

  // 1<=seed<=m
   int intrnd (int* seed) {
    #if LONG_MAX > (16807*2147483647)
      long int const seed_= *seed;
      long int const a    = 16807;      //ie 7**5
      int const m    = 2147483647; //ie 2**31-1
      *seed = (seed_ * a)%m;
      return *seed;
    #else
      int const a    = 16807;      //ie 7**5
      int const m    = 2147483647; //ie 2**31-1
      int temp = *seed * a;
      *seed = (int) (temp - m * floor ( temp / m ));
      return *seed;
    #endif
  }

   double drand_0to1()
          { return (double)(intrnd(&gpquick_seed)) / 2147483647.0; }

   int thread_safe_rnd(int * rnd_gen, int __num){ 
    return (intrnd(rnd_gen) % __num); 
  }

   int rnd(int __num){ 
    return (intrnd(&gpquick_seed) % __num); 
  }

   void rndomize(void){ 
    do {
      gpquick_seed = (unsigned) time(NULL);
    } while (gpquick_seed <= 0);
  }
  //end pch.h

   int myrnd(const int __num, int* prng_NK2) {
    const int save = gpquick_seed;
    gpquick_seed = *prng_NK2;
    const int ans = rnd(__num);
    *prng_NK2 = gpquick_seed;
    gpquick_seed = save;
    return ans;
  }

   double myrand_0to1(int* prng_NK1) {
    const int save = gpquick_seed;
    gpquick_seed = *prng_NK1;
    const double ans = drand_0to1();
    *prng_NK1 = gpquick_seed;
    gpquick_seed = save;
    return ans;
  }

   void thread_rand_mutation(const int N, char solution[N], int rand_array[N], const int rand_number_flips, int * rnd_gen){
    for(int v=0;v<rand_number_flips;v++){
      const int idmax = N-v-1;
      const int r = thread_safe_rnd(rnd_gen,idmax+1);
      const int tmp = rand_array[r];
      rand_array[r] = rand_array[idmax];
      rand_array[idmax] = tmp;
      solution[tmp] = (~solution[tmp]) & 1;
    }
  }
//}

//{ LANDSCAPE CREATION/PRINT/FREE
  int const_fit_stride;
  double* f_values;
  int* f_l;
  int* v_i;
  int* v_i_head;
  int* v_i_degree;
  int* vig_degree;
  int* vig_head;
  int* vig;

  void mk_landscape(const int M, const int N, const int K, const int Q, const int seed_functions, const int seed_values){
    const_fit_stride            = 1<<K;
    f_values                    = malloc(const_fit_stride*M*sizeof(double));
    f_l                         = malloc(M*K*sizeof(int));
    v_i_degree                  = calloc(N,sizeof(int));
    v_i_head                    = calloc(N,sizeof(int));
    v_i                         = malloc(M*K*sizeof(int));
    const int const_v_degree_upper_bound = 10*K*K;
    int * v_i_tmp = malloc(const_v_degree_upper_bound*N*sizeof(int));
    int seed1 = seed_values;
    int seed2 = seed_functions;
    int i,j;
    for (i=0; i<M; i++) {
      const int l = const_fit_stride*i;
      for (j=0; j<const_fit_stride; j++) {
        f_values[l+j] = myrand_0to1(&seed1);
      }
    } 
    for (i=0; i<M; i++) {
      const int l = i*K;
      for (j=0; j<K; j++) {
        const int idx = l+j;
        const int v = f_l[idx] = (i+1+myrnd(N-1,&seed2))%N;
        bool v_already_in_f = false;
        for (int k=0;k<j;k++) {
          if (v==f_l[l+k]) {
            v_already_in_f = true;
            break;
          }
        }
        if(!v_already_in_f){
          v_i_tmp[v*const_v_degree_upper_bound+(v_i_degree[v]++)]=i;
        }
      }
    }
    int tank=0;
    for (int v=0; v<N; v++) {
      v_i_head[v]=tank;
      const int start_idx = v*const_v_degree_upper_bound;
      for (int i=0; i<v_i_degree[v]; i++) {
        v_i[tank++] = v_i_tmp[start_idx+i]; 
      }
    }
    free(v_i_tmp);
  }

  void print_landscape(const int M, const int K){
    printf("Functions\n");
    for(int i=0;i<M;i++) {
      const int l = i*K;
      printf("f(%i) {%i",i,f_l[l]);
      for(int j=1;j<K;j++) {
        printf(",%i",f_l[l+j]);
      }
      printf("}\n");
    }
  }

  void print_vig(const int N){
    if (vig == NULL) {
      printf("VIG not defined\n");
      return;
    }
    else {
      printf("VIG\n");
      for(int v=0;v<N;v++) {
        int idx = vig_head[v];
        printf("v(%i) {",v);
        for(int j=0;j<vig_degree[v];j++) {
          printf(",%i",vig[idx++]);
        }
        printf("}\n");
      }
    }
  }

  void free_vig(){
    if (vig_head != NULL) {
      free(vig_head);
    }
    if (vig_degree != NULL) {
      free(vig_degree);
    }
    if (vig != NULL) {
      free(vig);
    }
  }

  void free_landscape(){
    if (f_values != NULL) {
      free(f_values);
    }
    if (f_l != NULL) {
      free(f_l);
    }
    if (v_i != NULL) {
      free(v_i);
    }
    if (v_i_degree != NULL) {
      free(v_i_degree);
    }
    if (v_i_head != NULL) {
      free(v_i_head);
    }
    free_vig();
  }
  
  void sequential_vig(const int M, const int N, const int K){
    vig_degree                  = calloc(N,sizeof(int));
    vig_head                    = calloc(N,sizeof(int));
    vig                         = malloc(M*K*K*sizeof(int));
    int tank=0;
    for (int v=0; v<N; v++) {
      vig_head[v]=tank;
      const int start_idx_l = v_i_head[v];
      for (int i=0; i<v_i_degree[v]; i++) {
        const int l = v_i[start_idx_l+i];
        for (int k=0; k<K; k++) {
          const int w = f_l[l*K+k];
          if(w == v){
            continue;
          }
          bool found = false;
          for (int j=vig_head[v]; j<tank; j++) {
            if (w == vig[j]) {
              found = true;
              break;
            }
          }
          if (!found) {
            vig[tank++]=w;
          }
        }
      }
      vig_degree[v]=tank-vig_head[v];
    }
  }

  

//}

//{ SCORE/CONTRIBUTION
   double subfunction_contribution(const int l, const int N, const int K, const char solution[N]){
    unsigned int L = 0;
    for(int k=0;k<K;k++) {
      const int w = f_l[l*K+k];
      const unsigned int s_v = solution[w];
      L = L + (s_v<<(k+1));
    }
    return f_values[const_fit_stride*l+L];
  }

   double compute_score(const int M, const int N, const int K, const char*solution){
    double sum = 0.0;
    for (int l=0; l<M; l++) {
      sum += subfunction_contribution(l, N, K, solution);
    }
    return sum;
  }
//}

//{ UTILITIES
   void shuffle(const int SIZE,int rand_array[SIZE], int * rnd_gen){
    for (int i=0; i<SIZE-1; i++) {
      int upper_bound = SIZE-i;
      const int r = thread_safe_rnd(rnd_gen,upper_bound);
      const int tmp = rand_array[r];
      rand_array[r] = rand_array[upper_bound-1];
      rand_array[upper_bound-1] = tmp;
    }
  }

   void activate_adjacent_variables(const int v, char * states){
    for (int i=vig_head[v]; i<vig_head[v]+vig_degree[v]; i++) {
      states[vig[i]] = 0;
    }
  }
//}

//{ SOLUTION
   int is_optimal(const int N, const int K, const char solution[N]){
    int cpt =0;
    for (int v=0; v<N; v++) {
      if (variable_score_flip(v,N,K,solution) > 0.0) {
        cpt++;
      }
    }
    return cpt;
  }

   void generate_random_solution(const int N, char solution[N], const int seed_solution){
    int seed = seed_solution;
    for (int v=0; v<N; v++) {
    	solution[v] = thread_safe_rand(&seed,2);
    }
  }
//}

//{ COLORING

  int * color_ascendant;
  int * number_colors;
  int * thread_works_bound;
  int * rand_color;

   int adjacent_colors_vig(const int v, int *C, int *v_color){
    int iterator_C = 0;
    const int start_idx = vig_head[v];
    for (int i=0; i<vig_degree[v]; i++) {
      const int w = vig[start_idx+i];
      bool found = true;
      for (int c=0; c<iterator_C; c++) {
        if (C[c] == v_color[w]) {
          found = false;
          break;
        }
      }
      if (found) {
        C[iterator_C++]=v_color[w];
      }
    }
    return iterator_C;
  }

   int get_color(const int v, int*C, int *v_color){
    const int iterator_C = adjacent_colors_vig(v, C, v_color);
    int color_min = 0;
    bool found = true;
    while (found) {
      found = false;
      for (int c = 0; c<iterator_C; c++) {
        if (C[c] == color_min) {
          found = true;
          color_min++;
          break;
        }
      }
    }
    return color_min;
  }

  int compare_colors(const void* a, const void* b, void *v_color){
    int * p = v_color;
    return (p[*(int*)a] - p[*(int*)b]);
  }

   void coloring_ordering(const int N, void *v_color, int * color_ascendant){
    for(int i=0;i<N;i++){
      color_ascendant[i]=i;
    }
    qsort_r(color_ascendant, N, sizeof(int), compare_colors, v_color);
  }

  void colors_bound_computation(const int N, const int t, const int nb_colors, int* colors_bound, int color_ascendant[N], int v_color[N]){
    colors_bound[0]=0;
    int j=0;
    for (int i=0; i<N; i++) {
      if (v_color[color_ascendant[i]] != j) {
        j++;
        colors_bound[j]=i;
      }
    }
    colors_bound[nb_colors] = N;
  }

  void dispatch_among_threads(const int N, const int T, int thread_works_bound[MAX_COLOR_OFFSET], const int nb_colors, int colors_bound[nb_colors+1]){
    int it = 0;
    thread_works_bound[0] = 0;
    for (int color=0; color<nb_colors; color++) {
      const int start_idx         = colors_bound[color];
      const int end_idx           = colors_bound[color+1];
      const int number_works      = end_idx-start_idx;
      const int works_per_thread  = number_works/T;
      const int carry             = number_works-(works_per_thread*T);
      for(int t = 0 ; t < T ; t++){
        const int idx       = (color*T) + t +1 ;
        if(t<carry){    
          thread_works_bound[idx] = it+= works_per_thread+1;
        }
        else {
          thread_works_bound[idx] = it+= works_per_thread;
        }
      }
    }
    thread_works_bound[T*nb_colors] = N;
  }

  void parallel_coloring(const int N, const int T, int * seed){
    int * seeds = malloc(T*sizeof(int));
    for(int t=0;t<T;t++){
      seeds[t] = thread_safe_rnd(seed, N);
    }
    color_ascendant = malloc(N*T*sizeof(int));
    thread_works_bound = malloc(MAX_COLOR_OFFSET*T*sizeof(int));
    number_colors = malloc(T*sizeof(int));
    rand_color = malloc(MAX_COLOR_OFFSET*T*sizeof(int));
    #pragma omp parallel num_threads(T)
    {
      const int t = omp_get_thread_num();
      int * v_color = calloc(N,sizeof(int));
      int * tmp_rand = malloc(N*sizeof(int));
      for(int i=0;i<N;i++){
        tmp_rand[i] = i;
      }
      shuffle(N,tmp_rand,&seeds[t]);
      int * tmp_colors_buffer = malloc(1000*sizeof(int));
      for (int i=0;i<N; i++) {
        const int v = tmp_rand[i];
        v_color[v] = get_color(v, tmp_colors_buffer, v_color);
      }
      free(tmp_colors_buffer);
      free(tmp_rand);
      int *colors_bound;
      coloring_ordering(N,v_color,&color_ascendant[t*N]);
      const int nb_colors = number_colors[t] = v_color[color_ascendant[(t+1)*N-1]]+1;
      int * p = &rand_color[t*MAX_COLOR_SMALL_OFFSET];
      for(int i=0;i<nb_colors;i++){
        p[i]=i;
      }
      colors_bound = malloc((nb_colors+1)*sizeof(int));
      colors_bound_computation(N,t,nb_colors,colors_bound,&color_ascendant[t*N],v_color);
      dispatch_among_threads(N,T,&thread_works_bound[MAX_COLOR_OFFSET*t],nb_colors,colors_bound);
      free(v_color);
      free(colors_bound);
    }
    free(seeds);
  }

  void free_coloration(){
    if (number_colors != NULL) {
      free(number_colors);
    }
    if (color_ascendant != NULL) {
      free(color_ascendant);
    }
    if (thread_works_bound != NULL) {
      free(thread_works_bound);
    }
    if (rand_color != NULL) {
      free(rand_color);
    }
  }

  void print_nb_colors_multiple_coloration(const int T){
    for (int t=0; t<T; t++) {
      printf("nb color %i => %i\n",t,number_colors[t]);
    }
  }

  void print_color_bounds(const int T){
    for (int t=0; t<T; t++) {
      printf("nb color %i => %i\n",t,number_colors[t]);
      for(int c = 0;c<number_colors[t];c++){
        printf("    color %i \n",c);
        for(int i=0;i<T;i++){
          printf("        t%i => {%i; %i}\n",i,thread_works_bound[MAX_COLOR_OFFSET*t+c*T+i],thread_works_bound[MAX_COLOR_OFFSET*t+c*T+i+1]);
        }
      }
    }
  }

//}

//{ PARALLEL HILL CLIMBER
  
  double init_parallel_hill_climber(const int M, const int N, const int K, const int T, const char solution[N], double scores[N]){
    double score = 0.0;
    #pragma omp parallel for num_threads(T) reduction(+: score) 
      for(int l=0;l<M;l++){
        score += scores[l] = subfunction_contribution(l, N, K, solution);
      }
    return score;
  }

  double parallel_hill_climber( const int M, const int N, const int K, const int T, char solution[N], double scores[N], char states[N], const int nb_colors, int rand_colors[nb_colors], 
                                int color_ascendant[N], int thread_works_bound[MAX_COLOR_OFFSET], double score, int * rnd_gen){
    double * scores_flip = malloc(M*sizeof(double));
    bool keep_going = true;
    double score_thread_update = 0.0;
    while (keep_going) {
      keep_going = false;
      shuffle(nb_colors, rand_colors, rnd_gen);
      #pragma omp parallel num_threads(T) reduction(|: keep_going) reduction(+: score_thread_update)
      {
        const int thread_id = omp_get_thread_num();
        for (int c=0; c<nb_colors; c++) {
          const int color_idx = rand_colors[c]*T; 
          for (int i=thread_works_bound[color_idx+thread_id]; i<thread_works_bound[color_idx+thread_id+1] ;i++) {
            const int v = color_ascendant[i];
            if (states[v] == 0){
              solution[v] = (~solution[v]) & 1;
              double new_sum = .0;
              double old_sum = .0;
              for (int j=v_i_head[v]; j<v_i_head[v] + v_i_degree[v]; j++) {
                const int l = v_i[j];
                new_sum += scores_flip[l] = subfunction_contribution(l, N, K, solution);
                old_sum += scores[l];
              }
              const double delta = new_sum - old_sum;
              if (delta > 0.0) {
                for (int j=v_i_head[v]; j<v_i_head[v] + v_i_degree[v]; j++) {
                  const int l = v_i[j];
                  scores[l] = scores_flip[l];
                }
                score_thread_update += delta;
                activate_adjacent_variables(v,states);
                keep_going = true;
              }
              else {
                solution[v] = (~solution[v]) & 1;
              }
              states[v] = 1;
            }
          }
          #pragma omp barrier
        }
      }
    }
    return score+score_thread_update;
  }

//}

//{ ILS
  double parallel_iterated_local_search(const int M, const int N, const int K, const int T, const int M_factor, const float perturbation, const double time, int *seed){
    double start,chrono;
    start = omp_get_wtime();
    const int rand_number_flips = (int)(perturbation*N);
    char * solution = calloc(N,1);
    char * best_solution = calloc(N, 1);
    int * rand_array = malloc(N*sizeof(int));
    for(int i=0;i<N;i++){
      rand_array[i] = i;
      best_solution[i] = solution[i] = thread_safe_rnd(seed,2);
    }
    sequential_vig(M, N, K);
    parallel_coloring(N,T,seed);
    double * scores = calloc(M, sizeof(double));
    char* states =calloc(N, 1);
    double fit;
    double best_fit = compute_score(M,N,K,best_solution);
    chrono = omp_get_wtime()-start;
    //printf("%i;%i;%i;%i;%f;%f;%f;par\n",M_factor,K,F,S,perturbation,chrono,best_fit);
    while (chrono<time) {
      start = omp_get_wtime();
      memset(states, 0, N);
      thread_rand_mutation(N, solution, rand_array, rand_number_flips,seed);
      fit = init_parallel_hill_climber(M, N, K, T, solution, scores);
      const int coloration = thread_safe_rnd(seed,T);
      fit = parallel_hill_climber(M, N, K, T, solution, scores, states, number_colors[coloration],&rand_color[coloration*MAX_COLOR_SMALL_OFFSET], &color_ascendant[N*coloration], &thread_works_bound[coloration*MAX_COLOR_OFFSET], fit, seed);
      if(fit>best_fit) {
        best_fit = fit;
        memcpy(best_solution, solution, N);
        chrono += omp_get_wtime()-start;
        //printf("%i;%i;%i;%i;%f;%f;%f;par\n",M_factor,K,F,S,perturbation,chrono,best_fit);
      }
      else{
        memcpy(solution, best_solution, N);
        chrono += omp_get_wtime()-start;
      }
    }
    //printf("%i;%i;%i;%i;%f;%f;%f;par\n",M_factor,K,F,S,perturbation,chrono,best_fit);
    free(best_solution);
    free(solution);
    free(rand_array);
    free(scores);
    free(states); 
    free_coloration();
    return best_fit;
  }
//}

void mpi_ils(const int M_factor, const int M, const int N, const int K, const int seed_functions, const int seed_fit, const int T_MPI, const float perturbation, const double time){
  int rank, size;
  double best_score = 0;
  double other_score;
  MPI_Init(NULL, NULL);
  MPI_Comm_rank (MPI_COMM_WORLD, &rank) ;
  MPI_Comm_size (MPI_COMM_WORLD, &size) ;
  if(rank != 0){
    int seed = rank;
    mk_landscape(M,N,K,0,seed_functions,seed_fit);
    other_score = parallel_iterated_local_search(M,N,K,T,M_factor,perturbation,time,&seed);
    //other_score = (double)rank;
    MPI_Send(&other_score, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    free_landscape();
  }
  else{
    for(int i=1;i<size;i++){
      MPI_Recv(&other_score, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      if(best_score<other_score){
        best_score = other_score;
      }
      printf("%f\n",best_score);
      free_landscape();
    }
  }
  MPI_Finalize();
}

int main(int argc, char *argv[]){
  int i;
  i=1;
  const int M_factor        = (argc>i && argv[i][0])? atoi(argv[i]) : 4;      i++;
  const int N               = (argc>i && argv[i][0])? atoi(argv[i]) : 100000; i++;
  const int K               = (argc>i && argv[i][0])? atoi(argv[i]) : 5;      i++;
  const int seed_functions  = (argc>i && argv[i][0])? atoi(argv[i]) : 1;      i++;
  const int seed_fit        = (argc>i && argv[i][0])? atoi(argv[i]) : 1;      i++;
  const int seed_solution   = (argc>i && argv[i][0])? atoi(argv[i]) : 1;      i++;
  const int T               = (argc>i && argv[i][0])? atoi(argv[i]) : 1;      i++;
  const float perturbation  = (argc>i && argv[i][0])? atof(argv[i]) : 0.05;      i++;
  const double time         = (argc>i && argv[i][0])? atof(argv[i]) : 0.05;      i++;

  const int M = N*M_factor;

  int seed = seed_solution;

  mk_landscape(M,N,K,0,seed_functions,seed_fit);
  //parallel_iterated_local_search(M,N,K,T,M_factor,seed_functions,seed_solution,perturbation,time,&seed);
  //test_parallel_hill_climber_mk(M,N,K,T,seed_functions,seed_solution,M_factor);
  //test_hamming_ball_distance_1(M,N,K,seed_functions,seed_solution,M_factor);
  //free_landscape();
  mpi_ils(M_factor,M,N,K,seed_functions,seed_fit,6,0.05,60);
  return 0;
}

// gcc -fopenmp -O2 -o main main.c
