int omp_get_max_threads(){
 return(1);
}

int omp_get_thread_num(){
 return(0);
}

void omp_set_num_threads(){
 error("OpenMP was compiled out.\n");
}
