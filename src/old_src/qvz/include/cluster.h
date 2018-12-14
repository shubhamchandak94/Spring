#ifndef _CLUSTER_H_
#define _CLUSTER_H_

// All of our structures are delcared in lines.h
// Because otherwise it makes this into a definition nightmare

#include "lines.h"

#define MAX_KMEANS_ITERATIONS 1000

// Memory management
struct cluster_list_t *alloc_cluster_list(struct quality_file_t *info);
void free_cluster_list(struct cluster_list_t *);
#endif
