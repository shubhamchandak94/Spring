#ifndef SPRING_QVZ_CLUSTER_H_
#define SPRING_QVZ_CLUSTER_H_

// All of our structures are delcared in lines.h
// Because otherwise it makes this into a definition nightmare

#include "qvz/include/lines.h"
namespace spring {
namespace qvz {
#define MAX_KMEANS_ITERATIONS 1000

// Memory management
struct cluster_list_t *alloc_cluster_list(struct quality_file_t *info);
void free_cluster_list(struct cluster_list_t *);

} // namespace qvz
} // namespace spring

#endif
