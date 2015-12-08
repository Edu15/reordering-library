#include <omp.h>
#include "protos.h"

#define NUM_THREADS 4
#define THREAD_ON 1
#define THREAD_OFF 0

#define iseven(n) ((n)%(2)==(0)?(1):(0))
#define isdivisor(d, n) ((n)%(d)==(0)?(1):(0))

typedef enum OPERATION {READ, WRITE} OPERATION;

typedef struct 
{
	int initial_prefix_sum;
	int curr_prefix_sum;
	int curr_total_sum;
	int last_prefix_sum;
	int last_total_sum;
} status_prefix_sum;


typedef struct
{
	int* mem;
	int size;
	int free_position;
} mem_write_next_level;


extern void	GRAPH_parallel_fixedpoint_bfs(MAT* adjacency, int root, int** levels);
extern void 	prefix_sum(const int* counts, int** sums, const int max_level);