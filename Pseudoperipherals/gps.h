#ifndef GPS_H
#define GPS_H

#include <assert.h>
#include <limits.h>
#include <float.h>
#include "../CommonFiles/graph_parallel.h"
#include "../CommonFiles/linked_list.h"

//typedef enum { START, END } PERIPHERAL_NODES;
typedef enum {FALSE = 0, TRUE} boolean; 

graph_diameter* GPS_get_pseudoperipherals(const METAGRAPH* meta_graph);
BFS*    	    GPS_build_BFS(const METAGRAPH* mgraph, int root);
void 		    GPS_BFS(const METAGRAPH* mgraph, int root, int** levels);
int GPS_count_nodes_by_level(const int* levels, const int n_nodes, int** counts);
#endif