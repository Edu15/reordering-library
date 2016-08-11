/*
 * Copyright 2016 Thiago Nascimento <nascimenthiago@gmail.com>
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */
#include <omp.h>
#include "graph.h"
#include "util_parallel.h"

#ifndef __GRAPH_PARALLEL_H__
#define __GRAPH_PARALLEL_H__

#define BFS_WORK_CHUNK 1024
#define BFS_PERCENT_CHUNK 0.5

#define INFINITY_LEVEL 2147483647
#define NON_VERTEX -1


typedef struct
{
	int num_children;
	GRAPH* children;
} genealogy;

typedef struct
{
	int start;
	int end;
} graph_diameter;


typedef struct
{
	MAT* mat;
	GRAPH* graph;
	int vertex_min_degree;
	int size;
} METAGRAPH; 

typedef struct 
{
	int height;
	int width;
	GRAPH** vertices_at_level;
	int* num_nodes_at_level;
} BFS;


typedef enum {
	HALF_SORTED,
	VERTEX_BY_DEGREE,
	FIVE_NON_ADJACENT
} Shrinking_Strategy;

void        	  GRAPH_parallel_fixedpoint_bfs             (MAT* adjacency, int root, int** levels, const float percent_chunk);
void 		  GRAPH_parallel_fixedpoint_BFS		    (const METAGRAPH* mgraph, int root, int** levels, const float percent_chunk);
METAGRAPH* 	  GRAPH_parallel_build_METAGRAPH            (MAT* mat);
void  	  	  GRAPH_parallel_destroy_METAGRAPH          (METAGRAPH* mgraph);
GRAPH*  	  GRAPH_shrinking_strategy_half_sorted      (GRAPH* nodes, int* length);
GRAPH* 	  	  GRAPH_shrinking_strategy_vertex_by_degree (GRAPH* nodes, int* length);
GRAPH* 	  	  GRAPH_shrinking_strategy_five_non_adjacent(GRAPH* nodes, int* length);
BFS*    	  GRAPH_parallel_build_BFS	       	    (const METAGRAPH* mgraph, int root);
void 	 	  GRAPH_parallel_destroy_BFS	            (BFS* bfs);
graph_diameter*   GRAPH_parallel_pseudodiameter             (const METAGRAPH* meta_graph, Shrinking_Strategy type);


#endif