#include "gps.h"


/*static GRAPH* (*strategy[3])(GRAPH* nodes, int* length) = {
	GRAPH_shrinking_strategy_half_sorted,
	GRAPH_shrinking_strategy_vertex_by_degree,
	GRAPH_shrinking_strategy_five_non_adjacent
};*/


graph_diameter* GPS_get_pseudoperipherals(const METAGRAPH* meta_graph)
{
	int iter = 3;
	int v;
	BFS* bfs;
	int bfs_height;
    int dimension = meta_graph->size;
    int root = meta_graph->vertex_min_degree;
	
	graph_diameter* peripherals = (graph_diameter*) calloc(1, sizeof(graph_diameter));
	graph_diameter* new_peripherals = (graph_diameter*) calloc(1, sizeof(graph_diameter));
    boolean* checked = (boolean*) malloc  (dimension * sizeof(boolean));
    for (int i = 0; i < dimension; i++) {
        checked[i] = FALSE;
    }

	ARRAY_LIST* exploration_ls;
	ARRAY_LIST_init(&exploration_ls);
	ARRAY_LIST_insert(&exploration_ls, root); // Adiciona no raiz a fila
	
	while (iter > 0 && !ARRAY_LIST_is_empty(&exploration_ls)) 
	{
//        ARRAY_LIST_print(exploration_ls);

		v = ARRAY_LIST_remove_first(&exploration_ls);
		if (v == NON_ELEMENT) {
			printf("Error: non element found\n");
			exit(1);
		}
		
		bfs = GRAPH_parallel_build_BFS(meta_graph, v);
		checked[v] = TRUE;
//		GRAPH_parallel_print_BFS(bfs);

		bfs_height = bfs->height;
		int num_nodes_last_level = bfs->num_nodes_at_level[bfs_height];
		GRAPH* last_level = (GRAPH*) bfs->vertices_at_level[bfs_height];
		for (int i = 0; i < num_nodes_last_level; i++) {
            int w = last_level[i].label;

//            if (ARRAY_LIST_contains(exploration_ls, w) == FALSE) {
//            	printf("nao contem %d\n", w);
//            } else {
//            	printf("contem %d\n", w);
//            }

            if (checked[w] == FALSE && ARRAY_LIST_contains(exploration_ls, w) == FALSE) {
                ARRAY_LIST_insert(&exploration_ls, w);
            }
        }


		if (num_nodes_last_level <= 0) {
            printf("Error: last level empty\n");
            exit(1);
		}
        new_peripherals->start = v;
		new_peripherals->end = last_level[0].label;
		new_peripherals->distance = bfs_height;
		if (new_peripherals->distance > peripherals->distance) {
			peripherals->start = new_peripherals->start;
			peripherals->end = new_peripherals->end;
			peripherals->distance = new_peripherals->distance;
		}
		
		iter--;
        GRAPH_parallel_destroy_BFS(bfs);
	}

	free(checked);
	free(new_peripherals);
	ARRAY_LIST_destroy(&exploration_ls);
	return peripherals;
}




/*graph_diameter* GPS_get_pseudoperipherals(const METAGRAPH* meta_graph)
{
	// Create two breadth first search engines
	BFS* forwardBFS;
	BFS* reverseBFS;
	graph_diameter* diameter;
	int local_diameter, size_cand_set, min_width, cand, candidate, swap;
	GRAPH* candidate_set;
	int shrink_type = FIVE_NON_ADJACENT;
	
	diameter = malloc(sizeof(graph_diameter));
	forwardBFS = reverseBFS = NULL;
	
	// Initialize start and end vertices of pseudo-diameter
	diameter->start = meta_graph->vertex_min_degree;
	diameter->end   = NON_VERTEX;
	
	do 
	{
		if (forwardBFS != NULL) GRAPH_parallel_destroy_BFS(forwardBFS);
		
		// do BFS starting at start node 'start'
		forwardBFS = GRAPH_parallel_build_BFS(meta_graph, diameter->start);
		
		// get candidate set of end nodes
		local_diameter = forwardBFS->height;
		candidate_set  = forwardBFS->vertices_at_level[local_diameter];
		size_cand_set  = 1;//forwardBFS->num_nodes_at_level[local_diameter];
		
		// shrink candidate set to a manageable number
		if (shrink_type == FIVE_NON_ADJACENT)
		{
			for (cand = 0; cand < size_cand_set; ++cand)
				candidate_set[cand].neighboors = GRAPH_adjacent(meta_graph->mat, candidate_set[cand].label);
		}
		
		candidate_set = strategy[shrink_type](candidate_set, &size_cand_set);
		
		min_width = INT_MAX;
		
		for (cand = 0; cand < size_cand_set; ++cand)
		{
			candidate = candidate_set[cand].label;
			
			// do BFS from each candidate
			if (reverseBFS != NULL) GRAPH_parallel_destroy_BFS(reverseBFS);
			reverseBFS = GRAPH_parallel_build_BFS(meta_graph, candidate);
			
			if (reverseBFS->width < min_width)
			{
				if (reverseBFS->height > local_diameter)
				{
					// reverseBFS is better than the forwardBFS
					// reset algorithm with candidate as new start
					diameter->start = candidate;
					diameter->end   = NON_VERTEX;
					cand = size_cand_set; // break
				}
				else
				{
					// reverseBFS is narrower than any others
					// make this new end node
					min_width = reverseBFS->width;
					diameter->end = candidate;
				}
			}
		}
		
		if (shrink_type == FIVE_NON_ADJACENT)
		{
			for (cand = 0; cand < size_cand_set; ++cand)
				free(candidate_set[cand].neighboors);
		}
		
		free(candidate_set);
		
	} while (diameter->end == NON_VERTEX);
	
	// swap start & end if the reverseBFS is narrower than forwardBFS
	if (forwardBFS->width > reverseBFS->width)
	{
		swap = diameter->start;
		diameter->start = diameter->end;
		diameter->end   = swap;
	}
	
	GRAPH_parallel_destroy_BFS(forwardBFS);
	GRAPH_parallel_destroy_BFS(reverseBFS);
	
	diameter->distance = local_diameter;
	return diameter;
}*/