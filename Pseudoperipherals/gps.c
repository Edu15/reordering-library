#include "gps.h"


/*static GRAPH* (*strategy[3])(GRAPH* nodes, int* length) = {
	GRAPH_shrinking_strategy_half_sorted,
	GRAPH_shrinking_strategy_vertex_by_degree,
	GRAPH_shrinking_strategy_five_non_adjacent
};*/


graph_diameter* GPS_get_pseudoperipherals(const METAGRAPH* meta_graph)
{
	//printf("\n");
	int iter = 2;
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
//        printf("root = %d\n", v);
		if (v == NON_ELEMENT) {
			//printf("Error: non element found\n");
			exit(1);
		}

//		bfs = GRAPH_parallel_build_BFS(meta_graph, v);
        bfs = GPS_build_BFS(meta_graph, v);
		checked[v] = TRUE;
//		GRAPH_parallel_print_BFS(bfs);

		bfs_height = bfs->height;
//		printf(" height = %d\n", bfs_height);
		int num_nodes_last_level = bfs->num_nodes_at_level[bfs_height];
		GRAPH* last_level = (GRAPH*) bfs->vertices_at_level[bfs_height];
//        printf("Ultimo nivel:\n");
		for (int i = 0; i < num_nodes_last_level; i++) {
            int w = last_level[i].label;
//            printf("%d, ", w);
//            if (ARRAY_LIST_contains(exploration_ls, w) == FALSE) {
//            	printf("nao contem %d\n", w);
//            } else {
//            	printf("contem %d\n", w);
//            }

            if (checked[w] == FALSE && ARRAY_LIST_contains(exploration_ls, w) == FALSE) {
                ARRAY_LIST_insert(&exploration_ls, w);
            }
        }
//        printf("\n");

		if (num_nodes_last_level <= 0) {
            //printf("Error: last level empty\n");
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

BFS* GPS_build_BFS(const METAGRAPH* mgraph, int root)
{
    int n_nodes, max_level;
    BFS* bfs;
    int* levels;
    int* counts;

    n_nodes = mgraph->size;
    bfs     = malloc(sizeof(BFS));
	levels = calloc(n_nodes, sizeof(int));

	GPS_BFS(mgraph, root, &levels);

    // Preeche o array "counts" com o num de nos por nivel. Os indices sao deslocados em uma unidade para a direita
    max_level = GPS_count_nodes_by_level(levels, n_nodes, &counts);

    // decrementing two levels added to max_level by count_nodes_by_level
    bfs->height = max_level - 2;
    bfs->width  = 0;

    GRAPH graph_node;
    int level, count_level, node;
    bfs->vertices_at_level  = calloc((bfs->height + 1), sizeof(GRAPH*));
    bfs->num_nodes_at_level = calloc((bfs->height + 1), sizeof(int));

    // Allocating memory for BFS structure
    for (count_level = 1; count_level < max_level; ++count_level)
    {
        level = count_level - 1;

        bfs->vertices_at_level[level]  = calloc(counts[count_level], sizeof(GRAPH));
        bfs->num_nodes_at_level[level] = counts[count_level];

        if (counts[count_level] > bfs->width)
            bfs->width = counts[count_level];
    }

    // Filling BFS structure with respective datas
    for (node = 0; node < n_nodes; ++node)
    {
        level       = levels[node];
        count_level = level + 1;

        graph_node.label  = node;
        graph_node.degree = mgraph->graph[node].degree;
        // [num do nivel][nume de vertices no nivel - 1]
        bfs->vertices_at_level[level][counts[count_level]-1] = graph_node; // ????? graph_node dinamicamente alocado
        counts[count_level]--;
    }

    free(levels);
    free(counts);

    return bfs;
}

int GPS_count_nodes_by_level(const int* levels, const int n_nodes, int** counts) {
    int node, count_thread, level;
    int* count;
    int max_level;

    max_level = 0;

    count = calloc(n_nodes, sizeof(int));

    for (node = 0; node < n_nodes; ++node)
    {
        ++count[levels[node]];
        max_level = maximun(max_level, levels[node]);
    }

    max_level += 2;

    *counts = calloc(max_level, sizeof(int));
    (*counts)[0] = 0;

    // Atribui ao array externo e desloca para direita
    for (level = max_level; level > 0; --level)
    {
        (*counts)[level] += count[level-1];
    }

    free(count);

//    printf("\nCONTAGEM DE NIVEIS: ");
//	for (int i = 0; i < max_level; i++) {
//		printf(" %d ", (*counts)[i]);
//	}
//	printf("\n");

    return max_level;
}


void GPS_BFS(const METAGRAPH* mgraph, int root, int** levels)
{
    //double time = omp_get_wtime();

    int node, n_nodes, tail, head, node_degree, active_node, level, count_nodes, adj_node, has_unreached_nodes;
    int* work_set;
    int ws_size;
    int* neighboors;
    BFS* bfs;

    n_nodes = mgraph->size;
    ws_size = n_nodes+1; // oversizing estimate

    bfs     = malloc(sizeof(BFS));

    for (node = 0; node < n_nodes; ++node) {
        (*levels)[node] = INFINITY_LEVEL;
    }

    (*levels)[root] = 0; // O num do nivel raiz eh "0"
    tail = head = 0;

    work_set = calloc(ws_size, sizeof(int));
    QUEUE_enque(&work_set, ws_size, &tail, root); // Adiciona o no "0" na fila para ser explorado
    has_unreached_nodes = n_nodes - 1; // Conta que um no ja foi verificado

    while (!QUEUE_empty(work_set, head, tail)) {
        active_node = QUEUE_deque(&work_set, ws_size, &head);
        //printf("Active node: %d\n", active_node);

        // Get neighboors
        node_degree = mgraph->graph[active_node].degree;
        neighboors  = GRAPH_neighboors(mgraph->mat, active_node, node_degree);
        level       = (*levels)[active_node] + 1; // The number of the new level
        //printf("%d vizinhos\n", node_degree);

        for (count_nodes = 0; count_nodes < node_degree; ++count_nodes) {
            // Para cada vizinho
            adj_node = neighboors[count_nodes];
            //printf(" %d ", adj_node);

            // Se vizinho nunca foi exploradog
            if (level < (*levels)[adj_node]) {
                if ((*levels)[adj_node] == INFINITY_LEVEL) {
                    --has_unreached_nodes;
                }

                (*levels)[adj_node] = level;

                // Adiciona vizinho em fila de cache
                QUEUE_enque(&work_set, n_nodes, &tail, adj_node);
            }
        }

        free(neighboors);
    }

    free(work_set);

//    printf("\nLEVELS: ");
//	for (int i = 0; i < n_nodes; i++) {
//		printf(" %d ", (*levels)[i]);
//	}
//	printf("\n");
    //printf("Tempo com BFS: %lf seg\n", (omp_get_wtime() - time));
}



/*
graph_diameter* GPS_get_pseudoperipherals(const METAGRAPH* meta_graph)
{
	printf("\n");
	int u, v;
	BFS* bfs;
	int bfs_height_u, bfs_height_v;
	int dimension = meta_graph->size;
	u = meta_graph->vertex_min_degree;

	graph_diameter* peripherals = (graph_diameter*) calloc(1, sizeof(graph_diameter));

	bfs = GRAPH_parallel_build_BFS(meta_graph, u);
	bfs_height_u = bfs->height;
	printf(" height = %d\n", bfs_height_u);
	GRAPH* last_level = (GRAPH*) bfs->vertices_at_level[bfs_height_u];
	v = last_level[0].label;
	GRAPH_parallel_destroy_BFS(bfs);

	bfs = GRAPH_parallel_build_BFS(meta_graph, v);
	bfs_height_v = bfs->height;
	printf(" height = %d\n", bfs_height_v);
	last_level = (GRAPH*) bfs->vertices_at_level[bfs_height_v];


	peripherals->start = v;
	peripherals->end = last_level[0].label;
	peripherals->distance = bfs_height_v;

	GRAPH_parallel_destroy_BFS(bfs);
	return peripherals;
}*/




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