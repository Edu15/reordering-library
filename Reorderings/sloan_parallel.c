/*
 * Copyright 2016 Thiago Nascimento nascimenthiago@gmail.com
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

#include "sloan_parallel.h"

/*----------------------------------------------------------------------------
 * SLOAN Reordering
 *--------------------------------------------------------------------------*/

void Parallel_Sloan_original(const METAGRAPH* mgraph, int** permutation, int start_node, int end_node)
{
	int num_nodes, next_id, min_priority, max_priority, num_prior_bags, 
	    num_threads, chunk_size, count_threads_on;
	int* distance;
	int** priority;
	int* status;
	int* size_bags;
	int* degree;
	
	num_nodes = mgraph->size;
	distance  = calloc(num_nodes, sizeof (int));
	GRAPH_parallel_fixedpoint_static_BFS(mgraph, end_node, &distance, BFS_PERCENT_CHUNK);
	
	#pragma omp parallel 
	{
		int node, vertex, vertex_degree, neighbor, ngb, prior_bag, 
		    neighbor_degree, far_ngb, far_neighbor, update_far;
		int* neighbors;
		int* far_neighbors;
		
		#pragma omp single 
		{
			num_threads = omp_get_num_threads();
			chunk_size = num_nodes / num_threads;
			chunk_size = chunk_size > 0 ? chunk_size : 1;
		}
		
		#pragma omp single nowait
		*permutation = calloc (num_nodes, sizeof (int));
		
		#pragma omp single nowait
		priority = calloc (num_nodes, sizeof (int*));
		
		#pragma omp single nowait
		status = calloc (num_nodes,sizeof (int));
		
		#pragma omp single nowait
		degree = calloc (num_nodes,sizeof (int));
		
		#pragma omp single nowait
		min_priority = INFINITY_LEVEL;
		
		#pragma omp barrier
		
		#pragma omp for schedule(static, chunk_size)
		for (node = 0; node < num_nodes; ++node)
		{
			status[node]   = INACTIVE;
			degree[node]   = mgraph->graph[node].degree;
			priority[node] = calloc(2, sizeof(int));
			priority[node][SLOAN_CURR_PRIOR] = SLOAN_W1*distance[node] - SLOAN_W2*(degree[node] + 1);
				
			#pragma omp critical
			if (priority[node][SLOAN_CURR_PRIOR] < min_priority) 
				min_priority = priority[node][SLOAN_CURR_PRIOR];
		}
		
		#pragma omp single nowait
		if (min_priority > 0) 
		{
			printf("*** [Error] Sloan minimum initial priority is higher than 0 ***\n");
			exit(1);
		}
		
		#pragma omp single
		{
			num_prior_bags = SLOAN_PRIORITY_FACTOR * 
				(priority[start_node][SLOAN_CURR_PRIOR] - min_priority);
			size_bags = calloc(num_prior_bags, sizeof(LIST*));
		}
		
		#pragma omp for schedule(static, chunk_size)
		for (node = 0; node < num_nodes; ++node) 
		{
			priority[node][SLOAN_CURR_PRIOR] -= min_priority;
			priority[node][SLOAN_NEW_PRIOR]   = priority[node][SLOAN_CURR_PRIOR];
			
			#pragma omp atomic
			++size_bags[priority[node][SLOAN_CURR_PRIOR]];
		}
		
		#pragma omp single nowait
		status[start_node] = PREACTIVE;
		
		#pragma omp single nowait
		next_id = 0;
		
		#pragma omp barrier
		
		/* ************************************
		 * ********** Processing nodes ********
		 * ************************************/
		
		while (next_id < num_nodes)
		{
			#pragma omp sections
			{
				// Chosing maximum priority == Defining the Logical Bag
				#pragma omp section
				{ 	
					for (prior_bag = num_prior_bags - 1; prior_bag >= 0; --prior_bag)
					{
						if (size_bags[prior_bag] > 0)
						{
							max_priority = prior_bag;
							prior_bag = 0; // stop seaching
						}
					}
				}
				
				#pragma omp section
				count_threads_on = 0;
			}
			
			#pragma omp for schedule(static, chunk_size)
			for (vertex = 0; vertex < num_nodes; ++vertex)
			{
				// Processing Logical Bag
				if ((status[vertex] != NUMBERED) && (priority[vertex][SLOAN_CURR_PRIOR] == max_priority))
				{
					neighbors     = mgraph->graph[vertex].neighboors;
					vertex_degree = mgraph->graph[vertex].degree;
					
					#pragma omp atomic
					count_threads_on++;
					
					for (ngb = 0; ngb < vertex_degree; ++ngb)
					{
						update_far = OFF;
						neighbor   = neighbors[ngb];
						
						if (status[vertex] == PREACTIVE)
						{
							if (status[neighbor] == ACTIVE)
							{
								#pragma omp atomic
								priority[neighbor][SLOAN_NEW_PRIOR] += SLOAN_W2;
							}
							else if ((status[neighbor] == INACTIVE) || (status[neighbor] == PREACTIVE))
							{
								#pragma omp atomic
								priority[neighbor][SLOAN_NEW_PRIOR] += SLOAN_W2;
								
								#pragma omp critical
								status[neighbor] = ACTIVE;
								
								update_far = ON;
							}
						}
						else if ((status[vertex] == ACTIVE) && (status[neighbor] == PREACTIVE))
						{
							#pragma omp atomic
							priority[neighbor][SLOAN_NEW_PRIOR] += SLOAN_W2;
							
							#pragma omp critical
							status[neighbor] = ACTIVE;
							
							update_far = ON;
						}
				
						if (update_far)
						{
							/* ***********************
							* Updating far neighbors
							* ***********************
							*/
							
							far_neighbors   = mgraph->graph[neighbor].neighboors;
							neighbor_degree = mgraph->graph[neighbor].degree;
							
							for (far_ngb = 0; far_ngb < neighbor_degree; ++far_ngb) 
							{
								far_neighbor = far_neighbors[far_ngb];
								
								if (status[far_neighbor] == INACTIVE)
								{
									#pragma omp critical
									status[far_neighbor] = PREACTIVE;
									
									#pragma omp atomic
									priority[far_neighbor][SLOAN_NEW_PRIOR] += SLOAN_W2;
								}
								
								#pragma omp atomic
								priority[neighbor][SLOAN_NEW_PRIOR] += SLOAN_W2;
							}
						}
					}
					
					#pragma omp atomic
					--count_threads_on;
					
					while (count_threads_on > 0); // Barrier
					
					#pragma omp critical
					{
						size_bags[priority[vertex][SLOAN_CURR_PRIOR]]--;
						(*permutation)[next_id++] = vertex;
						status[vertex] = NUMBERED;
					}
				}
			}
			
			#pragma omp for schedule(static, chunk_size)
			for (node = 0; node < num_nodes; ++node)
			{
				if ((status[node] != NUMBERED) && 
					(priority[node][SLOAN_CURR_PRIOR] != priority[node][SLOAN_NEW_PRIOR]))
				{
					#pragma omp atomic
					size_bags[priority[node][SLOAN_CURR_PRIOR]]--;
						
					#pragma omp atomic
					size_bags[priority[node][SLOAN_NEW_PRIOR]]++;
					
					priority[node][SLOAN_CURR_PRIOR] = priority[node][SLOAN_NEW_PRIOR];
				}
			}
		}
		
		#pragma omp barrier
		
		#pragma omp single nowait
		free(distance);
		
		#pragma omp single nowait
		free(status);
		
		#pragma omp single nowait
		free(size_bags);
		
		#pragma omp single nowait
		free(degree);
		
		#pragma omp for schedule(static, chunk_size)
		for (node = 0; node < num_nodes; ++node)
			free(priority[node]);
		
		#pragma omp single nowait
		free(priority);
	}
}


static inline void add_bag(int** elements, const int num_nodes, int* tail, const int neighbor)
{
	#pragma omp critical
	QUEUE_enque(elements, num_nodes, tail, neighbor);
}


double Parallel_Sloan_v1(const METAGRAPH* mgraph, int** permutation, int start_node, int end_node)
{
	int num_nodes, next_id, num_threads, chunk_size, min_priority, 
	    max_priority, num_prior_bags, count_threads_on, size_max_bag, 
	    c_node, max_processed_priority;
	int* distance;
	bag* bags_priority;
	int* status;
	int* degree;
	int* priority;
	int* max_bag;
	double time, time_pre, time_pos;
	
	num_nodes = mgraph->size;
	distance  = calloc(num_nodes, sizeof(int));
	GRAPH_parallel_fixedpoint_static_BFS(mgraph, end_node, &distance, BFS_PERCENT_CHUNK);
	
	#pragma omp parallel 
	{
		int node, index_vertex, vertex, vertex_degree, neighbor, ngb, neighbor_degree, 
		    far_ngb, far_neighbor, update_far, th_min_priority, th_max_priority, n_bag,
		    p_bag;
		int* neighbors;
		int* far_neighbors;
		
		#pragma omp single
		time_pre = omp_get_wtime();
		
		#pragma omp sections
		{
			#pragma omp section
			{
				num_threads = omp_get_num_threads();
				chunk_size = num_nodes / num_threads;
				chunk_size = chunk_size > 0 ? chunk_size : 1;
			}
			
			#pragma omp section
			*permutation = calloc(num_nodes, sizeof(int));
			
			#pragma omp section
			status = calloc(num_nodes, sizeof(int));
			
			#pragma omp section
			priority = calloc(num_nodes, sizeof(int));
			
			#pragma omp section
			degree = calloc(num_nodes, sizeof(int));
			
			#pragma omp section
			max_bag = calloc(num_nodes, sizeof(int));
			
			#pragma omp section
			min_priority = INFINITY_LEVEL;
			
			#pragma omp section
			max_priority = -INFINITY_LEVEL;
		}
		
		th_min_priority = INFINITY_LEVEL;
		th_max_priority = -INFINITY_LEVEL;
		
		/* ********************************
		 * ****** Pre-processing **********
		 * ********************************
		 */
		#pragma omp for schedule(static, chunk_size)
		for (node = 0; node < num_nodes; ++node)
		{
			status[node]   = INACTIVE;
			degree[node]   = mgraph->graph[node].degree;
			priority[node] = SLOAN_W1*distance[node] - SLOAN_W2*(degree[node] + 1);
			
			if (priority[node] > th_max_priority) th_max_priority = priority[node];
			if (priority[node] < th_min_priority) th_min_priority = priority[node];
		}
		
		// Defining initial minimum and maximum priority 
		#pragma omp critical
		{
			if (th_max_priority > max_priority) max_priority = th_max_priority;
			if (th_min_priority < min_priority) min_priority = th_min_priority;
		}
		
		#pragma omp barrier
		
		/* ********************************
		 * *** Generating Priority Bags ***
		 * ********************************
		 */
		
		#pragma omp single
		{
			// oversizing estimate
			num_prior_bags = SLOAN_PRIORITY_FACTOR * (priority[start_node] - min_priority); 
			bags_priority = malloc(num_prior_bags * sizeof(bag));
		}
		
		// Initializing priority bags
		#pragma omp for 
		for (n_bag = 0; n_bag < num_prior_bags; ++n_bag)
		{
			bags_priority[n_bag].elements = calloc(num_nodes, sizeof(int));
			bags_priority[n_bag].head = bags_priority[n_bag].tail = 0;
			omp_init_lock(&bags_priority[n_bag].lock);
		}
		
		// Loading priority bags
		#pragma omp for schedule(static, chunk_size)
		for (node = 0; node < num_nodes; ++node)
		{
			priority[node] -= min_priority;
			omp_set_lock(&(bags_priority[priority[node]].lock));
			QUEUE_enque(&(bags_priority[priority[node]].elements), num_nodes, &(bags_priority[priority[node]].tail), node);
			omp_unset_lock(&(bags_priority[priority[node]].lock));
		}
		
		#pragma omp sections
		{
			#pragma omp section
			max_priority -= min_priority;
			
			#pragma omp section
			status[start_node] = PREACTIVE; 
			
			#pragma omp section
			next_id = 0;
			
			#pragma omp section
			count_threads_on = 0;
		}
		
		#pragma omp barrier
		
		#pragma omp single
		{
			time_pre = (omp_get_wtime() - time_pre)/100.0;
			printf("Time pre-processing: %lf\n", time_pre);fflush(stdout);
		}
		
		/* ************************************
		 * ********** Processing nodes ********
		 * ************************************/
		
		#pragma omp single
		time = omp_get_wtime();
		
		// TODO: criar caches de bags priority para somente após o consumo da inteira bag corrente
		// as threads descarregarem na lista de priority bags
		while (next_id < num_nodes)
		{
// 			#pragma omp critical
// 			{
// 				printf("thead %d - running\n", omp_get_thread_num());fflush(stdout);
// 			}
			th_max_priority = -INFINITY_LEVEL;
			
			// Copying current max priority bag
			#pragma omp single	
			{
// 				printf("Priority Bag [%d]: ", max_priority);fflush(stdout);
				size_max_bag = bags_priority[max_priority].tail - bags_priority[max_priority].head;
				c_node = 0;
				
				for (index_vertex = bags_priority[max_priority].head; index_vertex < size_max_bag; ++index_vertex)
				{
					max_bag[c_node++] = bags_priority[max_priority].elements[index_vertex];
// 					printf("%d ", bags_priority[max_priority].elements[index_vertex]);
				}
// 				printf("\n");fflush(stdout);
			}
			
			#pragma omp sections
			{
				// Cleaning current processed max priority bag
				#pragma omp section
				bags_priority[max_priority].tail = 0;	
				
				#pragma omp section
				max_processed_priority = -INFINITY_LEVEL;
			}
			
			#pragma omp atomic
			++count_threads_on;
					
			// Processing maximum priority bag of nodes
			#pragma omp for 
			for (index_vertex = 0; index_vertex < size_max_bag; ++index_vertex)
			{
				vertex 	      = max_bag[index_vertex];
				neighbors     = mgraph->graph[vertex].neighboors;
				vertex_degree = mgraph->graph[vertex].degree;
// 				printf("[th %d]processing node %d(d=%d) with priority %d from position %d\n", omp_get_thread_num(), vertex, vertex_degree, priority[vertex], index_vertex);fflush(stdout);
				
				if (status[vertex] == NUMBERED) 
				{
// 					printf("[th %d]descarding NUMBERED node %d\n", omp_get_thread_num(), vertex);fflush(stdout);
					continue;
				}
				
				for (ngb = 0; ngb < vertex_degree; ++ngb)
				{
					neighbor   = neighbors[ngb];
// 					printf("[th %d]processing neighboor %d\n", omp_get_thread_num(), neighbor);fflush(stdout);
					
					update_far = OFF;
					
					if (status[vertex] == PREACTIVE)
					{
						if (status[neighbor] == ACTIVE)
						{
							omp_set_lock(&mgraph->lock_node[neighbor]);
							priority[neighbor] += SLOAN_W2;
							if (priority[neighbor] > th_max_priority) th_max_priority = priority[neighbor];
							add_bag(&bags_priority[priority[neighbor]].elements, num_nodes, &bags_priority[priority[neighbor]].tail, neighbor);
							omp_unset_lock(&mgraph->lock_node[neighbor]);
						}
						else if ((status[neighbor] == INACTIVE) || (status[neighbor] == PREACTIVE))
						{
							omp_set_lock(&mgraph->lock_node[neighbor]);
							priority[neighbor] += SLOAN_W2;
							if (priority[neighbor] > th_max_priority) th_max_priority = priority[neighbor];
							if ((status[neighbor] == INACTIVE) || (status[neighbor] == PREACTIVE)) status[neighbor] = ACTIVE;
							add_bag(&(bags_priority[priority[neighbor]].elements), num_nodes, &(bags_priority[priority[neighbor]].tail), neighbor);
							omp_unset_lock(&mgraph->lock_node[neighbor]);
							
							update_far = ON;
						}
					}
					else if ((status[vertex] == ACTIVE) && (status[neighbor] == PREACTIVE))
					{
						omp_set_lock(&mgraph->lock_node[neighbor]);
						priority[neighbor] += SLOAN_W2;
						if (priority[neighbor] > th_max_priority) th_max_priority = priority[neighbor];
						if ((status[vertex] == ACTIVE) && (status[neighbor] == PREACTIVE)) status[neighbor] = ACTIVE;
						add_bag(&(bags_priority[priority[neighbor]].elements), num_nodes, &(bags_priority[priority[neighbor]].tail), neighbor);
						omp_unset_lock(&mgraph->lock_node[neighbor]);
						
						update_far = ON;
					}
			
					if (update_far)
					{
						/* ***********************
						* Updating far neighbors
						* ***********************
						*/
						
						far_neighbors   = mgraph->graph[neighbor].neighboors;
						neighbor_degree = mgraph->graph[neighbor].degree;
						
						for (far_ngb = 0; far_ngb < neighbor_degree; ++far_ngb) 
						{
							far_neighbor = far_neighbors[far_ngb];
							
							if (far_neighbor == vertex) continue;
							
// 							printf("[th %d]processing far_neighboor %d\n", omp_get_thread_num(), far_neighbor);fflush(stdout);
							
							if (status[far_neighbor] == INACTIVE)
							{
								omp_set_lock(&mgraph->lock_node[far_neighbor]);
								priority[far_neighbor] += SLOAN_W2;
								if (priority[far_neighbor] > th_max_priority) th_max_priority = priority[far_neighbor];
								if (status[far_neighbor] == INACTIVE) status[far_neighbor] = PREACTIVE;
								add_bag(&(bags_priority[priority[far_neighbor]].elements), num_nodes, &(bags_priority[priority[far_neighbor]].tail), far_neighbor);
								omp_unset_lock(&mgraph->lock_node[far_neighbor]);
							}
							
							omp_set_lock(&mgraph->lock_node[far_neighbor]);
							priority[far_neighbor] += SLOAN_W2;
							if (priority[far_neighbor] > th_max_priority) th_max_priority = priority[far_neighbor];
							add_bag(&(bags_priority[priority[far_neighbor]].elements), num_nodes, &(bags_priority[priority[far_neighbor]].tail), far_neighbor);
							omp_unset_lock(&mgraph->lock_node[far_neighbor]);
						}
					}
				}
				
				// Placing vertex in permutation array
				#pragma omp critical
				{
					if (status[vertex] != NUMBERED)
					{
// 						printf("[th %d] --> permut: nextid %d, vertex %d\n", omp_get_thread_num(), next_id, vertex);fflush(stdout);
						(*permutation)[next_id++] = vertex;
						status[vertex] = NUMBERED;
					}
// 					else
// 					{
// 						printf("[th %d] not permut NUMBERED node %d\n", omp_get_thread_num(), vertex);fflush(stdout);
// 					}
				}
			}
			
			// Defining max processed priority 
			#pragma omp critical
			if (th_max_priority > max_processed_priority) max_processed_priority = th_max_priority;
			
			#pragma omp atomic
			--count_threads_on;
					
			while (count_threads_on > 0); // Barrier
			
			// Defining next max priority bag
			#pragma omp single
			{
// 				printf("------------------------------------------rrr\n");fflush(stdout);
				if (max_processed_priority >= max_priority)
				{
					max_priority = max_processed_priority;
				}
				else
				{
					// Searching nex max priority bag
					for (p_bag = max_priority; p_bag >= 0; --p_bag)
					{
						if (!QUEUE_empty(bags_priority[p_bag], bags_priority[p_bag].head, 
							bags_priority[p_bag].tail))
						{
							// Find out next max priority bag
							max_priority = p_bag;
							p_bag = 0;
						}
					}
				}
			}
		}
		
		#pragma omp barrier
		
		/* ********************************
		 * ****** Post-processing *********
		 * ********************************
		 */
		#pragma omp single
		{
			time = (omp_get_wtime() - time)/100.0;
			printf("Time processing: %lf\n", time);fflush(stdout);
		}
		
// 		#pragma omp single
// 		{
// 			printf("thead %d - permutation: ", omp_get_thread_num());fflush(stdout);
// 			for (node = 0; node < num_nodes; ++node) printf("%d ", (*permutation)[node]);fflush(stdout);
// 			printf("\n");
// 		}
// 		
// 		#pragma omp single
// 		{
// 			printf("thread %d - num_prior_bags %d\n", omp_get_thread_num(), num_prior_bags);
// 			fflush(stdout);
// 		}
		
		#pragma omp single
		time_pos = omp_get_wtime();
		
		#pragma omp for 
		for (n_bag = 0; n_bag < num_prior_bags; ++n_bag)
		{
			free(bags_priority[n_bag].elements);
			omp_destroy_lock(&bags_priority[n_bag].lock);
		}
		
		#pragma omp sections
		{
			#pragma omp section
			free(distance);
			
			#pragma omp section
			free(status);
			
			#pragma omp section
			free(degree);
			
			#pragma omp section
			free(priority);
			
			#pragma omp section
			free(bags_priority);
			
			#pragma omp section
			free(max_bag);
		}
		
		#pragma omp single
		{
			time_pos = (omp_get_wtime() - time_pos)/100.0;
			printf("Time pos-processing: %lf\n", time_pos);fflush(stdout);
		}
	}
	
	return time + time_pre + time_pos;
// 	printf("sloan is over!\n");fflush(stdout);
}



double Parallel_Sloan(METAGRAPH* mgraph, int** permutation, int start_node, int end_node)
{
	int num_nodes, next_id, num_threads, chunk_size, min_priority, max_priority, 
	    num_prior_bags, size_max_bag, max_processed_priority;
	bag* bags_priority;
	int* status;
	int* max_bag;
	double time, time_pre, time_pos;
// 	ARRAY_LIST* priority_queue;
	
	num_nodes = mgraph->size;
	GRAPH_parallel_fixedpoint_sloan_BFS(mgraph, end_node, BFS_PERCENT_CHUNK);
	
	#pragma omp parallel 
	{
		int node, index_vertex, vertex, vertex_degree, neighbor, ngb, neighbor_degree, 
		    far_ngb, far_neighbor, update_far, n_bag, p_bag, tail_bag_change, c_node;
		int* neighbors;
		int* far_neighbors;
		GRAPH* bag_change;
		
		#pragma omp single
		time_pre = omp_get_wtime();
		
		#pragma omp sections
		{
			#pragma omp section
			{
				num_threads = omp_get_num_threads();
				chunk_size = num_nodes / num_threads;
				chunk_size = chunk_size > 0 ? chunk_size : 1;
			}
			
			#pragma omp section
			*permutation = calloc(num_nodes, sizeof(int));
			
			#pragma omp section
			status = calloc(num_nodes, sizeof(int));
			
			#pragma omp section
			max_bag = calloc(num_nodes, sizeof(int));
			
			#pragma omp section
			min_priority = mgraph->min_sloan_priority;
			
			#pragma omp section
			max_priority = mgraph->max_sloan_priority;
			
// 			#pragma omp section
// 			ARRAY_LIST_init(&priority_queue);
		}
		
		/* ********************************
		 * ****** Pre-processing **********
		 * ** Generating Priority Bags **** 
		 * ********************************
		 */
		#pragma omp single
		{
			// Oversizing estimate
			num_prior_bags = SLOAN_PRIORITY_FACTOR * (mgraph->sloan_priority[start_node] - min_priority); 
			bags_priority = malloc(num_prior_bags * sizeof(bag));
		}
		
		// Initializing priority bags - So expensive loop!!!
		#pragma omp for 
		for (n_bag = 0; n_bag < num_prior_bags; ++n_bag)
		{
// 			bags_priority[n_bag].elements = calloc(num_nodes, sizeof(int)); // so so so expensive
// 			bags_priority[n_bag].head = bags_priority[n_bag].tail = 0;
			ARRAY_LIST_init(&bags_priority[n_bag].list_elements);
			omp_init_lock(&bags_priority[n_bag].lock);
		}
		
		// Loading priority bags
		#pragma omp for schedule(static, chunk_size)
		for (node = 0; node < num_nodes; ++node)
		{
			status[node]   = INACTIVE;
			mgraph->sloan_priority[node] -= min_priority;
			omp_set_lock(&(bags_priority[mgraph->sloan_priority[node]].lock));
// 			QUEUE_enque(&(bags_priority[mgraph->sloan_priority[node]].elements), num_nodes, 
// 				    &(bags_priority[mgraph->sloan_priority[node]].tail), node);
			ARRAY_LIST_insert(&bags_priority[mgraph->sloan_priority[node]].list_elements, node);
// 			ARRAY_LIST_add_desc_order(&priority_queue, mgraph->sloan_priority[node]);
			omp_unset_lock(&(bags_priority[mgraph->sloan_priority[node]].lock));
		}
		
		#pragma omp sections
		{
			#pragma omp section
			max_priority -= min_priority;
			
			#pragma omp section
			status[start_node] = PREACTIVE; 
			
			#pragma omp section
			next_id = 0;
		}
		
		// Bag that stores each node change per thread
		bag_change = malloc((num_nodes/2) * sizeof(GRAPH));
		
		#pragma omp single
		{
			time_pre = (omp_get_wtime() - time_pre)/100.0;
			printf("Time pre-processing: %lf\n", time_pre);fflush(stdout);
		}
		
		/* ************************************
		 * ********** Processing nodes ********
		 * ************************************/
		
		#pragma omp single
		time = omp_get_wtime();
		
		while (next_id < num_nodes)
		{
			
// 			#pragma omp atomic
// 			++count_threads_on;
			
			// Cloning current max priority bag
			#pragma omp single	
			{
				size_max_bag = bags_priority[max_priority].list_elements->size;
				c_node       = 0;
				
				while (bags_priority[max_priority].list_elements->size > 0)
				{
					max_bag[c_node++] = 
						ARRAY_LIST_remove_first(&bags_priority[max_priority].list_elements);
				}
				
// 				printf("----------------Priority Bag [%d]: ", max_priority);fflush(stdout);
// 				for (c_node = 0; c_node < size_max_bag; ++c_node)
// 					printf("%d ", max_bag[c_node]);
// 				printf("\n");fflush(stdout);
				
			}
			
			tail_bag_change = 0;
			
			#pragma omp single
			max_processed_priority = -INFINITY_LEVEL;
			
			// Processing maximum priority bag of nodes
			#pragma omp for 
			for (index_vertex = 0; index_vertex < size_max_bag; ++index_vertex)
			{
				vertex = max_bag[index_vertex];
				
				if (status[vertex] == NUMBERED) 
				{
// 					printf("[th %d]descarding NUMBERED node %d\n", omp_get_thread_num(), vertex);fflush(stdout);
					continue;
				}
				
				vertex_degree = mgraph->graph[vertex].degree;
				neighbors     = GRAPH_neighboors(mgraph->mat, vertex, vertex_degree);
				
// 				printf("[th %d]processing node %d(d=%d) with priority %d from position %d\n", omp_get_thread_num(), vertex, vertex_degree, mgraph->sloan_priority[vertex], index_vertex);fflush(stdout);
				
				for (ngb = 0; ngb < vertex_degree; ++ngb)
				{
					neighbor   = neighbors[ngb];
// 					printf("[th %d]processing neighboor %d\n", omp_get_thread_num(), neighbor);fflush(stdout);
					
					update_far = OFF;
					
					if (status[vertex] == PREACTIVE)
					{
						if (status[neighbor] == ACTIVE)
						{
							bag_change[tail_bag_change].label    = neighbor;
							bag_change[tail_bag_change].status   = ACTIVE;
							bag_change[tail_bag_change].distance = mgraph->sloan_priority[neighbor] + SLOAN_W2;
							++tail_bag_change;
						}
						else if ((status[neighbor] == INACTIVE) || (status[neighbor] == PREACTIVE))
						{
							bag_change[tail_bag_change].label    = neighbor;
							bag_change[tail_bag_change].status   = ACTIVE;
							bag_change[tail_bag_change].distance = mgraph->sloan_priority[neighbor] + SLOAN_W2;
							++tail_bag_change;
							
							update_far = ON;
						}
					}
					else if ((status[vertex] == ACTIVE) && (status[neighbor] == PREACTIVE))
					{
						bag_change[tail_bag_change].label    = neighbor;
						bag_change[tail_bag_change].status   = ACTIVE;
						bag_change[tail_bag_change].distance = mgraph->sloan_priority[neighbor] + SLOAN_W2;
						++tail_bag_change;
						
						update_far = ON;
					}
			
					if (update_far)
					{
						/* ***********************
						* Updating far neighbors
						* ***********************
						*/
						
						far_neighbors   = GRAPH_adjacent(mgraph->mat, neighbor);
						neighbor_degree = mgraph->graph[neighbor].degree;
						
						for (far_ngb = 0; far_ngb < neighbor_degree; ++far_ngb) 
						{
							far_neighbor = far_neighbors[far_ngb];
							
							if (far_neighbor == vertex) continue;
							
// 							printf("[th %d]processing far_neighboor %d\n", omp_get_thread_num(), far_neighbor);fflush(stdout);
							bag_change[tail_bag_change].label = far_neighbor;
							
							if (status[far_neighbor] == INACTIVE)
							{
								bag_change[tail_bag_change].status   = PREACTIVE;
								bag_change[tail_bag_change].distance = mgraph->sloan_priority[far_neighbor] + SLOAN_W2;
							}
							else
							{
								bag_change[tail_bag_change].status   = status[far_neighbor];
								bag_change[tail_bag_change].distance = mgraph->sloan_priority[far_neighbor];
							}
							
							bag_change[tail_bag_change].distance += SLOAN_W2;
							++tail_bag_change;
						}
						
						free(far_neighbors);
					}
				}
				
				free(neighbors);
				
				// Placing vertex in permutation array
// 				omp_set_lock(&mgraph->lock_node[vertex]);
				#pragma omp critical
				{
					if (status[vertex] != NUMBERED)
					{
	// 					printf(">>>permut id: %d, vertex %d\n", th_next_id, vertex);fflush(stdout);
						
						(*permutation)[next_id++] = vertex;
						status[vertex] = NUMBERED;
					}
				}
// 				omp_unset_lock(&mgraph->lock_node[vertex]);
			} // priority bag loop
			
			// Unloading each bag change
			if (tail_bag_change > 0)
			{
				#pragma omp critical
				{
					for (index_vertex = 0; index_vertex < tail_bag_change; ++index_vertex)
					{
						vertex = bag_change[index_vertex].label;
						
						if (bag_change[index_vertex].status > status[vertex])
							status[vertex] = bag_change[index_vertex].status;
						
//  						// Defining max processed priority
						if (bag_change[index_vertex].distance > max_processed_priority) 
							max_processed_priority = bag_change[index_vertex].distance;
						
						mgraph->sloan_priority[vertex] = bag_change[index_vertex].distance;
						
						ARRAY_LIST_insert(&bags_priority[mgraph->sloan_priority[vertex]].list_elements, vertex);
					}
				}
			}
			
// 			#pragma omp single
// 			max_priority = ARRAY_LIST_remove_first(&priority_queue);
			
// 			#pragma omp barrier
			
			// Defining next max priority bag
			#pragma omp single
			{
// 				printf("------------------------------------------rrr\n");fflush(stdout);
				if (max_processed_priority >= max_priority)
				{
					max_priority = max_processed_priority;
				}
				else
				{
					// Searching next max priority bag
					for (p_bag = max_priority; p_bag >= 0; --p_bag)
					{
						if (bags_priority[p_bag].list_elements->size > 0)
						{
							// Find out next max priority bag
							max_priority = p_bag;
							p_bag = 0;
						}
					}
				}
			}
		} // main loop - while
		
		free(bag_change);
// 		free(th_max_bag);
		
		#pragma omp barrier
		
		/* ********************************
		 * ****** Post-processing *********
		 * ********************************
		 */
		#pragma omp single
		{
			time = (omp_get_wtime() - time)/100.0;
			printf("Time processing: %lf\n", time);fflush(stdout);
		}
		
// 		#pragma omp single
// 		{
// 			printf("thead %d - permutation: ", omp_get_thread_num());fflush(stdout);
// 			for (node = 0; node < num_nodes; ++node) printf("%d ", (*permutation)[node]);fflush(stdout);
// 			printf("\n");
// 		}
		
		#pragma omp single
		time_pos = omp_get_wtime();
		
		#pragma omp for 
		for (n_bag = 0; n_bag < num_prior_bags; ++n_bag)
		{
			ARRAY_LIST_destroy(&bags_priority[n_bag].list_elements);
			omp_destroy_lock(&bags_priority[n_bag].lock);
		}
		
		#pragma omp sections
		{
			#pragma omp section
			free(status);
			
			#pragma omp section
			free(bags_priority);
			
			#pragma omp section
			free(max_bag);
			
// 			#pragma omp section
// 			ARRAY_LIST_destroy(&priority_queue);
		}
		
		#pragma omp single
		{
			time_pos = (omp_get_wtime() - time_pos)/100.0;
			printf("Time pos-processing: %lf\n", time_pos);fflush(stdout);
		}
	}
	
	return time + time_pre + time_pos;
// 	printf("sloan is over!\n");fflush(stdout);
}


