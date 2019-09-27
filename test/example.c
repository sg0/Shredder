/*---------------------------------------------------------------------------------
 *  This is a test program to check shredder functions
 *	Expected input: graph_type(0, 1, 2, 3), decompostion_type{0,1}
 *					nrows (>1), ncols (>1)
 *	Expected output: On success decompostion time
 *---------------------------------------------------------------------------------
 */

#include <stdio.h>

#include <stdlib.h>

#include "shredder.h"

int main()
{
	int nrows, ncols, gtype, dtype, nnz, count = 0;

	int *rowidx, *colptr, *decomposition = NULL, **triangles = NULL;

	printf("Graph Type?[0-3] Decompostion Type?{0,1} nrows? ncols?\n");
	scanf("%d%d%d%d", &gtype, &dtype, &nrows, &ncols);

	if (shredder_create_graph(gtype, nrows, ncols, &nnz, &rowidx, &colptr) != SHREDDER_TOKEN_SUCCESS) {
		
		printf("*** Graph Creation Failed ***\n");
		exit(EXIT_FAILURE);	
	}

	shredder_read_graph(gtype, nrows, ncols, nnz, rowidx, colptr);

	//shredder_print_graph();

	free(rowidx);

	rowidx = NULL;

	free(colptr);

	colptr = NULL;

	if (dtype == SHREDDER_DTYPE_STAR 
		&& shredder_stars(&count, &decomposition, &rowidx, &colptr) == SHREDDER_TOKEN_SUCCESS) {
		
		printf("Star Decompostion Time: %lf seconds\n", shredder_decomposition_time());
		
		//shredder_print_stars(count, decomposition, rowidx, colptr);
		
		if (shredder_delete_stars(&count, &decomposition, &rowidx, &colptr) != SHREDDER_TOKEN_SUCCESS ){
			
			printf("*** Stars deletion Failed ***\n");
			exit(EXIT_FAILURE);		
		}
	}
	else if (dtype == SHREDDER_DTYPE_CLIQUE 
		&& shredder_cliques(&count, &triangles) == SHREDDER_TOKEN_SUCCESS) {

		printf("Clique Decomposition Time: %lf seconds\n", shredder_decomposition_time());
		
		//shredder_print_cliques(count, triangles);

		if (shredder_delete_cliques(&count, &triangles) != SHREDDER_TOKEN_SUCCESS) {

			printf("*** Cliques deletion Failed ***\n");
			exit(EXIT_FAILURE);	
		}
	}
	else {

		printf("*** Decomposition Failed ***\n");
		exit(EXIT_FAILURE);	
	}

	if (shredder_delete_graph() != SHREDDER_TOKEN_SUCCESS) {

		printf("*** Graph deletion Failed ***\n");
		exit(EXIT_FAILURE);	
	}
	return 0;
}