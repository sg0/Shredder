/*---------------------------------------------------------------------------------
 *  This is a test program to check shredder functions: shredder_cliques
 *													    shredder_decomposition_time
 *													    shredder_delete_cliques
 *	Expected input: graph_type(0, 1, 2, 3), nrows (>1), ncols (>1)
 *	Expected output: On success nothing
 *---------------------------------------------------------------------------------
 */

#include <stdio.h>

#include <stdlib.h>

#include "shredder.h"

int main()
{
	int nrows, ncols, gtype, nnz, count = 0;

	int *rowidx, *colptr, **triangles = NULL;

	//printf("Graph Type?[0-3] nrows? ncols?\n");

	//scanf("%d%d%d", &gtype,&nrows,&ncols);

	gtype = 1; nrows = 50; ncols = 50;

	if (shredder_create_graph(gtype, nrows, ncols, &nnz, &rowidx, &colptr) != SHREDDER_TOKEN_SUCCESS) {
		
		printf("*** Graph Creation Failed ***\n");
		exit(EXIT_FAILURE);	
	}

	if (shredder_read_graph(gtype, nrows, ncols, nnz, rowidx, colptr) != SHREDDER_TOKEN_SUCCESS) {

		exit(EXIT_FAILURE);	
	}
	
	free(rowidx);

	rowidx = NULL;

	free(colptr);

	colptr = NULL;

	if (shredder_cliques(&count, &triangles) == SHREDDER_TOKEN_SUCCESS) {

		printf("Clique Decomposition Time: %lf seconds\n", shredder_decomposition_time());
		
		//shredder_print_cliques(count, triangles);
				
		if (shredder_delete_cliques(&count, &triangles) != SHREDDER_TOKEN_SUCCESS) {

			exit(EXIT_FAILURE);	
		}
	}
	else {
		
		shredder_delete_graph();	
		
		exit(EXIT_FAILURE);
	}

	if (shredder_delete_graph() != SHREDDER_TOKEN_SUCCESS) {

		exit(EXIT_FAILURE);
	}

	return 0;
}