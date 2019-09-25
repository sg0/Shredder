/*---------------------------------------------------------------------------------
 *  This is a test program to check shredder functions: shredder_stars
 *													    shredder_decomposition_time
 *													    shredder_delete_stars
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

	int *rowidx, *colptr, *decomposition = NULL;

	printf("Graph Type?[0-3] nrows? ncols?\n");

	scanf("%d%d%d", &gtype,&nrows,&ncols);

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

	if (shredder_stars(&count, &decomposition, &rowidx, &colptr) == SHREDDER_TOKEN_SUCCESS) {

		printf("Star Decompostion Time: %lf seconds\n", shredder_decomposition_time());
		
		
		if (shredder_delete_stars(&count, &decomposition, &rowidx, &colptr) != SHREDDER_TOKEN_SUCCESS) {

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