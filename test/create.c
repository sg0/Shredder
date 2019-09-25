/*---------------------------------------------------------------------------------
 *  This is a test program to check shredder function: shredder_create_graph
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

	int *rowidx, *colptr;

	printf("Graph Type?[0-3] nrows? ncols?\n");

	scanf("%d%d%d", &gtype,&nrows,&ncols);

	if (shredder_create_graph(gtype, nrows, ncols, &nnz, &rowidx, &colptr) != SHREDDER_TOKEN_SUCCESS) {
		
		printf("*** Graph Creation Failed ***\n");
		exit(EXIT_FAILURE);	
	}

	free(rowidx);

	rowidx = NULL;

	free(colptr);

	colptr = NULL;

	return 0;
}