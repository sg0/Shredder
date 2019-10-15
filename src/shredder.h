/*******************************************************************************
    This file is part of Shredder (a package for graph decomposition), which is 
    under its License protection. You should have received a copy of the License. 
    If not, see <https://github.com/neo8git/Shredder>
*******************************************************************************/


#ifndef SHREDDER_H
#define SHREDDER_H

#define SHREDDER_TOKEN_SUCCESS 0

#define SHREDDER_TOKEN_FAILURE -1

#define SHREDDER_TOKEN_INVALID -2

#define SHREDDER_TOKEN_UNKNOWN -3

#define SHREDDER_TOKEN_FALSE 0

#define SHREDDER_TOKEN_TRUE 1


#define SHREDDER_DTYPE_STAR 0
//Decompostion type stars/stencils

#define SHREDDER_DTYPE_CLIQUE 1
//Decompostion type cliques/triangles


#define SHREDDER_GTYPE_RG 0
//Graph type rectangular grid

#define SHREDDER_GTYPE_RG_NWSE 1
//Graph type rectangular grid with NWSE diagonals

#define SHREDDER_GTYPE_RG_SWNE 2
//Graph type rectangular grid with SWNE diagonals

#define SHREDDER_GTYPE_RG_D 3
//Graph type rectangular grid with both diagonals


typedef int status;
//type definition for returning functions' status


status shredder_read_graph(int gtype, int nrows, int ncols, int nnz, int* rowidx, int* colptr);
/*---------------------------------------------------------------------------------
 *  shredder_read_graph  -  reads a grid graph in CSR format into shredder's memory
 *							this will copy the graph into dynamically allocated
 *							memory, caller is reponsible for calling 
 *							shredder_delete_graph to clear the allocated memory
 *	
 *  - gtype : can only be one of SHREDDER_GTYPE_ defined above
 			  (e.g. SHREDDER_GTYPE_RG_NWSE)
 *
 *  - nrows : number of rows in the grid
 *
 *  - ncols : number of colums in the grid
 *
 *  - nnz :   number of edges in the grid, shredder_create_graph can be used for
 *			  regular grid, for grid with missing edges this need to be figured
 *			  out correctly, grid is undirectional, so nnz is twice the # edges
 *
 *  - rowidx : CSR egde pointers
 *
 *	- colptr : CSR edges
 *---------------------------------------------------------------------------------
 */

status shredder_stars(int* size, int** decomposition, int** rowidx, int** colptr);
/*---------------------------------------------------------------------------------
 *  shredder_stars  -  decompose the graph loaded into Shredder's memory into
 *					   stars, and saves in the arguments passed, the graph in
 *					   Shredder's memory will be kept intact, the caller is reposible
 *					   for freeing the memory used for storing the decompositions,
 *					   once the need is over. This can be conveniently done using 
 *					   shredder_delete_stars
 *	
 *  - size : Address of the variable where size of the decomposition 
 * 			 will be stored.
 *
 *  - decomposition : Address of a pointer where an array of decompositions 
 *		 			  will be stored. The array will contain indices of the 
 *					  rowidx corresponding the vertices picked
 *					  
 *  - rowidx : similar to CSR format consecutive entries will hold the range of 
 * 			   indices into colptr holding all the edges of a star,
 *			   caller will need to pass address of a pointer
 *
 *	- colptr : CSR edges holding only the edges of the stars picked as decomposition,
 *			   caller will need to pass address of a pointer
 *-----------------------------------------------------------------------------------
 */

status shredder_cliques(int* size, int*** decomposition);
/*---------------------------------------------------------------------------------
 *  shredder_cliques  -  decompose the graph loaded into Shredder's memory into
 *					     clqiues (right now only triangles), and saves in the 
 *						 arguments passed, the graph in Shredder's memory will 
 *						 be keet intact, the caller is reposible for freeing the 
 *						 memory used for storing the decompositions, once the need 
 *						 is over. This can be conveniently done using 
 *					   	 shredder_delete_cliques
 *	
 *  - size : Address of the variable where size of the decomposition 
 * 			 will be stored.
 *
 *  - decomposition : Address of a double pointer where a 2D array of decompositions 
 *		 			  will be stored. The array will contain indices of the 
 *					  rowidx corresponding the vertices picked from graph loaded 
 *					  Each row in the array will contain a clique.
 *-----------------------------------------------------------------------------------
 */

double shredder_decomposition_time();
/*---------------------------------------------------------------------------------
 *  shredder_decomposition_time  -  will return the time spend for decomposing a
 *									graph already loaded into Shredder's memoroy.
 *									This function can be used either for star
 *									or clique decomposition. The function should
 *									be called after the call of shredder_stars 
 *									or shredder_cliques. Timing information will
 *									be lost once shredder_delete_graph is called
 *									Unit of time is second.
 *---------------------------------------------------------------------------------
 */

status shredder_delete_graph();
/*---------------------------------------------------------------------------------
 *  shredder_delete_graph  -  clears shredder's memory, shredder can be reused for
 *							  loading a new graph after this operation
 *---------------------------------------------------------------------------------
 */

status shredder_delete_stars(int *size, int** decomposition, int** rowidx, int** colptr);
/*---------------------------------------------------------------------------------------
 *  shredder_delete_stars  -  clears memory allocated for stars through shreder_stars
 *
 *  - size : Address of the variable where size of the decomposition is stored
 *
 *  - decomposition : Address of a pointer where an array of decompositions 
 *		 			  is stored.
 *					  
 *  - rowidx : caller will need to pass address of a pointer
 *
 *	- colptr : caller will need to pass address of a pointer							
 *---------------------------------------------------------------------------------------
 */

status shredder_delete_cliques(int *size, int ***decomposition);
/*---------------------------------------------------------------------------------------
 *  shredder_delete_cliques  -  clears memory allocated for cliquess through shreder_cliques
 *
 *  - size : Address of the variable where size of the decomposition is stored
 *
 *  - decomposition : Address of a double pointer where a 2D array of decompositions 
 *		 			  is stored.
 *---------------------------------------------------------------------------------------
 */

status shredder_create_graph(int gtype, int nrows, int ncols, int* nnz, int** rowidx, int** colptr);
/*--------------------------------------------------------------------------------------------------
 *  shredder_create_graph  -  creates a regular grid graph in CSR format, which can be passed
 *							  to shredder_read_graph. This is a conveient routine for creating 
 *							  a regular grid graph of graph type defined above (GTYPE_...).
 *							  If the grid graph has missing edges, then this should not be used. 
 *	
 *  - gtype : can only be one of SHREDDER_GTYPE_ defined above (e.g. SHREDDER_GTYPE_RG_NWSE)
 *
 *  - nrows : number of rows in the grid
 *
 *  - ncols : number of colums in the grid
 *
 *  - nnz :   Address of an integer variable to store number of edges in the grid
 *
 *  - rowidx : Addres of a pointer to store CSR egde pointers
 *
 *	- colptr : Addres of a pointer to store CSR edges
 *---------------------------------------------------------------------------------
 */

status shredder_print_graph();
/*---------------------------------------------------------------------------------
 *  shredder_print_graph  -  prints graph loaded into shredder's memory
 *---------------------------------------------------------------------------------
 */

status shredder_print_stars(int size, int* decomposition, int* rowidx, int* colptr);
/*---------------------------------------------------------------------------------
 *  shredder_print_stars  -  prints all the stars and associated edges of the stars
 *---------------------------------------------------------------------------------
 */

status shredder_print_cliques(int size, int** decomposition);
/*---------------------------------------------------------------------------------
 *  shredder_print_stars  -  prints all the cliques of the decomposition array,
 *							 right now only prints traingles for each decomposition
 *---------------------------------------------------------------------------------
 */

#endif
