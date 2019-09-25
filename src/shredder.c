#include <stdio.h>

#include <stdlib.h>

#include <time.h>

#include "shredder.h"


static int graph_status = SHREDDER_TOKEN_INVALID;

static int graph_type = SHREDDER_TOKEN_UNKNOWN;

static int dimension_row = SHREDDER_TOKEN_INVALID;

static int dimension_column = SHREDDER_TOKEN_INVALID;

static int csr_nnz = SHREDDER_TOKEN_INVALID;

static int* csr_edge_pointers = NULL;

static int* csr_edges = NULL;

static double decomposition_time = SHREDDER_TOKEN_INVALID;


static int validate_read_graph(int gtype, int nrows, int ncols, int nnz, int* rowidx, int* colptr)
{
	if (gtype != SHREDDER_GTYPE_RG && gtype != SHREDDER_GTYPE_RG_NWSE 
		&& gtype != SHREDDER_GTYPE_RG_SWNE && gtype != SHREDDER_GTYPE_RG_D) {
		
		return SHREDDER_TOKEN_UNKNOWN;
	}

	if (nrows <= 1 || ncols <= 1 || nnz <= 0 || rowidx == NULL || colptr == NULL) {

		return SHREDDER_TOKEN_INVALID;

	}

	return SHREDDER_TOKEN_SUCCESS;
}

status shredder_read_graph(int gtype, int nrows, int ncols, int nnz, int* rowidx, int* colptr)
{
	int nvertices;

	graph_status = validate_read_graph(gtype, nrows, ncols, nnz, rowidx, colptr);

	if (graph_status!= SHREDDER_TOKEN_SUCCESS) {

		return graph_status;
	}

	graph_type = gtype;

	dimension_row = nrows;

	dimension_column = ncols;

	csr_nnz = nnz;

	nvertices = nrows * ncols;

	csr_edge_pointers = malloc(sizeof(int) * (nvertices + 1));

	if (!csr_edge_pointers) {

		return SHREDDER_TOKEN_FAILURE;
	}

	for (int i = 0; i <= nvertices; i++) {
		
		csr_edge_pointers[i] = rowidx[i];
	}

	csr_edges = malloc(sizeof(int) * nnz);

	if (!csr_edges) {

		return SHREDDER_TOKEN_FAILURE;
	}

	for (int i = 0; i < nnz; i++) {

		csr_edges[i] = colptr[i];
	}

	return SHREDDER_TOKEN_SUCCESS;
}

static int nstars(int gtype, int nrows, int ncols)
{
	int sz = SHREDDER_TOKEN_INVALID;

	if (gtype == SHREDDER_GTYPE_RG) {

		sz = (nrows * ncols) / 2;
	}
	else if (gtype == SHREDDER_GTYPE_RG_SWNE || gtype == SHREDDER_GTYPE_RG_NWSE) {

		sz = 2 * (ncols / 3) * nrows;

		if (ncols % 3) {

			sz++;

			if (nrows >= 3) {
				sz++;
			}
			if (nrows > 3) {
				sz += ((nrows - 2) / 3) + ((nrows - 3) / 3);
			}
		}
		if (ncols % 3 == 2) {

			sz++;

			if (nrows >= 3) {
				sz++;
			}
			if (nrows > 3) {
				sz += ((nrows - 1) / 3) + ((nrows - 3) / 3);
			}
		}
	}
	else if (gtype == SHREDDER_GTYPE_RG_D) {
		sz = (nrows / 2) * ncols + (nrows / 2 + nrows % 2) * (ncols / 2);
	}

	return sz;
}

static int nedges_stars(int gtype, int nrows, int ncols)
{
	int nedges = SHREDDER_TOKEN_INVALID;

	if (gtype == SHREDDER_GTYPE_RG) {
		nedges = csr_nnz / 2;
	}
	else if (gtype == SHREDDER_GTYPE_RG_NWSE) {

		nedges = 0;

		for (int i = 1; i <= nrows; i++) {
			for (int j = 1; j <= ncols; j++) {
				if (i % 3 != j % 3) {
					nedges += (i > 1 && j > 1) + (i > 1) + (j > 1) 
					+ (j < ncols) + (i < nrows) + (i < nrows && j < ncols);
				}
			}
		}
	}
	else if (gtype == SHREDDER_GTYPE_RG_SWNE) {

		nedges = 0;

		for (int i = 1; i <= nrows; i++) {
			for (int j = 1; j <= ncols; j++) {
				if (i % 3 != j % 3) {
					nedges += (i > 1 && j < ncols) + (i > 1) + (j > 1) 
					+ (j < ncols) + (i < nrows) + (i < nrows && j > 1);
				}
			}
		}
	}
	else if (gtype == SHREDDER_GTYPE_RG_D) {

		nedges = 0;

		for (int i = 1; i <= nrows; i++) {
			for (int j = 1; j <= ncols; j++) {
				if (i % 2 == 0 || j % 2 == 0) {
					nedges += (i > 1 && j > 1) + (i > 1) + (i > 1 && j < ncols) + (j > 1) 
					+ (j < ncols) + (i < nrows) + (i < nrows && j > 1) + (i < nrows && j < ncols);
				}
			}
		}
	}

	return nedges;
}

static void stars_RG(int* decomposition, int* rowidx, int* colptr)
{
	int index, index_d = 0, index_r = 0, nedges;

	rowidx[index_r] = 0;

	for (int i = 1; i <= dimension_row; i++) {
		
		for (int j = 1; j <= dimension_column; j++) {

			if (i % 2 != j % 2) {

				index = ((i-1) * dimension_column + j) - 1;

				decomposition[index_d++] = index;

				index_r++;

				nedges = csr_edge_pointers[index + 1] - csr_edge_pointers[index];

				rowidx[index_r] = rowidx[index_r - 1] + nedges;

				for (int k = 0; k < nedges; k++) {

					colptr[rowidx[index_r-1] + k] = csr_edges[csr_edge_pointers[index] + k];
				}			
			}
		}
	}
}

static void stars_RG_NWSE(int* decomposition, int* rowidx, int* colptr)
{
	int index, index_d = 0, index_r = 0, nedges;

	rowidx[index_r] = 0;

	for (int i = 1; i <= dimension_row; i++) {
		
		for (int j = 1; j <= dimension_column; j++) {

			if (i % 3 != j % 3) {

				index = ((i-1) * dimension_column + (dimension_column - j + 1)) - 1;

				decomposition[index_d++] = index;

				index_r++;

				nedges = csr_edge_pointers[index + 1] - csr_edge_pointers[index];

				rowidx[index_r] = rowidx[index_r - 1] + nedges;

				for (int k = 0; k < nedges; k++) {

					colptr[rowidx[index_r-1] + k] = csr_edges[csr_edge_pointers[index] + k];
				}			
			}
		}
	}
}

static void stars_RG_SWNE(int* decomposition, int* rowidx, int* colptr)
{
	int index, index_d = 0, index_r = 0, nedges;

	rowidx[index_r] = 0;

	for (int i = 1; i <= dimension_row; i++) {
		
		for (int j = 1; j <= dimension_column; j++) {

			if (i % 3 != j % 3) {

				index = ((i-1) * dimension_column + j) - 1;

				decomposition[index_d++] = index;

				index_r++;

				nedges = csr_edge_pointers[index + 1] - csr_edge_pointers[index];

				rowidx[index_r] = rowidx[index_r - 1] + nedges;

				for (int k = 0; k < nedges; k++) {

					colptr[rowidx[index_r-1] + k] = csr_edges[csr_edge_pointers[index] + k];
				}			
			}
		}
	}
}

static void stars_RG_D(int* decomposition, int* rowidx, int* colptr)
{
	int index, index_d = 0, index_r = 0, nedges;

	rowidx[index_r] = 0;

	for (int i = 1; i <= dimension_row; i++) {
		
		for (int j = 1; j <= dimension_column; j++) {

			if (i % 2 == 0 || j % 2 == 0) {

				index = ((i-1) * dimension_column + j) - 1;

				decomposition[index_d++] = index;

				index_r++;

				nedges = csr_edge_pointers[index + 1] - csr_edge_pointers[index];

				rowidx[index_r] = rowidx[index_r - 1] + nedges;

				for (int k = 0; k < nedges; k++) {

					colptr[rowidx[index_r-1] + k] = csr_edges[csr_edge_pointers[index] + k];
				}			
			}
		}
	}
}


status shredder_stars(int* size, int** decomposition, int** rowidx, int** colptr)
{
	clock_t time_start = clock();

	if (graph_type == SHREDDER_TOKEN_UNKNOWN) {

		return SHREDDER_TOKEN_FAILURE;
	}

	*size = nstars(graph_type, dimension_row, dimension_column);

	*decomposition = malloc(sizeof(int) * (*size));

	*rowidx = malloc(sizeof(int) * (*size + 1));

	*colptr = malloc(sizeof(int) * nedges_stars(graph_type, dimension_row, dimension_column));

	if (graph_type == SHREDDER_GTYPE_RG) { 
		
		stars_RG(*decomposition, *rowidx, *colptr);		
	}
	else if (graph_type == SHREDDER_GTYPE_RG_NWSE) { 
		
		stars_RG_NWSE(*decomposition, *rowidx, *colptr);		
	}
	else if (graph_type == SHREDDER_GTYPE_RG_SWNE) { 
		
		stars_RG_SWNE(*decomposition, *rowidx, *colptr);		
	}
	else if (graph_type == SHREDDER_GTYPE_RG_D) { 
		
		stars_RG_D(*decomposition, *rowidx, *colptr);		
	}

	clock_t time_end = clock();

	decomposition_time = (double)(time_end - time_start)/CLOCKS_PER_SEC;

	return SHREDDER_TOKEN_SUCCESS;
}

static void SG_NWSE_T(int n, int m, int sr, int sc, int dsr, int dsc, int** triangle_set, int* index)
{
	int der, dec;
	if (n > m) {
		der = dsr + m - 1;
		dec = m;

		for (int i = sr; i < der; i++) {
			//for (int j = sc; j <= min(i - dsr + 1, dec - 1); j++) {
			for (int j = sc; j <= i - dsr + 1; j++) {
				
				triangle_set[*index][0] = ((i-1) * m + j) - 1;
				triangle_set[*index][1] = (i * m + j) - 1;
				triangle_set[*index][2] = (i * m + j + 1) - 1;

				(*index)++;
			}
			//for (int j = max(i - dsr + 1, sc); j < dec; j++) {
			for (int j = i - dsr + 1; j < dec; j++) {

				triangle_set[*index][0] = ((i-1) * m + j) - 1;
				triangle_set[*index][1] = ((i-1) * m + j + 1) - 1;
				triangle_set[*index][2] = (i * m + j + 1) - 1;

				(*index)++;
			}
		}
	}
	else {
		der = n;
		dec = dsc + n - 1;

		for (int i = sr; i < der; i++) {
			//for (int j = sc; j < min(i + dsc, dec); j++) {
			for (int j = sc; j < i + dsc; j++) {

				triangle_set[*index][0] = ((i-1) * m + j) - 1;
				triangle_set[*index][1] = (i * m + j) - 1;
				triangle_set[*index][2] = (i * m + j + 1) - 1;

				(*index)++;
			}
			
			int tmax = ((i + dsc -1) > sc) ? (i + dsc -1) : sc; //max(i + dsc - 1, sc)

			for (int j = tmax; j < dec; j++) {
				
				triangle_set[*index][0] = ((i-1) * m + j) - 1;
				triangle_set[*index][1] = ((i-1) * m + j + 1) - 1;
				triangle_set[*index][2] = (i * m + j + 1) - 1;

				(*index)++;
			}
		}
	}
}

static void Triangulation_NWSE(int n, int m, int *size, int*** decomposition)
{
	int p, q, r;

	int** triangle_set, index = 0, sz;

	p = (n > m) ? n : m; //max(n, m);
	q = (n < m) ? n : m; //min(n, m);
	r = (p-1) - (q - 1) * ((p - 1) / (q - 1));

	sz = (p-1) * q;

	triangle_set = malloc(sizeof(int*) * sz);

	for (int i = 0; i < sz; i++) {
		triangle_set[i] = malloc(sizeof(int) * 3);
	}

	if (n > m) {
		for (int i = 1; i <= (p - 1) / (q - 1); i++) {
			
			SG_NWSE_T(n, m, (i - 1) * (m - 1) + 1, 1, (i - 1) * (m - 1) + 1, 1, 
				triangle_set, &index);
		}

	}
	else {
		for (int i = 1; i <= (p - 1) / (q - 1); i++) {
			
			SG_NWSE_T(n, m, 1, (i - 1) * (n - 1) + 1, 1, (i - 1) * (n - 1) + 1, 
				triangle_set, &index);
		}
	}

	if (n > m && r > 0) {

		SG_NWSE_T(n, m, n - r, 1, n - m + 1, 1, triangle_set, &index);

	}
	else if (m > n && r > 0) {

		SG_NWSE_T(n, m, 1, m - r, 1, m - n + 1, triangle_set, &index);
	}
	
	*size = sz;
	*decomposition = triangle_set;
}

static void SG_SWNE_T(int n, int m, int sr, int sc, int dsr, int dsc, int** triangle_set, int* index)
{
	int der, dec;
	if (n > m) {
		der = dsr + m - 1;
		dec = m;

		for (int i = sr; i < der; i++) {
			for (int j = sc; j <= dec - (i - dsr + 1); j++) {

				triangle_set[*index][0] = ((i-1) * m + j) - 1;
				triangle_set[*index][1] = ((i-1) * m + j + 1) - 1;
				triangle_set[*index][2] = (i * m + j) - 1;

				(*index)++;
			}
			int tmax = ((dec - (i - dsr + 1) + 1) > (sc + 1)) ? (dec - (i - dsr + 1) + 1) : (sc + 1); 
			//max(dec - (i - dsr + 1) + 1, sc + 1)

			for (int j = tmax; j <= dec; j++) {

				triangle_set[*index][0] = ((i-1) * m + j) - 1;
				triangle_set[*index][1] = (i * m + j) - 1;
				triangle_set[*index][2] = (i * m + j - 1) - 1;

				(*index)++;
			}
		}
	}
	else {
		der = n;
		dec = dsc + n - 1;

		for (int i = sr; i < der; i++) {
			for (int j = sc; j <= dec - i; j++) {

				triangle_set[*index][0] = ((i-1) * m + j) - 1;
				triangle_set[*index][1] = ((i-1) * m + j + 1) - 1;
				triangle_set[*index][2] = (i * m + j) - 1;

				(*index)++;
			}
			
			int tmax = ((dec - i + 1) > (sc + 1)) ? (dec - i + 1) : (sc + 1);
			//max(dec - i + 1, sc + 1)
			
			for (int j = tmax; j <= dec; j++) {

				triangle_set[*index][0] = ((i-1) * m + j) - 1;
				triangle_set[*index][1] = (i * m + j) - 1;
				triangle_set[*index][2] = (i * m + j - 1) - 1;

				(*index)++;
			}
		}
	}
}

static void Triangulation_SWNE(int n, int m, int *size, int*** decomposition)
{
	int p, q, r;

	int** triangle_set, index = 0, sz;

	p = (n > m) ? n : m; //max(n, m);
	q = (n < m) ? n : m; //min(n, m);
	r = (p-1) - (q - 1) * ((p - 1) / (q - 1));

	sz = (p-1) * q;

	triangle_set = malloc(sizeof(int*) * sz);

	for (int i = 0; i < sz; i++) {
		triangle_set[i] = malloc(sizeof(int) * 3);
	}

	if (n > m) {
		for (int i = 1; i <= (p - 1) / (q - 1); i++) {
			
			SG_SWNE_T(n, m, (i - 1) * (m - 1) + 1, 1, (i - 1) * (m - 1) + 1, 1, triangle_set, &index);
		}

	}
	else {
		for (int i = 1; i <= (p - 1) / (q - 1); i++) {
			
			SG_SWNE_T(n, m, 1, (i - 1) * (n - 1) + 1, 1, (i - 1) * (n - 1) + 1, triangle_set, &index);
		}

	}

	if (n > m && r > 0) {

		SG_SWNE_T(n, m, n - r, 1, n - m + 1, 1, triangle_set, &index);
	}
	else if (m > n && r > 0) {

		SG_SWNE_T(n, m, 1, m - r, 1, m - n + 1, triangle_set, &index);
	}

	*size = sz;
	*decomposition = triangle_set;
}

status shredder_cliques(int* size, int*** decomposition) {

	if (graph_type != SHREDDER_GTYPE_RG_NWSE && graph_type != SHREDDER_GTYPE_RG_SWNE
		&& graph_type != SHREDDER_GTYPE_RG_D) {

		printf("*** Invalid Graph Type for Clique Decomposition ***\n");

		return SHREDDER_TOKEN_INVALID;
	}
	if (!size || !decomposition) {

		return SHREDDER_TOKEN_INVALID;
	}

	clock_t time_start = clock();

	if (graph_type == SHREDDER_GTYPE_RG_NWSE) {

		Triangulation_NWSE(dimension_row, dimension_column, size, decomposition);
	}
	else if (graph_type == SHREDDER_GTYPE_RG_SWNE) {

		Triangulation_SWNE(dimension_row, dimension_column, size, decomposition);
	}
	else if (graph_type == SHREDDER_GTYPE_RG_D) {

		int size_NWSE, size_SWNE, **ts_NWSE, **ts_SWNE, sz, index = 0;

		Triangulation_NWSE(dimension_row, dimension_column, &size_NWSE, &ts_NWSE);

		Triangulation_SWNE(dimension_row, dimension_column, &size_SWNE, &ts_SWNE);

		sz = size_NWSE + size_SWNE;

		*decomposition = malloc(sizeof(int*) * sz);

		for (int i = 0; i < sz; i++) {
			(*decomposition)[i] = malloc(sizeof(int) * 3);
		}

		for (int i = 0; i < size_NWSE; i++) {

			for (int j = 0; j < 3; j++) {

				(*decomposition)[index][j] = ts_NWSE[i][j];
			}
			index++;
		}

		for (int i = 0; i < size_SWNE; i++) {

			for (int j = 0; j < 3; j++) {

				(*decomposition)[index][j] = ts_SWNE[i][j];
			}
			index++;
		}

		shredder_delete_cliques(&size_NWSE, &ts_NWSE);

		shredder_delete_cliques(&size_SWNE, &ts_SWNE);

		*size = sz;
	}
	else {
		return SHREDDER_TOKEN_INVALID;
	}

	clock_t time_end = clock();

	decomposition_time = (double)(time_end - time_start)/CLOCKS_PER_SEC;

	return SHREDDER_TOKEN_SUCCESS;
}


static void create_RG(int nrows, int ncols, int* rowidx, int* colptr)
{
	int index = 0, index_e = 0, count_edges;

	rowidx[index] = 0;	

	for (int i = 1; i <= nrows; i++) {

		for (int j = 1; j <= ncols; j++) {
			
			count_edges = 0;
			
			if (i > 1) {

				colptr[index_e++] = ((i-2) * ncols + j) - 1;
				count_edges++;

			}
			if (j > 1) {
				
				colptr[index_e++] = ((i-1) * ncols + (j-1)) - 1;
				count_edges++;
			}
			if (j < ncols) {

				colptr[index_e++] = ((i-1) * ncols + (j+1)) - 1;
				count_edges++;
			}
			if (i < nrows) {
				colptr[index_e++] = (i * ncols + j) - 1;
				count_edges++;
			}

			index++;
			rowidx[index] = count_edges + rowidx[index - 1];
		}
	}

}


static void create_RG_NWSE(int nrows, int ncols, int* rowidx, int* colptr)
{
	int index = 0, index_e = 0, count_edges;

	rowidx[index] = 0;	

	for (int i = 1; i <= nrows; i++) {

		for (int j = 1; j <= ncols; j++) {
			
			count_edges = 0;
			
			if (i > 1 && j > 1) {

				colptr[index_e++] = ((i-2) * ncols + (j - 1)) - 1;
				count_edges++;
			}
			if (i > 1) {

				colptr[index_e++] = ((i-2) * ncols + j) - 1;
				count_edges++;

			}
			if (j > 1) {
				
				colptr[index_e++] = ((i-1) * ncols + (j-1)) - 1;
				count_edges++;
			}
			if (j < ncols) {

				colptr[index_e++] = ((i-1) * ncols + (j+1)) - 1;
				count_edges++;
			}
			if (i < nrows) {
				colptr[index_e++] = (i * ncols + j) - 1;
				count_edges++;
			}
			if (i < nrows && j < ncols) {

				colptr[index_e++] = (i * ncols + (j + 1)) - 1;
				count_edges++;
			}

			index++;
			rowidx[index] = count_edges + rowidx[index - 1];
		}
	}

}


static void create_RG_SWNE(int nrows, int ncols, int* rowidx, int* colptr)
{
	int index = 0, index_e = 0, count_edges;

	rowidx[index] = 0;	

	for (int i = 1; i <= nrows; i++) {

		for (int j = 1; j <= ncols; j++) {
			
			count_edges = 0;
			
			if (i > 1) {

				colptr[index_e++] = ((i-2) * ncols + j) - 1;
				count_edges++;

			}
			if (i > 1 && j < ncols) {

				colptr[index_e++] = ((i-2) * ncols + (j + 1)) - 1;
				count_edges++;
			}
			if (j > 1) {
				
				colptr[index_e++] = ((i-1) * ncols + (j-1)) - 1;
				count_edges++;
			}
			if (j < ncols) {

				colptr[index_e++] = ((i-1) * ncols + (j+1)) - 1;
				count_edges++;
			}
			if (i < nrows && j > 1) {

				colptr[index_e++] = (i * ncols + (j - 1)) - 1;
				count_edges++;
			}
			if (i < nrows) {
				colptr[index_e++] = (i * ncols + j) - 1;
				count_edges++;
			}

			index++;
			rowidx[index] = count_edges + rowidx[index - 1];
		}
	}

}

static void create_RG_D(int nrows, int ncols, int* rowidx, int* colptr)
{
	int index = 0, index_e = 0, count_edges;

	rowidx[index] = 0;	

	for (int i = 1; i <= nrows; i++) {

		for (int j = 1; j <= ncols; j++) {
			
			count_edges = 0;
			
			if (i > 1 && j > 1) {

				colptr[index_e++] = ((i-2) * ncols + (j - 1)) - 1;
				count_edges++;
			}
			if (i > 1) {

				colptr[index_e++] = ((i-2) * ncols + j) - 1;
				count_edges++;

			}
			if (i > 1 && j < ncols) {

				colptr[index_e++] = ((i-2) * ncols + (j + 1)) - 1;
				count_edges++;
			}
			if (j > 1) {
				
				colptr[index_e++] = ((i-1) * ncols + (j-1)) - 1;
				count_edges++;
			}
			if (j < ncols) {

				colptr[index_e++] = ((i-1) * ncols + (j+1)) - 1;
				count_edges++;
			}
			if (i < nrows && j > 1) {

				colptr[index_e++] = (i * ncols + (j - 1)) - 1;
				count_edges++;
			}
			if (i < nrows) {
				colptr[index_e++] = (i * ncols + j) - 1;
				count_edges++;
			}
			if (i < nrows && j < ncols) {

				colptr[index_e++] = (i * ncols + (j + 1)) - 1;
				count_edges++;
			}

			index++;
			rowidx[index] = count_edges + rowidx[index - 1];
		}
	}

}

status shredder_create_graph(int gtype, int nrows, int ncols, int* nnz, int** rowidx, int** colptr)
{

	if ((gtype != SHREDDER_GTYPE_RG && gtype != SHREDDER_GTYPE_RG_SWNE 
		&& gtype != SHREDDER_GTYPE_RG_NWSE && gtype != SHREDDER_GTYPE_RG_D) 
		|| nrows <= 1 || ncols <= 1) {

		return SHREDDER_TOKEN_INVALID;
	}

	*rowidx = (int*)malloc(sizeof(int) * (nrows * ncols +1));

	if (!(*rowidx)) {
		return SHREDDER_TOKEN_FAILURE;
	}


	if (gtype == SHREDDER_GTYPE_RG) {

		*nnz = 2 * (2 * nrows * ncols - (nrows + ncols));
	}
	else if (gtype == SHREDDER_GTYPE_RG_NWSE || gtype == SHREDDER_GTYPE_RG_SWNE) {

		*nnz = 2 * (3 * nrows * ncols - 2 * (nrows + ncols) + 1);
	}
	else if (gtype == SHREDDER_GTYPE_RG_D){
		
		*nnz = 2 * (4 * nrows * ncols - 3 * (nrows + ncols) + 2);
	}

	*colptr = (int*)malloc(sizeof(int) * (*nnz));

	if (!(*colptr)) {
		free(*rowidx);
		return SHREDDER_TOKEN_FAILURE;
	}

	if (gtype == SHREDDER_GTYPE_RG) {

		create_RG(nrows, ncols, *rowidx, *colptr);
	}
	else if (gtype == SHREDDER_GTYPE_RG_NWSE) {

		create_RG_NWSE(nrows, ncols, *rowidx, *colptr);
	}
	else if (gtype == SHREDDER_GTYPE_RG_SWNE) {

		create_RG_SWNE(nrows, ncols, *rowidx, *colptr);
	}
	else if (gtype == SHREDDER_GTYPE_RG_D){
		
		create_RG_D(nrows, ncols, *rowidx, *colptr);
	}

	return SHREDDER_TOKEN_SUCCESS;
}

status shredder_print_graph()
{
	if (graph_status != SHREDDER_TOKEN_SUCCESS) {

		printf("\n*** Invalid Graph ***\n");

		return graph_status;
	}

	int nvertices = dimension_row * dimension_column;

	switch (graph_type) {
		case SHREDDER_GTYPE_RG: printf("\nGTYPE_RG");
				break;
		case SHREDDER_GTYPE_RG_NWSE: printf("\nGTYPE_RG_NWSE");
				break;
		case SHREDDER_GTYPE_RG_SWNE: printf("\nGTYPE_RG_SWNE");
				break;
		case SHREDDER_GTYPE_RG_D: printf("\nGTYPE_RG_D");
				break;
		default:
				printf("\nUnknown Graph Type");
				return SHREDDER_TOKEN_FAILURE;	
	}

	printf("\nRows: %d, Columns: %d\n", dimension_row, dimension_column);

	for (int i = 0; i < nvertices; i++) {

		printf("%d: ",i);

		for (int j = csr_edge_pointers[i]; j < csr_edge_pointers[i+1]; j++) {

			printf("%d ", csr_edges[j]);
		}

		printf("\n");
	}
	printf("\n");

	return SHREDDER_TOKEN_SUCCESS;
}


status shredder_print_stars(int size, int* decomposition, int* rowidx, int* colptr)
{
	if (size <= 0 || !decomposition || !rowidx || !colptr) {

		printf("\n*** Invalid Stars ***\n");

		return SHREDDER_TOKEN_INVALID;
	}

	printf("\n# Stars: %d\n", size);

	for (int i = 0; i < size; i++) {

		printf("%d: ",decomposition[i]);

		for (int j = rowidx[i]; j < rowidx[i+1]; j++) {

			printf("%d ", colptr[j]);
		}

		printf("\n");
	}
	printf("\n");

	return SHREDDER_TOKEN_SUCCESS;
}

status shredder_print_cliques(int size, int** decomposition)
{
	if (size <= 0 || !decomposition) {

		return SHREDDER_TOKEN_INVALID;
	}

	printf("\n*** Decomposition Size = %d ****\n", size);
	
	for (int i = 0; i < size; i++) {
		
		for (int j = 0; j < 3; j++) {
			
			printf("%d ", decomposition[i][j]);
		}

		printf("\n");
	}

	return SHREDDER_TOKEN_SUCCESS;
}

double shredder_decomposition_time()
{
	return decomposition_time;
}

status shredder_delete_graph()
{
	if (graph_status != SHREDDER_TOKEN_SUCCESS) {
		printf("\n*** No Valid Graph Stored in Shredder ***\n");

		return SHREDDER_TOKEN_UNKNOWN;
	}
	graph_status = SHREDDER_TOKEN_INVALID;

	graph_type = SHREDDER_TOKEN_UNKNOWN;

	dimension_row = SHREDDER_TOKEN_INVALID;

	dimension_column = SHREDDER_TOKEN_INVALID;

	csr_nnz = SHREDDER_TOKEN_INVALID;

	decomposition_time = SHREDDER_TOKEN_INVALID;

	free(csr_edge_pointers);
	csr_edge_pointers = NULL;

	free(csr_edges);
	csr_edges = NULL;

	return SHREDDER_TOKEN_SUCCESS;
}

status shredder_delete_stars(int *size, int** decomposition, int** rowidx, int** colptr)
{

	if (*size <= 0 || !decomposition || !rowidx || !colptr) {
		printf("\n*** Invalid Parameters for Stars Deletion ***\n");

		return SHREDDER_TOKEN_INVALID;
	}

	*size = 0;

	free(*decomposition);
	*decomposition = NULL;

	free(*rowidx);
	*rowidx = NULL;
	
	free(*colptr);
	*colptr = NULL;

	return SHREDDER_TOKEN_SUCCESS;
}


status shredder_delete_cliques(int *size, int ***decomposition)
{
	if (*size <= 0 || !decomposition) {
		printf("\n*** Invalid Parameters for Cliques Deletion ***\n");

		return SHREDDER_TOKEN_INVALID;
	}

	for (int i = 0; i < (*size); i++) {
		free((*decomposition)[i]);
	}
	
	free(*decomposition);
	
	*decomposition = NULL;
	
	*size = 0;
	
	return SHREDDER_TOKEN_SUCCESS;
}
