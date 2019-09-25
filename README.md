# Table of Contents
1. [Shredder](#shredder)
2. [Download and Installation](#download-install-shredder-instructions)  
3. [API](#api)
3. [Usage](#usage)
&nbsp;

# Shredder

Shredder is a 'c' library to decompose grid graphs into suitable decomposition units.

Classes of graphs considered right now are rectangular grid graphs (RG), rectangular grid graphs with NWSE diagonals (RG_NWSE), rectangular grid graphs with SWNE diagonals (RG_SWNE), rectangular grid graphs with both diagonals (RG_D).

Decomposition units are star graphs and cliques (triangles). Star graphs are 5 point stencil (S4), 7 point stencil (S6), and 9 point stencil (S8).

Each decompostion tries to minimize number of decomposition units, except RG_D to triangles. For decompostion of RG_D to triangles two decompositions (RG_SWNE to triangles and RG_SWNE to triangles) are simiply put together.


Downloading and Installating Shredder
=====================================
To download, and install Shredder run following instructions in terminals.

    cd              
    git clone https://github.com/neo8git/Shredder.git   #Download Shredder
    ./configure                   # generate make files
    make install       			  # install shredder into 
    make                          # compile the code

To check shredder distribution run followng instructions in terminal. This will give a distributable zipped copy of shredder. 

	cd              			  # cd into the shredder directory 
    make distcheck       		  # uninstall shredder

To uninstall shredder run followng instructions in terminal.  

	cd              			  # cd into the shredder directory 
    make uninstall       		  # uninstall shredder
   
To clean up files generated by configure run followng instruction in terminal.  

	cd              			  # cd into the shredder directory 
    make distclean       		  # uninstall shredder


# API

After installation, you can use shredder as an external library by including "shredder.h" in your program, and then compiling with -lshredder flag. API routines with their arguments are as follows.

status shredder_read_graph(int gtype, int nrows, int ncols, int nnz, int* rowidx, int* colptr)

	This function reads a grid graph in CSR format into shredder's memory (will copy the graph into dynamically allocated memory). Caller is reponsible for calling shredder_delete_graph to clear the allocated memory. Succesful invocation will return SHREDDER_TOKEN_SUCCESS.

	- gtype : can only be one of following SHREDDER_GTYPE_ 
								SHREDDER_GTYPE_RG (Rectangular grid graphs)
								SHREDDER_GTYPE_RG_NWSE (Rectangular grid graphs with NWSE diagonal)
								SHREDDER_GTYPE_RG_SWSE (Rectangular grid graphs with SWNE diagonal)
								SHREDDER_GTYPE_RG_D (Rectangular grid graphs with both diagonals)
	- nrows : number of rows in the grid
	- ncols : number of colums in the grid
	- nnz :   number of edges in the grid, shredder_create_graph can be used for regular grid,
			  for grid with missing edges this need to be figured out correctly, 
			  grid is undirectional, so nnz is twice the # edges
	- rowidx : CSR egde pointers
	- colptr : CSR edges

status shredder_stars(int* size, int** decomposition, int** rowidx, int** colptr)

	Decomposes the graph loaded into Shredder's memory into stars, and saves in the arguments passed, the graph in Shredder's memory will be kept intact, the caller is reposible for freeing the memory used for storing the decompositions, once the need is over. This can be conveniently done using shredder_delete_stars

	- size : Address of the variable where size of the decomposition will be stored.
	- decomposition : Address of a pointer where an array of decompositions will be stored. 
					  The array will contain indices of the rowidx corresponding the vertices picked
	- rowidx : similar to CSR format consecutive entries will hold the range of indices into 
			   colptr holding all the edges of a star, caller will need to pass address of a pointer
	- colptr : CSR edges holding only the edges of the stars picked as decomposition, caller will 			   need to pass address of a pointer



	$./Shredder -r <#rows> -c <#columns> [-g <graph_type>] [-d <methods>] [-v] ...

### DISPLAY HELP 
	$./Shredder

### OPTIONs 
		
	<#rows>:  	 Number of rows in the grid
	<#columns>:  Number of columns in the grid
	<graph_type>  :  RG
	               RG_NWSE,
	               RG_SWNE,
	               RG_D,
		       ...
	<methods>   :  S
	               T
	-v          :  # verbose for debug infomation

### EXAMPLES:
	
	./Shredder -r 5 -c 5 -g RG -v
    ./Shredder -r 5000 -c 2000 -g RG_NWSE -d S
    ./Shredder -r 200 -c 400 -g RG_D -d S
    ./Shredder -r 4 -c 6 -g RG_SWNE -d T -v
	
	
### EXAMPLE OUTPUT

	5 X 5 rectangular grid graph in csr: 

	1: 2 6 
	2: 1 3 7 
	3: 2 4 8 
	4: 3 5 9 
	5: 4 10 
	6: 1 7 11 
	7: 2 6 8 12 
	8: 3 7 9 13 
	9: 4 8 10 14 
	10: 5 9 15 
	11: 6 12 16 
	12: 7 11 13 17 
	13: 8 12 14 18 
	14: 9 13 15 19 
	15: 10 14 20 
	16: 11 17 21 
	17: 12 16 18 22 
	18: 13 17 19 23 
	19: 14 18 20 24 
	20: 15 19 25 
	21: 16 22 
	22: 17 21 23 
	23: 18 22 24 
	24: 19 23 25 
	25: 20 24 

	#rows: 5 #cols: 5 Decomposition Size: 12 Decomposition Time: 7e-06

	5 X 5 rectangular grid graph decomposed to S4 (5 point stencil): 

	2: 1 3 7 
	4: 3 5 9 
	6: 1 7 11 
	8: 3 7 9 13 
	10: 5 9 15 
	12: 7 11 13 17 
	14: 9 13 15 19 
	16: 11 17 21 
	18: 13 17 19 23 
	20: 15 19 25 
	22: 17 21 23 
	24: 19 23 25 

	#rows: 5000 #cols: 2000 Decomposition Size: 6666666 Decomposition Time: 0.141501

	#rows: 200 #cols: 400 Decomposition Size: 60000 Decomposition Time: 0.001285

	4 X 6 rectangular grid graph with SWNE diagonal in csr: 

	1: 2 7 8 
	2: 1 3 8 9 
	3: 2 4 9 10 
	4: 3 5 10 11 
	5: 4 6 11 12 
	6: 5 12 
	7: 1 8 13 14 
	8: 1 2 7 9 14 15 
	9: 2 3 8 10 15 16 
	10: 3 4 9 11 16 17 
	11: 4 5 10 12 17 18 
	12: 5 6 11 18 
	13: 7 14 19 20 
	14: 7 8 13 15 20 21 
	15: 8 9 14 16 21 22 
	16: 9 10 15 17 22 23 
	17: 10 11 16 18 23 24 
	18: 11 12 17 24 
	19: 13 20 
	20: 13 14 19 21 
	21: 14 15 20 22 
	22: 15 16 21 23 
	23: 16 17 22 24 
	24: 17 18 23 

	#rows: 4 #cols: 6 Decomposition Size: 20 Decomposition Time: 0.000151

	4 X 6 rectangular grid graph with SWNE diagonal decomposed to triangles: 

	1 2 7 
	2 3 8 
	3 4 9 
	4 10 9 
	7 8 13 
	8 9 14 
	9 15 14 
	10 16 15 
	13 14 19 
	14 20 19 
	15 21 20 
	16 22 21 
	4 5 10 
	5 6 11 
	6 12 11 
	10 11 16 
	11 17 16 
	12 18 17 
	17 23 22 
	18 24 23


&nbsp;  
&nbsp;  
&nbsp;
