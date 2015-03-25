#include <stdio.h>
#include <stdlib.h>

#define HEURISTIC 0

#define NMAX 1024
#define MMAX ((NMAX +31)/ 32)

typedef enum status_t { SUCCESS, SMALL_GRAPH_ERR, LARGE_GRAPH_ERR, EOF_GRAPH_ERR,
VERT_BOUNDS_ERR, MISSED_VERTEX, ASYMMETRIC_GRAPH_ERR
} status_t;

void print_graph(int n, int [NMAX][MMAX]);
int set_size(int, int []);
void print_set(int, int []);
void handle_error(status_t, int, int, int, int, int);
int set_size(int, int []);
void print_graph(int, int [NMAX][MMAX]);
void print_set(int, int []);
void read_graph (int *, int *, int [NMAX][MMAX], status_t *, int *, int *, int *, int *);
void set_u_blue(int G[][MMAX], int n, int u, int num_dominated[], int n_dominated, int dom[],
				 int* min_size, int min_dom[], int num_choice[], int level, int size, int delta, int num_tests);
void set_u_red(int G[][MMAX], int n, int u, int num_dominated[], int n_dominated, int dom[],
				 int* min_size, int min_dom[], int num_choice[], int level, int size, int delta, int num_tests);
void find_dom_set(int G[][MMAX], int n, int num_dominated[], int n_dominated, int dom[],
				 int* min_size, int min_dom[], int num_choice[], int level, int size, int delta, int num_tests);
void print_solution(int status, int size, int dom_set[]);

/***** COMPRESSED ADJACENCY MATRIX CODE *****/

#define SETWD(pos) ((pos)>>5)    /* number of longword containing bit pos */
#define SETBT(pos) ((pos)&037)   /* position within longword of bit pos */
#define ADD_ELEMENT(setadd,pos) ((setadd)[SETWD(pos)] |= bit[SETBT(pos)])
#define DEL_ELEMENT(setadd,pos) ((setadd)[SETWD(pos)] &= ~bit[SETBT(pos)])
#define IS_ELEMENT(setadd,pos) ((setadd)[SETWD(pos)] & bit[SETBT(pos)])

/* number of 1-bits in longword x */
#define POP_COUNT(x) bytecount[(x)>>24 & 0377] + bytecount[(x)>>16 & 0377] \
                        + bytecount[(x)>>8 & 0377] + bytecount[(x) & 0377]

int bytecount[] =        {0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,
                          1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
                          1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
                          2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
                          1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
                          2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
                          2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
                          3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
                          1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
                          2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
                          2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
                          3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
                          2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
                          3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
                          3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
                          4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8};

int bit[] = {020000000000,010000000000,04000000000,02000000000,01000000000,
              0400000000,0200000000,0100000000,040000000,020000000,010000000,
              04000000,02000000,01000000,0400000,0200000,0100000,040000,020000,
              010000,04000,02000,01000,0400,0200,0100,040,020,010,04,02,01};

/*********************************************/


/******* Test functions for test function table *******/
char n_extra_test(int n, int n_dominated, int* min_size, int num_choice[], int size, int delta)
{

	int undominated, n_extra;

	undominated = n - n_dominated;
	n_extra = (undominated + delta) / (delta + 1);

	if(size + n_extra >= *min_size){
		return 0;
	}

	return 1;
}

char zero_choice_test(int n, int n_dominated, int* min_size, int num_choice[], int size, int delta)
{
	int i;

	for(i = 0; i < n; i++){
		if(num_choice[i] == 0){
			return 0;
		}
	}

	return 1;
}

// The domination number of a graph without isolated vertices is at most n/2. (Ore)
char ore_test(int n, int n_dominated, int* min_size, int num_choice[], int size, int delta)
{
	if(size + 1 > (n>>1))
		return 0;

	return 1;
}

// The domination number of a graph with minimum degree >= 2 and n >= 8 is at most 2n/5. (Blank)
char blank_test(int n, int n_dominated, int* min_size, int num_choice[], int size, int delta)
{
	if(size + 1 > ((n<<1)/5))
		return 0;

	return 1;	
}

// The domination number of a graph with minimum degree >=3 is at most 3n/8. (Reed)
char reed_test(int n, int n_dominated, int* min_size, int num_choice[], int size, int delta)
{
	if(size + 1 > ((3*n)>>3))
		return 0;

	return 1;		
}

/* Function table for the tests. This way, we don't waste time checking whether or not the problem is some kind of graph
at every recursive call. */
char (*tests[3])(int n, int n_dominated, int* min_size, int num_choice[], int size, int delta);


/***********************************************/

// Status 0 means the dominating set is the best so far.
// Status 1 means the dominating set is a minimum dominating set.
// Code created by Wendy Myrvold.
void print_solution(int status, int size, int dom_set[])
{
    int i;

    printf("%1d %3d ", status, size);

    for (i=0; i < size; i++)
        printf(" %3d", dom_set[i]);
    printf("\n");
	
    fflush(stdout);
}

/*
* Set the vertex u to blue. That is, at this step, we've decided that u is not in the dominating set.
* Calls find_dom set and then cleans up after it's done.
*
* Note that set_u_blue and set_u_red are inlined-- this is a hint to the compiler that we should
* not actually perform a function call at the definition of set_u_blue and set_u_red, but instead
* put the function definition in the body of find_dom_set. This way, I can employ lots of heuristics
* in the main algorithm without copying code and (hopefully, as long as the compiler agrees) not take a
* hit on the maximum depth of recursion.
*/
inline void set_u_blue(int G[][MMAX], int n, int u, int num_dominated[], int n_dominated, int dom[],
				 int* min_size, int min_dom[], int num_choice[], int level, int size, int delta, int num_tests)
{
	int i, j;
	char undominatable = 0;

	/***** Try setting u to blue. That is, decide that u is not in the dominating set. *****/
	for(i = 0; i < n; i++){
		if(IS_ELEMENT(G[u], i)){
			if(--(num_choice[i]) == 0){
				undominatable = 1;
				break;
			}
		}
	}

	if (!undominatable)
		find_dom_set(G, n, num_dominated, n_dominated, dom, min_size, min_dom, num_choice, level + 1, size, delta, num_tests);

	// Clean up what we changed so that we can try setting u to red next.
	for(j = 0; j < i; j++){
		if(IS_ELEMENT(G[u], j)){
			++num_choice[j];
		}
	}
}

/*
* Set the vertex u to red. That is, at this step, we've decided that u must be in the dominating set.
* Calls find_dom_set and then cleans up after it's done.
*/
inline void set_u_red(int G[][MMAX], int n, int u, int num_dominated[], int n_dominated, int dom[],
				 int* min_size, int min_dom[], int num_choice[], int level, int size, int delta, int num_tests)
{
	int i;

	/***** Try setting u to red. That is, we decide that u is in the dominating set. *****/
	for(i = 0; i < n; i++){
		if(IS_ELEMENT(G[u], i)){
			if(num_dominated[i] == 0){ // This neighbour wasn't dominated before this step.
				n_dominated++;
			}
			num_dominated[i]++;
		}
	}

	dom[size] = u;
	find_dom_set(G, n, num_dominated, n_dominated, dom, min_size, min_dom, num_choice, level+1, size+1, delta, num_tests);

	// Clean up what we changed-- when we return from this level, we shouldn't have altered anything.
	// We only have to clean up the changes to arrays.
	for(i = 0; i < n; i++){
		if(IS_ELEMENT(G[u], i)){
			num_dominated[i]--;
		}
	}
}

/*
* Checks if we've found a dominating set at this level. If we have and it's better, we print it and return
* 1 signifying that we're done with this branch.
* Otherwise, returns 0 and says to continue.
*/
char found_dom_set(int level, int n_dominated, int n, int min_dom[], int *min_size, int dom[], int size)
{
	int i;

	/*** We dominated every vertex, so we found some dominating set ***/
	if(n_dominated == n || level == n) {

		if (size < *min_size){	// We found a better dominating set than the one we had before.

			for(i = 0; i < size; i++){
				min_dom[i] = dom[i];
			}

			*min_size = size;
			print_solution(0, *min_size, min_dom);
		}

		return 1;
	}

	// Default: Causes the find_dom_set routine to continue.
	return 0;
}

/*
* Exponential backtracker to find the minimum dominating set of a simple graph G.
*/
void find_dom_set(int G[][MMAX], int n, int num_dominated[], int n_dominated, int dom[],
				 int* min_size, int min_dom[], int num_choice[], int level, int size, int delta, int num_tests)
{

	int u, i, deg_u;

	u = level;

	// We went all the way down the branch, but didn't get every vertex.
	if(level == n && n_dominated < n)
		return;

	// We found a dominating set, so we should return.
	if (found_dom_set(level, n_dominated, n, min_dom, min_size, dom, size) == 1)
		return;

	/* Run a series of tests in order to check whether or not we've found a dominating set.
	Tests are stored in an array of functions pass_test. Before running the main program, we make some observations about
	the graph and then decide which tests are appropriate.
	This method allows us to avoid if/else chains in this routine. */
	for (i = 0; i < num_tests; i++){
		if((*tests[i])(n, n_dominated, min_size, num_choice, size, delta) == 0)
			return;
	}

	deg_u = set_size(n, G[u]);

	// If the current vertex is isolated, we need to put it in the dominating set.
	if(deg_u == 1){
		set_u_red(G, n, u, num_dominated, n_dominated, dom, min_size, min_dom, num_choice, level, size, delta, num_tests);
		return;
	}

	// At this point, u isn't in the set but num_choice[u] == 1. This means that the best we can do on this branch is
	// set it to red. Setting it to blue would just cause us to return.
	else if(num_choice[u] == 1){
		set_u_red(G, n, u, num_dominated, n_dominated, dom, min_size, min_dom, num_choice, level, size, delta, num_tests);
		return;
	}

#if HEURISTIC
	/** Greedy branching heuristic
	* If the degree of our current vertex is high, and there aren't "many" ways to dominate it left, we'll try putting it in
	* first.
	*
	* Otherwise, we'll just try setting it blue first, then setting it red.
	* 
	* Motivation: A kind of "human" heuristic. I thought back to how people were reasoning the ice-cream truck example in class.
	*/

	else if(deg_u >= delta>>1 && num_choice[u] <= deg_u>>1){
		set_u_red(G, n, u, num_dominated, n_dominated, dom, min_size, min_dom, num_choice, level, size, delta, num_tests);
		set_u_blue(G, n, u, num_dominated, n_dominated, dom, min_size, min_dom, num_choice, level, size, delta, num_tests);
		return;
	}
#endif

	else {
		set_u_blue(G, n, u, num_dominated, n_dominated, dom, min_size, min_dom, num_choice, level, size, delta, num_tests);
		set_u_red(G, n, u, num_dominated, n_dominated, dom, min_size, min_dom, num_choice, level, size, delta, num_tests);
		return;
	}
}


/* Prints an error to stderr based off an integer code. Exit if necessary. */
void handle_error(status_t err, int graph_number, int n, int m, int cert_len, int missed)
{
	switch(err){
		case SMALL_GRAPH_ERR:
			fprintf(stderr,"Error- Graphs must have n > 0 vertices. \n");
			fprintf(stderr,"Graph %3d: ERROR in graph.\nTerminating.\n", graph_number);
			fflush(stderr);
			exit(err);
			break;
		case LARGE_GRAPH_ERR:
			fprintf(stderr,"Error- n (%d) > NMAX (%d). Please increase the size of NMAX.\n", n, NMAX);
			fprintf(stderr,"Graph %3d: ERROR in graph.\nTerminating.\n", graph_number);
			fflush(stderr);
			exit(err);
			break;
		case VERT_BOUNDS_ERR:
			fprintf(stderr,"Error- Vertex %d's degree was out of the range [%d ... %d].\n", missed, 0, n-1);
			fprintf(stderr,"Graph %3d: ERROR in graph.\nTerminating.\n", graph_number);
			fflush(stderr);
			exit(err);
			break;
		case EOF_GRAPH_ERR:
			fprintf(stderr,"Error- Reached EOF before reading the entire graph.\n");
			fprintf(stderr,"Graph %3d: ERROR in graph.\nTerminating.\n", graph_number);
			fflush(stderr);
			exit(err);
			break;
		case MISSED_VERTEX:
			fprintf(stderr,"Error- Vertex %d is not dominated.\n", missed);
			fprintf(stderr,"Graph %3d: ERROR in certificate.\n", graph_number);
			fflush(stderr);
			break;
		case ASYMMETRIC_GRAPH_ERR:
			fprintf(stderr,"Graph %3d: Non-simple input graph. If (u,v) is in G, then (v,u) must be in G.\n", graph_number);
			fflush(stderr);
			exit(err);
			break;
		default:
			break;
	}
}

// Compute the size of a set.
int set_size(int n, int set[])
{
    int j, m, d;

    m = (n+31)/ 32;

    d=0;

    for (j=0; j < m; j++)
    {
       d+= POP_COUNT(set[j]);
    }
    return(d);
}

// Prints a graph stored in the compressed adjacency format.
void print_graph(int n, int G[NMAX][MMAX])
{
    int i, d;

    for (i = 0; i < n; i++)
    {
        d = set_size(n, G[i]);
        printf("%3d(%1d):", i, d);
        print_set(n, G[i]);
    }
}

// Prints a set.
void print_set(int n, int set[])
{
   int i;

   for (i=0; i < n; i++)
       if (IS_ELEMENT(set, i))
           printf("%3d", i);
   printf("\n");
}

/* 
*   Reads in the graph from stdin and puts it in the compressed adjacency
*   matrix format.
*
* 	Catches some errors along the way. Validate graph is called after to remove
*	any errors unrelated to reading in a graph.
*/
void read_graph (int *n, int *m, int G[NMAX][MMAX], status_t *stat, int *bad, int *min_degree, int *max_degree, int *maxv)
{
	int i, j;
	int deg, v;

	*bad = -1; // a flag storing a bad vertex or degree for error checking
	*stat = SUCCESS; // right now, we have no error

	/*** should we quit on this graph early? ***/
	// there are too few vertices for a valid graph (less than 1 makes no sense)
	if (*n < 1) {
		*stat = SMALL_GRAPH_ERR;
		return;
	}

	// there are too many vertices for the graph
	if (*n > NMAX) {
		*stat = LARGE_GRAPH_ERR;
		return;
	}

	*m = (*n+31)/32; // want to define edges by 32 bit ints
	*min_degree = *n;
	*max_degree = 0;

	/***  graph initialization. all edges are 0. ***/
	for (i = 0; i < *n; i++) {
        for (j = 0; j < *m; j++) {
            G[i][j]= 0;
        }
    }

	/*** reading the graph in and checking for errors as we go. ***/
	for (i = 0; i < *n; i++) {

		// Error: We hit EOF before actually reading the graph despite reading in a graph size?
		if ((scanf("%d", &deg)) == EOF) {
			*stat = EOF_GRAPH_ERR;
			return;
		}

		// Error: Vertex degree out of range.
		if (deg < 0 || deg > NMAX-1) { 
			*stat = VERT_BOUNDS_ERR;
			*bad = deg;
			return;
		}

		if (deg < *min_degree){
			*min_degree = deg;
		}

		if (deg > *max_degree){
			*max_degree = deg;
			*maxv = i;
		}

		// Store vertices in the graph
		for (j = 0; j < deg; j++) {

			// Hit EOF while reading the vertices.
			if (scanf("%d", &v) != 1) {
				*stat = EOF_GRAPH_ERR;
				return;
			}

			// Found a vertex of unsuitable degree.
			if (v < 0 || v > *n-1) {
				*bad = v;
				*stat = VERT_BOUNDS_ERR;
			}

			ADD_ELEMENT(G[i], (int)v);
		}
	}
}

/*
* Checks if an input graph is an undirected graph.
* If a graph isn't valid, we should give up because something could have been
* read incorrectly.
*/
void validate_graph(int G[NMAX][MMAX], int n, status_t *stat)
{
	int i, j;

	// there was already an error in read_graph so we shouldn't bother.
	if (*stat != SUCCESS){
		return;
	}

	for (i = 0; i < n; i++){
		for (j = 0; j < n; j++)
			if (IS_ELEMENT(G[i], j)){
				if (IS_ELEMENT(G[j], i)){
					continue;
				}
				else{
					*stat = ASYMMETRIC_GRAPH_ERR;
					return;
				}
			}
	}
}

void initialize_dom_set(int G[][MMAX], int n, int num_choice[], int num_dominated[], int min_dom[], int *delta)
{
	int i;

	for(i = 0; i < n; i++){
		ADD_ELEMENT(G[i], i); // every vertex dominates itself.
		num_choice[i] = set_size(n, G[i]) + 1;

		// Delta should be the maximum number of vertices we can dominate
		if (num_choice[i] > *delta) {
			*delta = num_choice[i];
		}

		num_dominated[i] = 0;
		min_dom[i] = i;
	}
}

/*
* main
*/
int main(void) {
	status_t stat;
	int m, n, delta, min_size;
	int missed = -1;
	int num_tests = 2;
	int graphs_processed = 0;
	int min_degree, max_degree, maxv;

	tests[0] = n_extra_test;
	tests[1] = zero_choice_test;

	do {
		if(scanf("%d",&n) == EOF) {
			break;
		}

		int G[NMAX][MMAX];

		read_graph(&n, &m, G, &stat, &missed, &min_degree, &max_degree, &maxv);
		validate_graph(G, n, &stat);

		if (stat != SUCCESS) { //read_graph failed
			handle_error(stat, graphs_processed, n, m, 0, missed); 
		}

		graphs_processed++;

		int dom[n];
		int num_dominated[n];
		int num_choice[n];
		int min_dom[n];

		min_size = n;
		delta = 0;

		// There's a vertex which dominates every other vertex. We can just print the solution here.
		// also allows us to avoid the C4 restriction on the blank test.
		if (max_degree == (n - 1) && min_degree != 0){
			min_dom[0] = maxv;
			print_solution(0, 1, min_dom);
			print_solution(1, 1, min_dom);
			continue;
		}


		// Set up the tests that will let us break out of a branch early.

		if (min_degree > 0){

			if (min_degree == 3){ 
				tests[num_tests] = reed_test;
				num_tests++;
			}

			else if (min_degree == 2 && n != 7){ // doesn't always hold for 7 vertices
				tests[num_tests] = blank_test;
				num_tests++;
			}

			else {
				tests[num_tests] = ore_test;
				num_tests++;
			}
		}


		initialize_dom_set(G, n, num_choice, num_dominated, min_dom, &delta);

		printf("%3d %3d\n",graphs_processed, n);
		find_dom_set(G, n, num_dominated, 0, dom, &min_size, min_dom, num_choice, 0, 0, delta, num_tests);
		print_solution(1, min_size, min_dom);

		num_tests = 2;

	} while(1);

	return 0;
}