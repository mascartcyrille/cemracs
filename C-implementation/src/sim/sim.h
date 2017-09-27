#ifndef _SIM_H_
#define _SIM_H_

/* PERSONNAL TYPES */
/* A time interval structure, with an upper and a lower bound, used for storing in which interval the spike in searched. */
typedef struct {
	long double lower_bound;	/* Lower bound, shall be lower than upper_bound */
	long double upper_bound;	/* Upper bound, shall be greater than lower_bound */
} time_interval_t;

/* A connection type structure, for knowing which kind of graph and graph reading method will be used.	*/
typedef enum {
	RECONSTRUCTION,	/* Very heavy graphs (very high number of nodes), which cannot enter in memory		*/
	COMPLETE,		/* Fully connected graphs, except self												*/
	FULL,			/* Fully connected graphs, including self											*/
	INDEPENDENT,	/* No connection graphs																*/
	RANDOM			/* Untyped graphs																	*/
} conn_type;

/* PROTOTYPES */
/* Reads the given line and sets the variable accordingly */
void param( const void* var, const char* const line, const unsigned int size, const unsigned char sep );

/* Checks the return value after a dynamic allocation, prints an error if the value is NULL */
void check_alloc( const void* const var, const char * const var_name );

/* Initialises the value of all global variables */
void init(void);

/* Allocates memory for the global variables */
void create(void);

/* Stores the parameters characteristics of the simulation to a binary file */
void save(void);

/* Frees the memory used for the global variables */
void destroy(void);

/* Computes the maxima and cumulative sum of the maxima of the "probability" of spike */
void compute_m_i_s( const time_interval_t* const time_int );

/* The simulation */
void simulate(void);

/* Returns the default interaction weight, given a presynaptic and a postsynaptic neuron */
long double interaction( const unsigned int presynaptic, const unsigned int postsynaptic );

/* Returns the "probability" of spike for given neuron index at current time */
long double probability( const unsigned int neuron );

/* Returns the maximum "probability" of spike for given neuron and given time interval */
long double max_prob( const unsigned int neuron, const time_interval_t* const time_int );

/* Returns which neuron is the spiking one */
unsigned int which_spiking(void);

/* Returns an integer between 0 and given limit, included */
unsigned int get_int( const unsigned int limit );

/* Returns an integer between 0 and given limit, included, and using the rng of index ind */
unsigned int get_int_using_rng( const unsigned int limit, const unsigned int ind_rng );

/* FOR FULL PARALLEL */
/*  */
void* init_rec_fn( void * ind );

/*  */
void* reconstruction_fn( void * ind );

#endif /* _SIM_H_ */
