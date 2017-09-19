#ifndef _SIM_H_
#define _SIM_H_

/* PERSONNAL TYPES */
typedef struct {
	long double lower_bound;
	long double upper_bound;
} time_interval_t;

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

/* Returns the default interaction weight */
long double interaction(void);

/* Returns the "probability" of spike for given neuron index at current time */
long double probability( const unsigned int neuron );

/* Returns the maximum "probability" of spike for given neuron and given time interval */
long double max_prob( const unsigned int neuron, const time_interval_t* const time_int );

/* Returns which neuron is the spiking one */
unsigned int which_spiking(void);

/* Returns an integer between 0 and given limit, included */
unsigned int get_int( const unsigned int limit );

#endif /* _SIM_H_ */