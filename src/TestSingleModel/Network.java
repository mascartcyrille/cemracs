/*
 * EVENT: Centre d’Eté Mathématique de Recherche Avancée en Calcul Scientifique (CEMRACS)
 * DATE: 2017
 * PROJECT: Network of interacting neurons with random synaptic weights.
 * AUTHOR: C.MASCART
 */
package TestSingleModel;

import static TestSingleModel.Simulator._FALSESPIKES;
import static TestSingleModel.Simulator._INFLUENCED;
import static TestSingleModel.Simulator._connectionProbability;
import static TestSingleModel.Simulator._falseSpike;
import static TestSingleModel.Simulator._neuronNumber;
import static TestSingleModel.Simulator._trueSpike;
import Util.Interaction;
import Util.SingleSimulator;
import java.util.TreeMap;
import java.util.TreeSet;
import umontreal.iro.lecuyer.randvar.*;
import umontreal.iro.lecuyer.rng.MRG32k3a;

/**
 * This is an atomic model for simulating large networks of neurons modeled
 * using a partial differential equation (PDE, that are differential equations
 * containing not entirely defined functions) derived from an integrate and fire
 * model.
 *
 * @author cmascart
 */
public final class Network {
	/**
	 * Constrains the intensity of the interaction between two neurons. The
	 * intensity is indeed {@link _ALPHA}/{@link _neuronNumber}.
	 */
	private static final double _ALPHA = 1;
	/**
	 * An array of containing all integers from 0 to {@link _neuronNumber} -
	 * 1. It is used to efficiently generate the graph of interactions.
	 */
	private static final int[] _NEURONS_IDS = new int[ _neuronNumber ];
	/**
	 * The neuron spiking at that time. Initialized to - so that no neuron can
	 * be taken as spiking neuron at the beginning of the simulation.
	 */
	private static int _spikingNeuron = -1;
	/**
	 * The base intensity (i.e. intensity at time t = 0) of every neuron.
	 */
	private static final double[] _BASE_POTENTIAL = new double[ _neuronNumber ];
	/**
	 * The Lambda_max value of every neuron in interval [t_k, t_{k+1}].
	 */
	private static final double[] _MAX = new double[ _neuronNumber ];
	/**
	 * The sigma value (for Brownian motion) for every neuron.
	 */
	private static final double[] _SIG = new double[ _neuronNumber ];
	/**
	 * The intensity value at current time for every neuron.
	 */
	private static final double[] _POTENTIAL = new double[ _neuronNumber ];
	/**
	 * The relevant interactions between neurons are stored there. By
	 * "relevant interaction" one must understand an interaction that may
	 * affect the state of the destination neuron whenever the source neuron
	 * is spiking. The information is stored in a dictionary fashion, that is
	 * to stay [sourceNeuron -> (dest_1, dest_2, ..., dest_n)]. The
	 * information "sourceNeuron" and "destinationNeuron" are integers,
	 * corresponding to the neuron id (from 0 to #_neuronNumber - 1).
	 */
	private static final TreeMap<Integer, TreeSet<Interaction>> _INTERACTIONS = new TreeMap<>();
	/**
	 * The coefficient of the linear function b:R->R, b(y)=-lambda*(y-a).
	 */
	private static final double[] _LAMBDA = new double[ _neuronNumber ];
	/**
	 * The delay of the linear function b:R->R, b(y)=-lambda*(y-a).
	 */
	private static final double[] _A = new double[ _neuronNumber ];
	/**
	 * The last spiking times of all neurons.
	 */
	private static final double[] _LAST_SPIKING_TIMES = new double[ _neuronNumber ];
	/**
	 * The random number generator, from which are derived random number
	 * generators fitting different distributions. It is the same for all the
	 * generators used later.
	 */
	private static final MRG32k3a _MRG = new MRG32k3a( "Network" );
	/**
	 * A uniform generator, giving numbers in [0,1].
	 */
	private static final UniformGen _UNIFORM_GEN = new UniformGen( _MRG, 0, 1 );
	/**
	 * A binomial generator,
	 * X~B({@link #_neuronNumber}, {@link #_CONNECTION_PROBABILITY}).
	 */
	private static final BinomialGen _BINOMIAL_GEN = new BinomialGen( _MRG, _neuronNumber, _connectionProbability );
	/**
	 * A tree map to store events (that is spiking) as they appear in the
	 * system.
	 */
	public static final TreeMap<Integer, TreeMap<Double, String>> _eventsMap = new TreeMap<>();
	/**
	 * Time interval between the last event that happen and and the next
	 * internal event occurrence date. Or the time
	 * from now to the next internal event occurrence date.
	 */
	protected double _timeAdvance;
	/**
	 * A reference to the simulator.
	 */
	protected SingleSimulator _simulator;
	/**
	 *
	 */
	public static double endTime = System.nanoTime();

	//////////////////////////////////////////////////////////////////////////
	// FACTORY AND CONSTRUCTORS								//
	//////////////////////////////////////////////////////////////////////////
	/**
	 * A factory that creates a new network. It also initializes the static
	 * attributes that have been declared, instantiated but not filled yet.
	 * The static arrays are of the size of the number of neurons to be
	 * simulated, thus this method can take a bit of time to finish for a high
	 * number of neurons.
	 *
	 * @return The atomic model that has been created and initialized.
	 */
	public static Network create() {

		Network network = new Network();

		for( int id = 0; id < _neuronNumber; ++id ) {
			_NEURONS_IDS[ id ] = id;
		}
		for( int i = 0; i < _neuronNumber; ++i ) {
			_BASE_POTENTIAL[ i ] = 1.1;
			_MAX[ i ] = ProbSpike.max( _BASE_POTENTIAL[ i ] );
			_SIG[ i ] = 0.3;
			_POTENTIAL[ i ] = _BASE_POTENTIAL[ i ];	// Setting intensity at t = 0.
			_LAST_SPIKING_TIMES[ i ] = 0.0;		// No spiking before the begining of the simulation.
			generateInteractions( i );			// The interactions are generated here.
			/*
			 * The b function arguments are set here.
			 */
			_LAMBDA[ i ] = 50.0;
			_A[ i ] = 10;
			_eventsMap.put( i, new TreeMap<>() );
		}

		return network;
	}

	/**
	 * Generating the graph of interactions. The method considers the
	 * interaction graph is quite sparse.
	 *
	 * @param neuron the neuron that is the source of the interaction.
	 */
	private static void generateInteractions( final int preSynapticNeuron ) {
		int max = _neuronNumber - 1, nbCouplings = _BINOMIAL_GEN.nextInt(), influencee, id;
		_INTERACTIONS.put( preSynapticNeuron, new TreeSet<>() );	// At creation, the table of interactions is empty, and every neuron must have a list of interactions, may be empty.
		TreeSet<Interaction> influencees = _INTERACTIONS.get( preSynapticNeuron );
		for( int postSynapticNeuron = 0; postSynapticNeuron < nbCouplings; ++postSynapticNeuron ) {
			id = ( new UniformIntGen( _MRG, 0, max ) ).nextInt();
			influencee = _NEURONS_IDS[ id ];
			_NEURONS_IDS[ id ] = _NEURONS_IDS[ max ];
			_NEURONS_IDS[ max ] = influencee;
			--max;
			influencees.add( new Interaction( influencee, interaction( preSynapticNeuron, postSynapticNeuron ) ) );
		}
	}

	/**
	 * Computes the intensity of the interaction from a given pre-synaptic
	 * neuron to a post-synaptic neuron. Self-interactions are allowed.
	 *
	 * @param preSynapticNeuron  The spiking neuron.
	 * @param postSynapticNeuron The receiving neuron.
	 * @return {@link #_ALPHA}/{@link #_neuronNumber}.
	 */
	private static double interaction( int preSynapticNeuron, int postSynapticNeuron ) {
		return _ALPHA / _neuronNumber;
	}

	/**
	 * A simple, empty default constructor.
	 */
	protected Network() {
	}

	//////////////////////////////////////////////////////////////////////////
	//	GETTERS & SETTERS									//
	//////////////////////////////////////////////////////////////////////////
	/**
	 * Sets the simulator of this model.
	 *
	 * @param simulator
	 */
	public void simulator( final SingleSimulator simulator ) {
		_simulator = simulator;
	}

	/**
	 * A getter for the array of intensities.
	 *
	 * @return The array of intensities, that is to say a reference to the
	 *         array.
	 */
	public double[] theInts() {
		return _POTENTIAL;
	}

	/**
	 * A getter for the post synaptic interaction set of the current spiking
	 * neuron. Returns null if no interaction tree has been found.
	 *
	 * @return A set of the neurons influenced by the current spiking neuron,
	 *         or null if the {@link #_spikingNeuron} number is not valid.
	 */
	public TreeSet<Interaction> influenced() {
		return ( _spikingNeuron >= 0 && _spikingNeuron < _neuronNumber ) ? _INTERACTIONS.get( _spikingNeuron ) : null;
	}

	/**
	 * A getter for the current spiking neuron.
	 *
	 * @return The current spiking neuron.
	 */
	public int spiking() {
		return _spikingNeuron;
	}

	/**
	 * A getter for the events map.
	 *
	 * @return The events map.
	 */
	public TreeMap<Integer, TreeMap<Double, String>> eventsMap() {
		return _eventsMap;
	}

	/**
	 * Gets the time from the last event to the next one.
	 *
	 * @return
	 */
	public double timeAdvance() {
		return _timeAdvance;
	}

	//////////////////////////////////////////////////////////////////////////
	// Simulation										//
	//////////////////////////////////////////////////////////////////////////
	/**
	 * Computes the time to the next real spike of the system. If a false
	 * spike is generated, then generation of the time to the next event is
	 * reiterated, adding to the time already generated. Also updates the
	 * state of the spiking and influenced neurons.
	 * <p>
	 */
	public void nextEvent() {
		double timeOfSpike = 0.0;
		double sumMax = 0.0;
		for( int i = 0; i < _neuronNumber; ++i ) { // Computing the sum, for all neurons, of the maximum value of the probability of spike function, over an interval [ , ].
			sumMax += _MAX[ i ];
		}
		do {
			timeOfSpike = nextTime( sumMax );	// Generates the next time of spike (can be true or false spike at this stage).
			spikingNeuron( sumMax );		// Determines which neuron is actually spiking.
		} while( !isRealSpike( timeOfSpike ) );	// Loops until a real spike is generated.

		_timeAdvance = timeOfSpike;
	}

	/**
	 * Generates the time of the next spike of the system. The spike thus
	 * generated can be either true or false, it has to be determined latter.
	 * The time is a random number, generated according to an exponential
	 * distribution of parameter sumMax.
	 *
	 * @param sumMax The sum of the maximum values of the probability of spike
	 *               function, over an interval [ , ].
	 * @return The time generated.
	 */
	private double nextTime( final double sumMax ) {
		return ( new ExponentialGen( _MRG, sumMax ) ).nextDouble();
	}

	/**
	 * Computes which of the neurons is actually spiking in the system. Each
	 * neuron has a probability of spiking P(I = 1) = lambda_{max}^i /
	 * Sum_{all neurons} lambda_{max}.
	 *
	 * @param sumMax The sum of the maximum values of the probability of spike
	 *               function, over an interval [ , ].
	 */
	private void spikingNeuron( final double sumMax ) {
		double sum = 0;
		double u = _UNIFORM_GEN.nextDouble();
		for( int i = 0; i < _neuronNumber; ++i ) {
			sum += _MAX[ i ];
			if( ( u * sumMax ) <= sum ) {
				_spikingNeuron = i;
				break;
			}
		}
	}

	/**
	 * Determines whether the spike previously generated by the previously
	 * determined neuron is actually a real spike or a false one. This method
	 * also updates the state of the spiking neuron, regardless of the
	 * "trueness" of the spike. If the spike is a true one however it also
	 * updates the state of the neurons under the influence of the spiking
	 * neuron.
	 *
	 * @param timeOfSpike The time elapsed since the last true spike of the
	 *                    system.
	 * @return A boolean stating whether the spike is a true one or not
	 *         ("true" for a real spike, "false" for a false one). The
	 *         {@link #nextEvent()} method will generate a new spike if the
	 *         result is "false".
	 */
	private boolean isRealSpike( final double timeOfSpike ) {
		UniformGen ug = new UniformGen( _MRG, 0, 1 );
		double u = ug.nextDouble();
		updateSpikingState( timeOfSpike );	// The potential of the spiking neuron changes, whatever the flavor of the spike.
		if( ProbSpike.prob( _POTENTIAL[ _spikingNeuron ] ) > ProbSpike.max( _POTENTIAL[ _spikingNeuron ] ) ) {
			System.out.println( "The maximum value for f(V_t) has been exceeded." );
			System.exit( -1 );
		}
		boolean isRealSpike = ( u <= ProbSpike.prob( _POTENTIAL[ _spikingNeuron ] ) / _MAX[ _spikingNeuron ] );
		if( _FALSESPIKES || isRealSpike ) {			// If false spikes must be recorded or the spike is a ture one, record the spike.
			_eventsMap.get( _spikingNeuron ).put( timeOfSpike + _simulator.tN(), "spiking-" + isRealSpike );
		}
		if( isRealSpike ) {
			++_trueSpike;					// Increment the number of true spikes.
			updateInfluencedState( timeOfSpike );	// If the spike is real, the neurons influenced by the spiking one will see their potential change.
		} else {
			++_falseSpike;					// Increment the number of false spikes.
		}

		return isRealSpike;
	}

	/**
	 * Update the value of the potential for the spiking neuron. The current
	 * spike can be real or false, the update algorithm still is the same in
	 * both cases.
	 *
	 * @param timeOfSpike The time elapsed since the last real spike of the
	 *                    system.
	 */
	private void updateSpikingState( final double timeOfSpike ) {
		double nextTime = _simulator.tN() + timeOfSpike;
		double elapsedTime = nextTime - _LAST_SPIKING_TIMES[ _spikingNeuron ];	// Computes the time elapsed since the last spike of the neuron (true or false spike).
		// The variance of the Brownian motion
		double variance = _SIG[ _spikingNeuron ] * _SIG[ _spikingNeuron ]
					* ( 1 - Math.exp( -2 * _LAMBDA[ _spikingNeuron ] * elapsedTime ) / ( 2 * _LAMBDA[ _spikingNeuron ] ) );

		_LAST_SPIKING_TIMES[ _spikingNeuron ] = nextTime;				// Makes the current spike the last spike of the neuron.
		NormalGen ng = new NormalGen( _MRG, 0, variance );				// RNG of normal distribution.

		_POTENTIAL[ _spikingNeuron ] = _A[ _spikingNeuron ]
							 + Math.exp( -_LAMBDA[ _spikingNeuron ] * elapsedTime ) * ( _POTENTIAL[ _spikingNeuron ] - _A[ _spikingNeuron ] ) // The primitive of the b function. The reset is made here.
							 + ng.nextDouble();	// The random normal value generated for the brownian motion.
		_MAX[ _spikingNeuron ] = ProbSpike.max( _POTENTIAL[ _spikingNeuron ] );
	}

	/**
	 * Update the value of the potential for all neurons influenced by the
	 * current spiking neuron. This method is only triggered in a real spike
	 * case.
	 *
	 * @param timeOfSpike The time elapsed since the last real spike of the
	 *                    system.
	 */
	private void updateInfluencedState( final double timeOfSpike ) {
		double nextTime = _simulator.tN() + timeOfSpike;
		_INTERACTIONS.get( _spikingNeuron ).stream().forEach( ( postSynInteraction ) -> {
			int neuron = postSynInteraction._postSynapticNeuron;
			if( neuron == _spikingNeuron ) {
				if( _INFLUENCED ) {
					_eventsMap.get( neuron ).put( timeOfSpike + _simulator.tN(), _eventsMap.get( neuron ).get( timeOfSpike + _simulator.tN() ) + "-influenced" );
				}
				_POTENTIAL[ _spikingNeuron ] += postSynInteraction._postSynapticInteraction;
				_MAX[ _spikingNeuron ] = ProbSpike.max( _POTENTIAL[ _spikingNeuron ] );
			} else {
				if( _INFLUENCED ) {
					_eventsMap.get( neuron ).put( timeOfSpike + _simulator.tN(), "influenced" );
				}
				double elapsedTime = nextTime - _LAST_SPIKING_TIMES[ neuron ];	// Computes the time elapsed since the last spike of the neuron (true or false spike).
				_LAST_SPIKING_TIMES[ neuron ] = nextTime;					// Makes the current spike the last spike of the neuron.
				double variance = _SIG[ neuron ] * _SIG[ neuron ] * ( 1 - Math.exp( -2 * _LAMBDA[ neuron ] * elapsedTime ) / ( 2 * _LAMBDA[ neuron ] ) ); // The variance for the brownian motion.
				NormalGen ng = new NormalGen( _MRG, 0, variance );
				_POTENTIAL[ neuron ] = _A[ neuron ] + Math.exp( -_LAMBDA[ neuron ] * elapsedTime ) * ( _POTENTIAL[ neuron ] - _A[ neuron ] ) // The primitive of the b function. The reset is made here.
							     + ng.nextDouble() // The random normal value generated for the brownian motion.
							     + postSynInteraction._postSynapticInteraction;	// The interaction between the spiking neuron and this neuron.
				_MAX[ neuron ] = ProbSpike.max( _POTENTIAL[ neuron ] );
			}
		} );
	}
}
