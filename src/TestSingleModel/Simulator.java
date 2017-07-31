/*
 * EVENT: Centre d’Eté Mathématique de Recherche Avancée en Calcul Scientifique (CEMRACS)
 * DATE: 2017
 * PROJECT: Network of interacting neurons with random synaptic weights.
 * AUTHOR: C.MASCART
 */
package TestSingleModel;

import Util.SingleSimulator;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.NavigableSet;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.imageio.ImageIO;

/**
 * Main class of the package CEMRACS.TestSingleModel, that launches the
 * simulation then renders the results. The class extends
 * {@link SingleSimulator}, as the system (a network of neurons) is modeled
 * using a single atomic model.
 *
 * @author cmascart
 */
public class Simulator extends SingleSimulator {
	/**
	 * The number of iterations, or true spikes, to be simulated.
	 */
	public static int _simulationIterations;
	/**
	 * Total number of neurons in the system.
	 */
	public static int _neuronNumber;
	/**
	 * The probability of a directed connection i->j between the two arbitrary
	 * neurons i and j (including self).
	 */
	public static double _connectionProbability;
	/**
	 * A boolean stating whether the post-synaptic neuron's potential change
	 * event must be recorded or not.
	 */
	public static final boolean _INFLUENCED = false;
	/**
	 * A boolean stating whether false spikes events must be recorded or not.
	 */
	public static final boolean _FALSESPIKES = false;
	/**
	 * The line separator of the current system. For use for all the package.
	 */
	public static final String _L_S = System.getProperty( "line.separator" );
	/**
	 * The path separator of the file system.
	 */
	public static final String _F_S = System.getProperty( "file.separator" );
	/**
	 * The network to be simulated.
	 */
	private static Network _net;
	/**
	 * A reference to an instance of this simulator;
	 */
	private static Simulator _simulator;
	/**
	 * The total cumulated sizes of all the event trees. For computing the
	 * average size of the event trees.
	 */
	private static double _totalSizes = 0.0;
	/**
	 * the path to the directory that contains the files storing the results
	 * of the simulation.
	 */
	private static final String _RESULT_DIRECTORY_PATH = _F_S + "dev" + _F_S + "shm" + _F_S + "CEMRACS" + _F_S;
	/**
	 * The path, as a String, to the file used to store the results of the
	 * simulation (that is to say the time of spikes).
	 */
	private static final String _RESULT_FILE_PATH = _RESULT_DIRECTORY_PATH + _F_S + "spikeTimes.csv";
	/**
	 * A counter storing the number of true and false spike. They are stored
	 * as real numbers (while obviously being pure integers) for simplifying
	 * the computation of the proportion between the two later on.
	 */
	public static double _trueSpike = 0, _falseSpike = 0;
	/**
	 * The order of magnitude for the measure of the taken by elements of the
	 * simulation.
	 */
	private static final short _TIME_ORDER_OF_MAGNITUDE = 6;
	/**
	 * A set of events sorted in increasing (time) order.
	 */
	private static final TreeSet<Double> _EVENTS_LINE = new TreeSet<>();
	/**
	 * The average time elapsed between two spikes. Equals the time of last
	 * event divided by the number of iterations (or true spikes).
	 */
	private static double _avgdt;
	/**
	 *
	 */
	private static double _avg;

	/**
	 * Main method of the package, as it launches the simulation.
	 *
	 * @param args
	 */
	public static void main( String[] args ) {
		////////////////////////////////////////////////////////////////////
		//	PARSING ARGUMENTS								//
		////////////////////////////////////////////////////////////////////
		switch( args.length ) {
			case 3:
				_neuronNumber = Integer.valueOf( args[ 0 ] );
				_connectionProbability = Double.valueOf( args[ 1 ] );
				_simulationIterations = Integer.valueOf( args[ 2 ] );
				break;
			default:
				System.out.println( "The number of arguments is 3:\n\tSimulator neuronNumber connectionProbability SimulationIterations" );
				System.exit( 0 );
				break;
		}

		////////////////////////////////////////////////////////////////////
		//	SIMULATION									//
		////////////////////////////////////////////////////////////////////
		long beforeTestCreation = System.nanoTime();
		_net = Network.create();					// Instantiates the network.
		long afterTestCreation = System.nanoTime();
		_simulator = new Simulator( _net );				// Instantiates the simulator.
		long afterCoordinatorCreation = System.nanoTime();
		_simulator.initialize();					// Initialisation.
		long afterInitialisation = System.nanoTime();
		_simulator.simulate( _simulationIterations );		// Simulates a number of iterations. The system iterates from real spike to real spike.
		long afterSimulation = System.nanoTime();

		////////////////////////////////////////////////////////////////////
		//	MEMORY USAGE								//
		////////////////////////////////////////////////////////////////////
		int mb = 1024 * 1024;
		Runtime runtime = Runtime.getRuntime();
		System.out.println( "##### Heap utilization statistics [Mibit] #####" );
		System.out.println( "# Used Memory: " + ( runtime.totalMemory() - runtime.freeMemory() ) / mb + " Mibit" );
		System.out.println( "###############################################\n\n" );

		////////////////////////////////////////////////////////////////////
		//	SIMULATION TIME								//
		////////////////////////////////////////////////////////////////////
		long magnitude;
		String unit;
		if( _TIME_ORDER_OF_MAGNITUDE == 9 ) {
			magnitude = (long) Math.pow( 10, 9 );
			unit = "ns";
		} else if( _TIME_ORDER_OF_MAGNITUDE == 6 ) {
			magnitude = (long) Math.pow( 10, 9 );
			unit = "\u00B5s";
		} else if( _TIME_ORDER_OF_MAGNITUDE == 3 ) {
			magnitude = (long) Math.pow( 10, 9 );
			unit = "ms";
		} else {
			magnitude = (long) Math.pow( 10, 9 );
			unit = "s";
		}
		System.out.println( "##### The different elapsed times [" + unit + "] #####" );
		System.out.println( "# Test creation: " + ( afterTestCreation - beforeTestCreation ) / magnitude + " " + unit );
		System.out.println( "# Coordinator creation: " + ( afterCoordinatorCreation - afterTestCreation ) / magnitude + " " + unit );
		System.out.println( "# Initialisation: " + ( afterInitialisation - afterCoordinatorCreation ) / magnitude + " " + unit );
		System.out.println( "# Simulation: " + ( afterSimulation - afterInitialisation ) / magnitude + " " + unit );
		System.out.println( "############################################\n\n" );

		////////////////////////////////////////////////////////////////////
		//	RENDERING RESULTS								//
		////////////////////////////////////////////////////////////////////
		// Printing the size of the event storage map.
		storeResults();
		spikeAnalysis();
		makeSpikeTrains();
		System.out.println( "##### System state #####" );
		System.out.println( "# Time of the system after simulation: " + _simulator.tN() );
		System.out.println( "# Size of the events map: " + _net.eventsMap().size() );
		System.out.println( "# Average size of the events: " + ( _totalSizes / _neuronNumber ) );
		System.out.println( "# Proportion of true spikes: " + ( 100 * _trueSpike / ( _trueSpike + _falseSpike ) ) + "%" );
		System.out.println( "# Average number of spikes through time: " + ( _trueSpike / _simulator.tL() ) + " spikes/s" );
		System.out.println( "# Average number of spikes on " + _avgdt + ": " + _avg );
		System.out.println( "########################" );
	}

	/**
	 * Stores the content of the events tree map in a file.
	 * Also sorts the events by increasing time.
	 */
	private static void storeResults() {
		try {
			Files.createDirectories( Paths.get( _RESULT_DIRECTORY_PATH ) );
			Path resultFilePath = Paths.get( _RESULT_FILE_PATH );
			Files.write( resultFilePath, "".getBytes() );	// Flushes the content of the file away.
			String event;
			String line;
			for( int i = 0; i < _neuronNumber; ++i ) {
				TreeMap<Double, String> events = _net.eventsMap().get( i );
				NavigableSet<Double> navigableKeySet = events.navigableKeySet();
				_totalSizes += events.size();
				line = i + "";
				for( Double time : navigableKeySet ) {
					line += "," + time;
					_EVENTS_LINE.add( time );
				}
				Files.write( resultFilePath, ( line + _L_S ).getBytes(), StandardOpenOption.APPEND );
				line = "";
				for( Double time : navigableKeySet ) {
					event = events.get( time );
					line += "," + event;
				}
				Files.write( resultFilePath, ( line + _L_S ).getBytes(), StandardOpenOption.APPEND );
			}
		} catch( IOException ex ) {
			Logger.getLogger( Simulator.class.getName() ).log( Level.SEVERE, null, ex );
		}
	}

	/**
	 * Computes analytics about the spikes, essentially rates.
	 * Computes the average time elapsing between two true spikes
	 * {@link _avgdt}, then the average number of spikes in an interval of
	 * length {@link _avgdt}.
	 */
	private static void spikeAnalysis() {
		_avgdt = _simulator.tL() / ( _simulationIterations + 1 );	// There is one true spike (exactly) per simulation iteration, plus the one of the initialization.
		ArrayList<Integer> avgSpike = new ArrayList<>();		// Each node i contains the number of spikes on interval [_avgdt * (i - 1), _avgdt * i].
		avgSpike.add( 0 );
		int nbSpike = 0;
		double avg = 0.0;
		int i = 1;
		for( Double time : _EVENTS_LINE ) {
			if( time <= _avgdt * i ) {					// There is a spike in interval [_avgdt * (i - 1), _avgdt * i].
				++nbSpike;
				avgSpike.set( avgSpike.size() - 1, avgSpike.get( avgSpike.size() - 1 ) );
			} else {								// There are no more spikes in previous interval, a new interval for this spike must be found.
				avg += ( nbSpike / _avgdt );
				do {
					++i;
				} while( time > _avgdt * i );
				nbSpike = 0;
			}
		}
		_avg = avg / ( i );
	}

	/**
	 * Creates a gray scale image (PNG) of the spikes ordered by increasing
	 * time on the real line.
	 */
	private static void makeSpikeTrains() {
		int imageWidth = (int) ( _simulator.tN() * 10000000 ) + 1, imageHeight = 100;
		BufferedImage bi = new BufferedImage( imageWidth, imageHeight, BufferedImage.TYPE_BYTE_BINARY );

		for( int i = 0; i < imageWidth; ++i ) {
			for( int j = 0; j < imageHeight; ++j ) {
				bi.setRGB( i, j, -1 );
			}
		}

		for( int i = 0; i < imageWidth; ++i ) {
			for( int j = 0; j < 1; ++j ) {
				bi.setRGB( i, j + imageHeight / 2 - 1, 1 );
			}
		}
		for( Double eventTime : _EVENTS_LINE ) {
			for( int j = 0; j < imageHeight / 2; ++j ) {
				bi.setRGB( (int) ( eventTime * 10000000 ), j + imageHeight / 4, 1 );
			}
		}
		File f = new File( "fily.png" );
		try {
			ImageIO.write( bi, "PNG", f );
		} catch( IOException ex ) {
			Logger.getLogger( Simulator.class.getName() ).log( Level.SEVERE, null, ex );
		}
	}

	//////////////////////////////////////////////////////////////////////////
	//	CONSTRUCTOR										//
	//////////////////////////////////////////////////////////////////////////
	/**
	 * Simple constructor that calls the {@link #super} one.
	 *
	 * @param network The model to simulate.
	 */
	public Simulator( Network network ) {
		super( network );
	}
}
