/*
 * EVENT: Centre d’Eté Mathématique de Recherche Avancée en Calcul Scientifique (CEMRACS)
 * DATE: 2017
 * PROJECT: Network of interacting neurons with random synaptic weights.
 * AUTHOR: C.MASCART
 */
package TestGUI;

import TestSingleModel.Network;
import TestSingleModel.Simulator;
import static TestSingleModel.Simulator._connectionProbability;
import static TestSingleModel.Simulator._neuronNumber;
import static TestSingleModel.Simulator._simulationIterations;

/**
 * This class only serves as an interface for creating the simulators and GUIs.
 * The only method is the {@link #main(java.lang.String[])}.
 *
 * @author cmascart
 */
public class execGUI {
	/**
	 * Instantiates the network, then the simulator, then creates the windowed
	 * representation of the network and starts the simulation.
	 *
	 * @param args Unused at the moment.
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

		Network net = Network.create();			// Instantiates the network.
		Simulator simulator = new Simulator( net );	// Instantiates the simulator.

		myGUI window = new myGUI( simulator );
		window.initialize();					// Initializes the window. Displays the window with system after initialization.
		window.start();						// Starts simulation and painting of the window.
	}
}
