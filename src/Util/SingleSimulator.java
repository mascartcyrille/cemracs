/*
 * EVENT: Centre d’Eté Mathématique de Recherche Avancée en Calcul Scientifique (CEMRACS)
 * DATE: 2017
 * PROJECT: Network of interacting neurons with random synaptic weights.
 * AUTHOR: C.MASCART
 */
package Util;

import TestSingleModel.Network;

/**
 * Simulator for only one simple model.
 * This class acts as a skeleton for the simulator, so that all the logic is
 * kept hidden here.
 *
 * @author cmascart
 */
public abstract class SingleSimulator {
	//////////////////////////////////////////////////////////////////////////
	//	ATTRIBUTES										//
	//////////////////////////////////////////////////////////////////////////
	/**
	 * Last event occurrence date.
	 */
	protected Double _tL;
	/**
	 * Next event occurrence date.
	 */
	protected Double _tN;
	/**
	 * Simulation model. Represents a network of interacting neurons.
	 */
	protected Network _network;

	//////////////////////////////////////////////////////////////////////////
	//	CONSTRUCTORS									//
	//////////////////////////////////////////////////////////////////////////
	/**
	 * Specialize constructor that set the name to "SingleSimulator" and take
	 * a model to simulate.
	 *
	 * @param network
	 */
	public SingleSimulator( Network network ) {
		_network = network;
		_tN = _tL = 0.0;
	}

	//////////////////////////////////////////////////////////////////////////
	//	GETTERS & SETTERS									//
	//////////////////////////////////////////////////////////////////////////
	/**
	 * Last event occurrence date.
	 */
	public Double tL() {
		return _tL;
	}

	/**
	 * Next event occurrence date.
	 */
	public Double tN() {
		return _tN;
	}

	/**
	 * Simulation model.
	 */
	public Network network() {
		return _network;
	}

	//////////////////////////////////////////////////////////////////////////
	//	CLASS' DYNAMIC & METHODS							//
	//////////////////////////////////////////////////////////////////////////
	/**
	 * Initialize the simulator and its model, setting the model's simulator
	 * to this and computing the time of occurrence of the first event.
	 */
	public void initialize() {
		_network.simulator( this );
		_network.nextEvent();
		timeAdvance();
	}

	/**
	 * Marks the end of the current iteration, thus advancing the system to
	 * its next event time.
	 */
	public void timeAdvance() {
		_tL = _tN;
		_tN = _tL + _network.timeAdvance();
	}

	/**
	 * Launches simulations, stopping after a number of iterations have
	 * passed.
	 */
	public void simulate( Integer itNb ) {
		int it = 0;
		while( it++ < itNb ) {
			simulate();
		}
	}

	/**
	 * Simulates an event, computing the evolution of the system (time and
	 * potentials).
	 */
	public void simulate() {
		_network.nextEvent();
		timeAdvance();
	}
}
