/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Util;

import TestSingleModel.Network;

/**
 * Simulator for only one simple model.
 *
 * @author cmascart
 */
public abstract class SingleSimulator {
	// ATTRIBUTES
	/**
	 * Last event occurrence date.
	 */
	protected Double _tL;
	/**
	 * Next event occurrence date.
	 */
	protected Double _tN;
	/**
	 * Simulation model.
	 */
	protected Network _network;

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

	// GETTERS & SETTERS
	/**
	 * {@inheritDoc}
	 */
	public Double tL() {
		return _tL;
	}

	/**
	 * {@inheritDoc}
	 */
	public Double tN() {
		return _tN;
	}

	/**
	 * {@inheritDoc}
	 */
	public Network network() {
		return _network;
	}

	// CLASS' DYNAMIC METHODS
	/**
	 * {@inheritDoc}
	 */
	public void initialize() {
		_network.simulator( this );
		_network.nextEvent();
		timeAdvance();
	}

	/**
	 * {@inheritDoc}
	 */
	public void timeAdvance() {
		_tL = _tN;
		_tN = _tL + _network.timeAdvance();
	}

	/**
	 * {@inheritDoc}
	 */
	public void simulate( Integer itNb ) {
		int it = 0;
		while( it++ < itNb ) {
			simulate();
		}
	}

	/**
	 * {@inheritDoc}
	 */
	public void simulate() {
		_network.nextEvent();
		timeAdvance();
	}
}
