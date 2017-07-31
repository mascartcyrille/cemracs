/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Util;

/**
 *
 * @author cmascart
 */
@SuppressWarnings( "EqualsAndHashcode" )
public final class Interaction implements Comparable<Interaction> {
	public int _postSynapticNeuron;
	public double _postSynapticInteraction;

	public Interaction( int postSynapticNeuron, double postSynapticInteraction ) {
		_postSynapticNeuron = postSynapticNeuron;
		_postSynapticInteraction = postSynapticInteraction;
	}

	@Override
	public final int compareTo( final Interaction inter ) {
		return ( _postSynapticNeuron < inter._postSynapticNeuron ) ? -1 : ( _postSynapticNeuron == inter._postSynapticNeuron ) ? 0 : 1;
	}

	@Override
	public final boolean equals( final Object inter ) {
		return ( inter instanceof Interaction ) ? ( _postSynapticNeuron == ( (Interaction) inter )._postSynapticNeuron ) : false;
	}

	@Override
	public String toString() {
		return "" + _postSynapticNeuron;
	}
}
