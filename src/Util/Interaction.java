/*
 * EVENT: Centre d’Eté Mathématique de Recherche Avancée en Calcul Scientifique (CEMRACS)
 * DATE: 2017
 * PROJECT: Network of interacting neurons with random synaptic weights.
 * AUTHOR: C.MASCART
 */
package Util;

/**
 * A ++C equivalent of a pair, which instances can be compared to one another.
 * Interactions are compared by their {@link _postSynapticNeuron} attribute,
 * that is the id of the post-synaptic neuron in an interaction.
 * This class stores the interaction between a given pre-synaptic neuron and one
 * of its post-synaptic fellow neuron, that is to say the identity of the
 * post-synaptic neuron and the potential change of the interaction.
 *
 * @author cmascart
 */
@SuppressWarnings( "EqualsAndHashcode" )
public final class Interaction implements Comparable<Interaction> {
	/**
	 * The identity (number) of the post-synaptic neuron in this interaction.
	 */
	public int _postSynapticNeuron;
	/**
	 * The potential change in this interaction.
	 */
	public double _postSynapticInteraction;

	/**
	 * A simple constructor that initializes the two attributes.
	 *
	 * @param postSynapticNeuron      An identity value.
	 * @param postSynapticInteraction A potential.
	 */
	public Interaction( int postSynapticNeuron, double postSynapticInteraction ) {
		_postSynapticNeuron = postSynapticNeuron;
		_postSynapticInteraction = postSynapticInteraction;
	}

	/**
	 * Used as a comparison operator (&lt;, &gt;, =, etc.).
	 *
	 * @param inter Another interaction for comparison.
	 * @return a negative integer, zero, or a positive integer as this object
	 *         is less than, equal to, or greater than the specified object.
	 */
	@Override
	public final int compareTo( final Interaction inter ) {
		return ( _postSynapticNeuron < inter._postSynapticNeuron ) ? -1 : ( _postSynapticNeuron == inter._postSynapticNeuron ) ? 0 : 1;
	}

	/**
	 * Equality operator.
	 *
	 * @param inter An object for equality testing.
	 * @return true if the object is of type Interaction and its
	 *         {@link _postSynapticNeuron} attribute is the same as this
	 *         object.
	 */
	@Override
	public final boolean equals( final Object inter ) {
		return ( inter instanceof Interaction ) ? ( _postSynapticNeuron == ( (Interaction) inter )._postSynapticNeuron ) : false;
	}

	/**
	 * Creates and returns a String representation of this instance. The
	 * String corresponds to the String representation of the
	 * {@link _postSynapticNeuron} attribute.
	 *
	 * @return A String representation of this object.
	 */
	@Override
	public String toString() {
		return "" + _postSynapticNeuron;
	}
}
