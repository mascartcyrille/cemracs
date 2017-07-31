/*
 * EVENT: Centre d’Eté Mathématique de Recherche Avancée en Calcul Scientifique (CEMRACS)
 * DATE: 2017
 * PROJECT: Network of interacting neurons with random synaptic weights.
 * AUTHOR: C.MASCART
 */
package TestSingleModel;

/**
 * This class aims at giving a skeleton for the probability of spike function.
 * It defines two public methods, {@link #max(double, double, double)} and
 * {@link #prob(double)}, that are used in the {@link Network} class.
 * <p>
 * This class represents the mathematical function f:R->R+,
 * f(x) = K * (Pos(x - 1))^n, with n = {@link _EXPONENT} and K = {@link _CONT}.
 *
 * @author cmascart
 */
public class ProbSpike {
	/**
	 * The exponent of the function.
	 */
	private static final double _EXPONENT = 2.0;
	/**
	 * The multiplicative constant of the function.
	 */
	private static final double _CONST = 100000;

	/**
	 * Computes a maximum value for the probability of spike, given a neuron's
	 * potential. This maximum is a "soft" one, as it may actually be exceeded
	 * in practice, due to the unboundedness of the Brownian term in the
	 * potential calculus. This "maximum" is computed considering the
	 * potential only has a very small chance of evolving more than a certain
	 * percentage of itself.
	 *
	 * @param potential The potential of the neuron the maximum of the
	 *                  probability of spike must be computed.
	 * @return The maximum value of the probability of spike.
	 */
	public static double max( double potential ) {
		double x = ( 1.1 * potential ) - 1;

		return ( x > 0 ) ? ( Math.pow( x, _EXPONENT ) * _CONST ) : 0.0;
	}

	/**
	 * Computes the actual probability of spike given the value a neuron's
	 * potential. The probability is exact (not counting rounding issues).
	 *
	 * @param potential A neuron's potential.
	 * @return The probability of spike.
	 */
	public static double prob( double potential ) {
		double x = potential - 1;

		return ( x > 0 ) ? ( Math.pow( x, _EXPONENT ) * _CONST ) : 0.0;
	}
}
