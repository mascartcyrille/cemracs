/*
 * This is an open source project, implemented as part of a neuroscience project during the 2017 session of CEMRACS in CIRM, Marseilles.
 */
package TestSingleModel;

/**
 * This class aims at giving a skeleton for the probability of spike function.
 * It defines two public methods, {@link #max(double, double, double)} and
 * {@link #prob(double)}, that are used in the {@link Network} class.
 *
 * @author cmascart
 */
public class ProbSpike {
	private static final double _EXPONENT = 2.0;
	private static final double _CONST = 100000;

	public static double max() {
		return Math.pow( 30, _EXPONENT ) * _CONST;
	}

	public static double max( double potential ) {
		double x = ( 1.1 * potential ) - 1;

		return ( x > 0 ) ? ( Math.pow( x, _EXPONENT ) * _CONST ) : 0.0;
	}

	public static double prob( double potential ) {
		double x = potential - 1;

		return ( x > 0 ) ? ( Math.pow( x, _EXPONENT ) * _CONST ) : 0.0;
	}
}
