/**
 * This package intends to simulate a network of neurons modeled as Poisson processes.
 * It is composed of three classes:
 *	- {@link Simulator}: The main class of the package, that launches the
 * simulation and computes the results.
 * - {@link Network}: The network model, what is executed at each iteration
 * turn.
 * - {@link ProbSpike}: The probability of spike function. It deserves its own
 * class as the behavior depends on what information is known by the modeler.
 * <p>
 * The main class of this package is the entry for a command line execution of
 * the simulation. For a more graphical application, see
 * {@link TestGUI.execGUI}.
 */
package TestSingleModel;
