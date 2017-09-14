/*
 * EVENT: Centre d’Eté Mathématique de Recherche Avancée en Calcul Scientifique (CEMRACS)
 * DATE: 2017
 * PROJECT: Network of interacting neurons with random synaptic weights.
 * AUTHOR: C.MASCART
 */
package TestGUI;

import static TestGUI.ScalePanel._NUMBER_OF_STEPS;
import static TestGUI.myGUI.potential;
import TestSingleModel.Simulator;
import static TestSingleModel.Simulator._neuronNumber;
import Util.Interaction;
import java.awt.Color;
import java.awt.Graphics;
import javax.swing.JPanel;

/**
 * Contains a square representation of the network, either as a set of
 * potentials or as spiking, influenced and idle neurons.
 *
 * @author cmascart
 */
public class myPanel extends JPanel {
	/**
	 * The highest potential ever recorded since the beginning of the
	 * simulation.
	 */
	public static double _highestPotential = 0.0;
	/**
	 * The lowest potential ever recorded since the beginning of the
	 * simulation.
	 */
	public static double _lowestPotential = 0.0;
	/**
	 * The simulator.
	 */
	private Simulator _sim;
	/**
	 * The width in pixels of the square representing the simulation.
	 */
	private static final int _PIXEL_WIDTH = (int) ( 400 / Math.sqrt( _neuronNumber ) );
	/**
	 * The height in pixels of the square representing the simulation.
	 */
	private static final int _PIXEL_HEIGHT = (int) ( 400 / Math.sqrt( _neuronNumber ) );
	/**
	 * The horizontal offset for the square representation of the simulation.
	 */
	private static int _OFFSET_X;
	/**
	 * The vertical offset for the square representation of the simulation.
	 */
	private static int _OFFSET_Y;

	/**
	 * A simple, empty constructor.
	 */
	public myPanel() {
	}

	/**
	 * Sets the simulator.
	 *
	 * @param sim
	 */
	public void setSim( Simulator sim ) {
		_sim = sim;
	}

	/**
	 * The painting method, where the network is put as a square and the color
	 * of the points are chosen depending on what information the user wants
	 * to see.
	 *
	 * @param g
	 */
	@Override
	public void paintComponent( Graphics g ) {
		g.setColor( Color.white );
		g.fillRect( 0, 0, getWidth(), getHeight() );
		_OFFSET_X = (int) ( ( getWidth() - ( Math.sqrt( _neuronNumber ) * _PIXEL_WIDTH ) ) / 2 );
		_OFFSET_Y = (int) ( ( getHeight() - ( Math.sqrt( _neuronNumber ) * _PIXEL_HEIGHT ) ) / 2 );
		double[] pots = _sim.network().theInts();
		double side = Math.floor( Math.sqrt( pots.length ) );
		int k = 0;
		for( int i = 0; i < side; ++i ) {
			for( int j = 0; j < side; ++j ) {
				_highestPotential = ( _highestPotential >= pots[ k ] ) ? _highestPotential : pots[ k ];
				_lowestPotential = ( _lowestPotential <= pots[ k ] ) ? _lowestPotential : pots[ k ];
				double step = ( _highestPotential - _lowestPotential ) / _NUMBER_OF_STEPS;
				int t = 0;
				if( potential ) {
					if( pots[ k ] >= _highestPotential ) {
						g.setColor( Color.darkGray );
					} else if( pots[ k ] >= _highestPotential - ( step * ++t ) ) {
						g.setColor( Color.lightGray );
					} else if( pots[ k ] >= _highestPotential - ( step * ++t ) ) {
						g.setColor( Color.red );
					} else if( pots[ k ] >= _highestPotential - ( step * ++t ) ) {
						g.setColor( Color.magenta );
					} else if( pots[ k ] >= _highestPotential - ( step * ++t ) ) {
						g.setColor( Color.orange );
					} else if( pots[ k ] >= _highestPotential - ( step * ++t ) ) {
						g.setColor( Color.yellow );
					} else if( pots[ k ] >= _highestPotential - ( step * ++t ) ) {
						g.setColor( Color.white );
					} else if( pots[ k ] >= _highestPotential - ( step * ++t ) ) {
						g.setColor( Color.green );
					} else if( pots[ k ] >= _highestPotential - ( step * ++t ) ) {
						g.setColor( Color.cyan );
					} else {
						g.setColor( Color.blue );
					}
				} else if( k == _sim.network().spiking() ) {
					g.setColor( Color.black );
				} else if( _sim.network().influenced() != null && !_sim.network().influenced().isEmpty() && _sim.network().influenced().contains( new Interaction( k, 0 ) ) ) {
					g.setColor( Color.darkGray );
				} else {
					g.setColor( Color.yellow );
				}
				++k;
				g.fillRect( _PIXEL_WIDTH * i + _OFFSET_X, _PIXEL_HEIGHT * j + _OFFSET_Y, _PIXEL_WIDTH, _PIXEL_HEIGHT );
			}
		}
		g.setColor( Color.black );	// Sets the text in black.
		g.drawString( "Current time of simulation: " + _sim.tL(), _OFFSET_X, (int) ( _OFFSET_Y + Math.sqrt( _neuronNumber ) * _PIXEL_HEIGHT ) + 10 );
		g.drawString( "Highest potential: " + _highestPotential, _OFFSET_X, (int) ( _OFFSET_Y + Math.sqrt( _neuronNumber ) * _PIXEL_HEIGHT ) + 25 );
		g.drawString( "Lowest potential: " + _lowestPotential, _OFFSET_X, (int) ( _OFFSET_Y + Math.sqrt( _neuronNumber ) * _PIXEL_HEIGHT ) + 40 );
	}
}
