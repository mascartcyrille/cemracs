/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package TestGUI;

import static TestGUI.ScalePanel._NUMBER_OF_STEPS;
import static TestGUI.myGUI.potential;
import TestSingleModel.Simulator;
import static TestSingleModel.Simulator._NEURON_NUMBER;
import Util.Interaction;
import java.awt.Color;
import java.awt.Graphics;
import javax.swing.JPanel;

/**
 *
 * @author cmascart
 */
public class myPanel extends JPanel {
	public static double _highestPotential = 0.0;
	public static double _lowestPotential = 0.0;
	private Simulator _sim;
	private static final int _PIXEL_WIDTH = (int) ( 400 / Math.sqrt( _NEURON_NUMBER ) );
	private static final int _PIXEL_HEIGHT = (int) ( 400 / Math.sqrt( _NEURON_NUMBER ) );
	private static int _OFFSET_X = 10;
	private static int _OFFSET_Y = 10;

	public myPanel() {
	}

	public void setSim( Simulator sim ) {
		_sim = sim;
	}

	@Override
	public void paintComponent( Graphics g ) {
		g.clearRect( 0, 0, getWidth(), getHeight() );
		g.setColor( Color.white );
		g.fillRect( 0, 0, getWidth(), getHeight() );
		_OFFSET_X = (int) ( ( getWidth() - ( Math.sqrt( _NEURON_NUMBER ) * _PIXEL_WIDTH ) ) / 2 );
		_OFFSET_Y = (int) ( getHeight() - ( Math.sqrt( _NEURON_NUMBER ) * _PIXEL_HEIGHT ) ) / 2;
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
		g.setColor( Color.black );
		g.drawString( "Current time of simulation: " + _sim.tL(), _OFFSET_X, (int) ( _OFFSET_Y + Math.sqrt( _NEURON_NUMBER ) * _PIXEL_HEIGHT ) + 10 );
		g.drawString( "Highest potential: " + _highestPotential, _OFFSET_X, (int) ( _OFFSET_Y + Math.sqrt( _NEURON_NUMBER ) * _PIXEL_HEIGHT ) + 25 );
		g.drawString( "Lowest potential: " + _lowestPotential, _OFFSET_X, (int) ( _OFFSET_Y + Math.sqrt( _NEURON_NUMBER ) * _PIXEL_HEIGHT ) + 40 );
	}
}
