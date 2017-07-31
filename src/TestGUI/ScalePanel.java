/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package TestGUI;

import static TestGUI.myGUI.potential;
import static TestGUI.myPanel._highestPotential;
import static TestGUI.myPanel._lowestPotential;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Graphics;
import java.util.ArrayList;
import javax.swing.JPanel;

/**
 *
 * @author cmascart
 */
public final class ScalePanel extends JPanel {
	public static final int _NUMBER_OF_STEPS = 9;
	private static final Colors _COLORS = new Colors();

	public ScalePanel() {
		_COLORS.setPreferredSize( new Dimension( 150, 150 ) );
		add( _COLORS );
	}

	static class Colors extends JPanel {
		private static final ArrayList<Color> _COLOR_BY_THRESHOLD = new ArrayList<>();

		public Colors() {
			int i = 0;
			double step = ( _highestPotential - _lowestPotential ) / _NUMBER_OF_STEPS;
			_COLOR_BY_THRESHOLD.add( Color.darkGray );
			_COLOR_BY_THRESHOLD.add( Color.lightGray );
			_COLOR_BY_THRESHOLD.add( Color.red );
			_COLOR_BY_THRESHOLD.add( Color.magenta );
			_COLOR_BY_THRESHOLD.add( Color.orange );
			_COLOR_BY_THRESHOLD.add( Color.yellow );
			_COLOR_BY_THRESHOLD.add( Color.white );
			_COLOR_BY_THRESHOLD.add( Color.green );
			_COLOR_BY_THRESHOLD.add( Color.blue );
			_COLOR_BY_THRESHOLD.add( Color.cyan );
		}

		@Override
		public void paintComponent( Graphics g ) {
			g.setColor( Color.white );
			g.fillRect( 0, 0, getWidth(), getHeight() );
			if( potential ) {
				double step = ( _highestPotential - _lowestPotential ) / _NUMBER_OF_STEPS;
				int offset = ( getHeight() - _COLOR_BY_THRESHOLD.size() * 51 ) / 2;
				int j = 0, i = 0;
				for( Color color : _COLOR_BY_THRESHOLD ) {
					g.setColor( Color.black );
					g.drawString( "Potential >= " + ( _highestPotential - ( step * i++ ) ), 60, j * 50 + 30 + offset + j );
					g.setColor( color );
					g.fillRect( 5, j * 50 + offset + j++, 50, 50 );
				}
			} else {
				int offset = ( getHeight() - 3 * 51 ) / 2;
				g.setColor( Color.black );
				g.drawString( "Pre-synaptic spiking neuron", 60, 30 + offset );
				g.fillRect( 5, offset, 50, 50 );

				g.setColor( Color.black );
				g.drawString( "Post-synaptic neuron", 60, 50 + 30 + offset + 1 );
				g.setColor( Color.darkGray );
				g.fillRect( 5, 50 + offset + 1, 50, 50 );

				g.setColor( Color.yellow );
				g.drawString( "Non activated neuron", 60, 100 + 30 + offset + 2 );
				g.setColor( Color.yellow );
				g.fillRect( 5, 100 + offset + 2, 50, 50 );
			}
		}
	}

	@Override
	public final Component add( Component c ) {
		return super.add( c );
	}
}
