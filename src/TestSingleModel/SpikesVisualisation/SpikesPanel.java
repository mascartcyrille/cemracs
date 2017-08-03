/*
 * Copyright (C) 2017 maitre
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */
package TestSingleModel.SpikesVisualisation;

import static TestSingleModel.SpikesVisualisation.MainFrame._TIME_SLIDER;
import java.awt.Graphics;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.JPanel;

/**
 *
 * @author cmascart
 */
public class SpikesPanel extends JPanel {
	private static final Path _SPIKES = Paths.get( "Outputs/spikeTimes.csv" );
	private static ArrayList<Double> _spikesTimes;
	private static double _last, dx;

	public SpikesPanel() {
		super();
		try {
			List<String> readAllLines = Files.readAllLines( _SPIKES );
			_spikesTimes = new ArrayList<>( readAllLines.size() );
			for( String t : readAllLines ) {
				_spikesTimes.add( Double.valueOf( t ) );
			}
			_last = _spikesTimes.get( _spikesTimes.size() - 1 );
		} catch( IOException ex ) {
			Logger.getLogger( SpikesPanel.class.getName() ).log( Level.SEVERE, null, ex );
		}
	}

	public int spikes() {
		return _spikesTimes.size() - 1;
	}

	@Override
	public void paintComponent( Graphics g ) {
		dx = _last / ( MainFrame._SCALE_SLIDER.getValue() );
		g.drawLine( 0, getHeight() / 2 - 1, getWidth(), getHeight() / 2 - 1 );	// Draws the time axis.
		int value = _TIME_SLIDER.getValue();
		double min = value * dx, x, max = ( _last < ( value + 1 ) * dx ) ? _last : ( value + 1 ) * dx;
		int indMin = 0, indMax = _spikesTimes.size() - 1, j = ( indMin + indMax ) / 2;
		while( j != 0 && !( _spikesTimes.get( j - 1 ) < min && _spikesTimes.get( j ) >= min ) ) {
			if( _spikesTimes.get( j ) < min ) {
				indMin = j;
			} else {
				indMax = j;
			}
			j = ( indMin + indMax ) / 2;
		}
		while( j < _spikesTimes.size() && _spikesTimes.get( j ) < max ) {
			x = 600 * ( _spikesTimes.get( j ) - min ) / ( max - min );
			g.drawLine( (int) x, getHeight() / 4, (int) x, getHeight() * 3 / 4 );
			++j;
		}
	}
}
