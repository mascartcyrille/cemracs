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

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JSlider;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

/**
 *
 * @author cmascart
 */
public class MainFrame extends JFrame {
	public static final double _SCALE = 100;
	private static SpikesPanel _SPIKES;
	public static JSlider _TIME_SLIDER;
	public static JSlider _SCALE_SLIDER;
	private static final JPanel _MAINPAN = new JPanel();

	public static void main( String[] args ) {
		MainFrame window = new MainFrame();
		window.initialize();
		window.start();
	}

	public MainFrame() {
		_SPIKES = new SpikesPanel();
		_SPIKES.setPreferredSize( new Dimension( 600, 500 ) );

		_SCALE_SLIDER = new JSlider( JSlider.VERTICAL, 1, 100, 10 );
		_SCALE_SLIDER.setPreferredSize( new Dimension( 50, 400 ) );
		_SCALE_SLIDER.addChangeListener( new ScaleListener() );

		_TIME_SLIDER = new JSlider( 0, _SCALE_SLIDER.getValue() - 1, 0 );
		_TIME_SLIDER.setPreferredSize( new Dimension( 500, 50 ) );
		_TIME_SLIDER.setMinorTickSpacing( 1 );
		_TIME_SLIDER.setPaintTicks( true );

		_MAINPAN.add( _SCALE_SLIDER, BorderLayout.WEST );
		_MAINPAN.add( _SPIKES, BorderLayout.CENTER );
		_MAINPAN.add( _TIME_SLIDER, BorderLayout.SOUTH );
	}

	private void initialize() {
		setTitle( "Title" );
		// Window size
		setSize( 800, 600 );
		// Set window to center
		setLocationRelativeTo( null );
		// Ends the processus
		setDefaultCloseOperation( JFrame.EXIT_ON_CLOSE );
		// Makes size static
		setResizable( true );
		// Keeps window on top
		setAlwaysOnTop( true );
		// Removes decorations
		setUndecorated( false );

		setContentPane( _MAINPAN );

		// Set window visible
		setVisible( true );

		_MAINPAN.repaint();
	}

	@SuppressWarnings( "SleepWhileInLoop" )
	private void start() {
		while( true ) {
			_MAINPAN.repaint();		// Repainting of the frame for acknowledging the changes.
			try {
				Thread.sleep( 500 );	// A few milliseconds sleep for the repainting
			} catch( InterruptedException e ) {
				Logger.getLogger( MainFrame.class.getName() ).log( Level.SEVERE, "The main thread has been interrupted during its sleeping in method myGUI.start().", e );
			}
		}
	}

	class ScaleListener implements ChangeListener {
		public ScaleListener() {
		}

		@Override
		public void stateChanged( ChangeEvent e ) {
			_TIME_SLIDER.setMaximum( _SCALE_SLIDER.getValue() - 1 );
			_TIME_SLIDER.updateUI();
			_TIME_SLIDER.repaint();
			_MAINPAN.repaint();
		}

	}
}
