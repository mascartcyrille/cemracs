/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package TestGUI;

import TestSingleModel.Simulator;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JPanel;

/**
 * Creates the GUI, organizing the different parts in layers and initializing
 * and staring the simulation, repainting the GUI at each step.
 *
 * @author cmascart
 */
public class myGUI extends JFrame {
	myPanel _pan = new myPanel();
	JPanel _scale = new ScalePanel();
	JPanel _drawing = new JPanel();
	JPanel _simPan = new JPanel();
	JPanel _buttonPanel = new JPanel();
	JButton _pauseButton = new JButton( "Pause" );
	JButton _displayButton = new JButton( "Show potential" );
	Simulator _sim;
	boolean pause = true;
	public static boolean potential = false;

	public myGUI( Simulator simulator ) {
		_sim = simulator;
		_pan.setSim( simulator );

		_pauseButton.addActionListener( new PauseListener() );
		_buttonPanel.add( _pauseButton );

		_displayButton.addActionListener( new DisplayListener() );
		_buttonPanel.add( _displayButton );

		_simPan.setLayout( new BorderLayout() );
		_simPan.add( _pan, BorderLayout.CENTER );
		_simPan.add( _buttonPanel, BorderLayout.SOUTH );

		_buttonPanel.setBackground( Color.white );
		_simPan.setBackground( Color.white );
		_drawing.setBackground( Color.white );

		_drawing.setLayout( new BoxLayout( _drawing, BoxLayout.LINE_AXIS ) );
		_simPan.setPreferredSize( new Dimension( 500, 600 ) );
		_drawing.add( _simPan );
		_scale.setLayout( new BoxLayout( _scale, BoxLayout.LINE_AXIS ) );
		_scale.setPreferredSize( new Dimension( 300, 600 ) );
		_drawing.add( _scale );
	}

	public void initialize() {
		////////////////////////////////////////////////////////////////////
		// Base window									//
		////////////////////////////////////////////////////////////////////
		// Window title
		setTitle( "HeatMap" );
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

		setContentPane( _drawing );

		// Set window visible
		setVisible( true );

		_drawing.repaint();
		try {
			Thread.sleep( 2 );
		} catch( InterruptedException e ) {
			Logger.getLogger( myGUI.class.getName() ).log( Level.SEVERE, "The main thread has been interrupted during its sleeping in method myGUI.initialize()." );
		}
	}

	@SuppressWarnings( "SleepWhileInLoop" )
	public void start() {
		while( pause ) {
			try {
				Thread.sleep( 500 );
			} catch( InterruptedException e ) {
				Logger.getLogger( myGUI.class.getName() ).log( Level.SEVERE, "The main thread has been interrupted during its sleeping in method myGUI.start()." );
			}
		}
		_sim.initialize();	// Initialization of the simulation.
		_drawing.repaint();
		try {
			Thread.sleep( 2 );
		} catch( InterruptedException e ) {
			Logger.getLogger( myGUI.class.getName() ).log( Level.SEVERE, "The main thread has been interrupted during its sleeping in method myGUI.start()." );
		}
		for( ;; ) {
			while( pause ) {
				try {
					Thread.sleep( 500 );
				} catch( InterruptedException e ) {
					Logger.getLogger( myGUI.class.getName() ).log( Level.SEVERE, "The main thread has been interrupted during its sleeping in method myGUI.start()." );
				}
			}
			_sim.simulate();
			_drawing.repaint();
			try {
				Thread.sleep( 2 );
			} catch( InterruptedException e ) {
				Logger.getLogger( myGUI.class.getName() ).log( Level.SEVERE, "The main thread has been interrupted during its sleeping in method myGUI.start()." );
			}
		}
	}

	class PauseListener implements ActionListener {
		@Override
		public void actionPerformed( ActionEvent e ) {
			pause = !pause;
		}

	}

	class DisplayListener implements ActionListener {
		@Override
		public void actionPerformed( ActionEvent e ) {
			potential = !potential;
			if( potential ) {
				_displayButton.setText( "Show spiking neurons" );
			} else {
				_displayButton.setText( "Show potential" );
			}
			_drawing.repaint();
		}

	}
}
