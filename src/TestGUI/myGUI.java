/*
 * EVENT: Centre d’Eté Mathématique de Recherche Avancée en Calcul Scientifique (CEMRACS)
 * DATE: 2017
 * PROJECT: Network of interacting neurons with random synaptic weights.
 * AUTHOR: C.MASCART
 */
package TestGUI;

import TestSingleModel.Simulator;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JPanel;

/**
 * Creates the GUI, organizing the different parts in layers and initializing
 * and starting the simulation, repainting the GUI at each step.
 *
 * @author cmascart
 */
public class myGUI extends JFrame {
	/**
	 * If true, displays the neurons' potential, else displays the spiking
	 * neurons.
	 */
	public static boolean potential = false;
	/**
	 * The panel is a square representation of the network.
	 */
	private final myPanel _pan;
	/**
	 * The color scale information, displayed at the left of the window.
	 */
	private final JPanel _scale;
	/**
	 * Contains the window relevant content (network representation, buttons
	 * and color scale).
	 */
	private final JPanel _drawing;
	/**
	 * Contains the network representation and the buttons.
	 */
	private final JPanel _simPan;
	/**
	 * Contains the buttons.
	 */
	private final JPanel _buttonPanel;
	/**
	 * A button for pausing the simulation.
	 */
	private final JButton _pauseButton;
	/**
	 * A button for changing the display, either the spiking neurons or their
	 * potential.
	 */
	private final JButton _displayButton;
	/**
	 * The actual simulator.
	 */
	private final Simulator _sim;
	/**
	 * Is the simulation paused or not?
	 */
	boolean pause = true;

	/**
	 * Constructor that instantiates the attributes and orders the components
	 * in several layouts.
	 *
	 * @param simulator The simulator this frame is a representation of.
	 */
	public myGUI( Simulator simulator ) {
		_pan = new myPanel();
		_scale = new ScalePanel();
		_drawing = new JPanel();
		_simPan = new JPanel();
		_buttonPanel = new JPanel();
		_pauseButton = new JButton( "Pause" );
		_displayButton = new JButton( "Show potential" );

		_sim = simulator;
		_pan.setSim( simulator );

		_pauseButton.addActionListener( ( ActionEvent e ) -> {
			pause = !pause;
		} );
		_buttonPanel.add( _pauseButton );

		_displayButton.addActionListener( ( ActionEvent e ) -> {
			potential = !potential;
			if( potential ) {
				_displayButton.setText( "Show spiking neurons" );
			} else {
				_displayButton.setText( "Show potential" );
			}
			_drawing.repaint();
		} );
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

	/**
	 * Initializes the frame, setting the window's title, size, location, and
	 * other options. Finally makes the window visible and paints the main
	 * frame.
	 */
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
	}

	/**
	 * Starts the simulation, repainting the frame as the time advances and
	 * the potentials evolve. If the pause button is pressed, the main thread
	 * goes to an idle state for 500 ms. It loops on this idle state until the
	 * user clicks on the pause button again. The simulation starts in pause.
	 */
	@SuppressWarnings( "SleepWhileInLoop" )
	public void start() {
		while( pause ) {
			try {
				Thread.sleep( 500 );	// Waiting for the user to click on the pause button for the simulation to start.
			} catch( InterruptedException e ) {
				Logger.getLogger( myGUI.class.getName() ).log( Level.SEVERE, "The main thread has been interrupted during its sleeping in method myGUI.start().", e );
			}
		}
		_sim.initialize();			// Initialization of the simulation.
		_drawing.repaint();			// Repainting of the frame for acknoledging the changes.
		try {
			Thread.sleep( 2 );		// A few milliseconds sleep for the repainting
		} catch( InterruptedException e ) {
			Logger.getLogger( myGUI.class.getName() ).log( Level.SEVERE, "The main thread has been interrupted during its sleeping in method myGUI.start().", e );
		}
		for( ;; ) {
			while( pause ) {
				try {
					Thread.sleep( 500 );
				} catch( InterruptedException e ) {
					Logger.getLogger( myGUI.class.getName() ).log( Level.SEVERE, "The main thread has been interrupted during its sleeping in method myGUI.start().", e );
				}
			}
			_sim.simulate();			// One iteration (advancing to the next true spike) of the simulation.
			_drawing.repaint();		// Repainting of the frame for acknoledging the changes.
			try {
				Thread.sleep( 2 );	// A few milliseconds sleep for the repainting
			} catch( InterruptedException e ) {
				Logger.getLogger( myGUI.class.getName() ).log( Level.SEVERE, "The main thread has been interrupted during its sleeping in method myGUI.start().", e );
			}
		}
	}
}
