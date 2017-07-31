/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package TestGUI;

import TestSingleModel.Network;
import TestSingleModel.Simulator;

/**
 *
 * @author cmascart
 */
public class execGUI {
	public static void main( String[] args ) {
		Network net = Network.create();	// Instantiates the network.
		Simulator simulator = new Simulator( net );	// Instantiates the simulator.

		myGUI window = new myGUI( simulator );
		window.initialize();					// Initializes the window. Displays the window with system after initialization.
		window.start();						// Starts simulation and painting of the window.
	}
}
