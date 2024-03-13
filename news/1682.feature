Enhanced the "_move_parcel_downstream" function in "network_sediment_transporter.py" to improve sediment movement in steep upstream links. 
The update recalculates parcel velocity at each link as it moves, adjusts movement based on remaining time due to the timestep length and velocity, and dynamically updates the transport capacity of links.
