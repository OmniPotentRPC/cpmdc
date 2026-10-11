`cpmdc_adopt_calculator_comm` reads rank and size from the communicator the host passed.
When `ranks_per_calc` is positive that size must match, and the group root broadcasts the calculator index.
A later call replaces the communicator while no session exists.
