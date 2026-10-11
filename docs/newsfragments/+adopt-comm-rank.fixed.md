`cpmdc_adopt_calculator_comm` reads rank and size from the communicator the host passed.
When `ranks_per_calc` is positive that size must match, and the index is the lowest world rank of one aligned block divided by that size.
A later call replaces the communicator while no session exists.
