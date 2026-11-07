An opt-in MPI setup failure now returns on every rank together, and an FFT setup failure releases the buffers it already took.
The optional CUDA compile line passes the scaffold include directory, and this change does not run CUDA.
