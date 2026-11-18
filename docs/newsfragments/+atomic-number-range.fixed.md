The evaluation calls of the C library return a failed `CPMDCResult` with the message `atomic numbers must be in 1..118` for an atomic number outside 1..118.
A session used to accept such a number on its first step, because the range was only checked when a later step was compared with the first, and then refused the identical step with `topology change requires a new session`.
The check runs before a session fixes its topology, so a refused step leaves the session usable.
