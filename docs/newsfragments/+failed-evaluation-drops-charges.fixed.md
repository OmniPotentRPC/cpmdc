A failed evaluation drops the Hirshfeld charges stored on the session.
`cpmdc_session_last_charges()` then returns the same not-available status as
before any evaluation.
Energy components, the gradient, and stress from that session are dropped with
the charges.
