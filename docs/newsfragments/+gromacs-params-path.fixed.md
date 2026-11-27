`qmi_input` writes a PARAMS path that survives GROMACS replacing every whole word equal to the `COORD_FILE_NAME`
value and then printing the stored input as a printf format.
A percent sign in the path is written doubled, so that printf restores the path
the library was given.
A path that contains a space, a tab, a hash, an exclamation mark, or a double
quote is one double-quoted token, with a backslash escaping a backslash or a
quote, and a newline is refused.
The coordinate file name stays `gromacs.pdb` unless that whole word occurs in
the written token, in which case the writer picks the next name that does not
occur.
`CPMDCCalculator` treats a one-component magnetic moment column as one moment
per atom.
A wider array with a nonzero entry still raises `ValueError`.
