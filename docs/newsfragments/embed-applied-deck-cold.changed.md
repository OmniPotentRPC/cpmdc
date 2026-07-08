Cold OpenCPMD SCF prefers Cap'n-rendered `applied_input_deck` when it includes
real `&ATOMS` PP lines. Method-only decks (empty C-render `&ATOMS` placeholder)
keep typed `&CPMD`/`&SYSTEM`/`&DFT` text and merge geometry atoms from the C
arrays instead of rebuilding a minimal BLYP deck.
