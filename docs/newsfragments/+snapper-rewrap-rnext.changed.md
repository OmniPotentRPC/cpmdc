Documentation prose is now rewrapped across the Org and Markdown source tree at
80 columns with updated `llms-full.txt` output.
The generator now joins only well-formed wrapped Org links during plain-text
conversion, leaves Cog markers and their multiline blocks as written, and does
not treat inline `*Fix:*` emphasis as headings.
The repository also wires `snapper-check` tooling with a pinned `snapper-fmt`
0.11.10 runtime in a minimal `prose` Pixi environment and a CI gate that runs on
documentation-only changes too.
The gate fails on files that would reformat, fused lines, and render-check
refusals, and reports the `long` and `wrap` diagnostics with counts per file.
`tools/check_snapper.sh --write` refuses symlinks and paths outside the
repository and replaces each file with a rename in its own directory.
