`CITATION.cff`, `codemeta.json`, and `.zenodo.json` record how to cite cpmdc.
The version in `CITATION.cff`, `codemeta.json`, and `.zenodo.json` follows
`version.txt`.
`.zenodo.json` carries a `version` field.
Zenodo lists `version` as optional deposit metadata.
`version-sync` writes that field from `version.txt`, and `release-assert` fails
when it differs.
`CITATION.cff` references CPMD with the citation in the CPMD manual.
The CPMD reference credits the CPMD consortium as its author.
`CITATION.cff` lists the Zenodo concept DOI for all versions.
