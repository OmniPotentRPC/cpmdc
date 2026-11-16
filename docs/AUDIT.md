# Documentation audit

This audit checks the cpmdc documentation at commit `cf75065`.
cpmdc is the C layer that embeds the OpenCPMD release of the Car-Parrinello Molecular Dynamics (CPMD) code.
The audit grades the nine items of the [One Good Tutorial](https://onegoodtutorial.org) checklist (version 1).
It also sorts the pages under `tutorials/`, `howto/`, `reference/`, and `explanation/`, plus the generated `api/` pages, into the four [Diataxis](https://diataxis.fr) quadrants.
Its sources are `README.md`, `LICENSE`, `meson.build`, `pixi.toml`, the Org files under `docs/orgmode/`, and the Sphinx configuration under `docs/source/`.
The `mkrst` task exports the Org files to reStructuredText pages under `docs/source/` at build time.
Git tracks the Org files but not those exported pages, so the audit reads the Org files.

## Checklist

Each row gives one checklist item a grade of pass, partial, missing, or not applicable.
Eight rows name a file behind the grade.
The evidence for those eight rows follows under the table, one paragraph per row in table order.
The Acknowledgments row is not applicable and names no file.

| Item | Status | Evidence |
| --- | --- | --- |
| Synopsis | pass | `README.md` and `docs/orgmode/index.org` each open with a synopsis. |
| Tutorial | partial | `docs/orgmode/tutorials/first-energy.org` is a lesson whose default half fails on a fresh run. |
| Contact information | missing | `README.md` and the pages under `docs/orgmode/` give no place for questions or bug reports. |
| Install instructions | pass | `docs/orgmode/howto/install.org` builds and installs the library from source. |
| Citation instructions | missing | The tree has no `CITATION.cff`, and no page explains how to cite. |
| Contribution statement | partial | `docs/orgmode/contributing/index.org` never says that contributions are welcome. |
| Reference material | pass | The pages under `docs/orgmode/reference/` and the generated `api/` pages form the reference. |
| Licensing statement | partial | `LICENSE` names the license, but `README.md` and the built site do not. |
| Acknowledgments | not applicable | |

For the synopsis, `README.md` opens with a paragraph of two sentences.
They call cpmdc the C application binary interface (ABI) layer for embedding OpenCPMD in OmniPotentRPC tools, and they list CPMD setup, geometry steps, feature discovery, and result messages.
The overview in `docs/orgmode/index.org` opens with a paragraph of three sentences.
The third sentence says that hosts pass method setup once and pass geometry for each calculation step.

`docs/orgmode/tutorials/first-energy.org` is a complete lesson, but its default half fails on a fresh run.
Its link step calls `$CC` and `$FC`.
`libcpmdc.so` depends on `libmpi.so.40` inside the Pixi shell that the tutorial opens.
The link then stops on undefined Message Passing Interface (MPI) symbols such as `MPI_Allreduce`.
Its transcripts show `0.958703657646`, `888` bytes, and `cpmdc/0.2.0`.
A run at `cf75065` prints `0.958703877646`, `816` bytes, and `cpmdc/0.3.0`.
`docs/orgmode/tutorials/eon-rgpot.org` shows no expected output for the eOn minimisation.

The tutorials name no place to ask a question or report a bug, and neither does any other page under `docs/orgmode/` or `README.md`.
The only link to the repository is a GitHub icon in the site header, set by `github_url` in `docs/source/conf.py`.

The install guide in `docs/orgmode/howto/install.org` builds the default library with `pixi run test-stub` or with three Meson commands.
The guide then installs the library and checks it with `pkg-config`.
`README.md` repeats the default build.
No page names a binary package, so every install that the docs describe builds from source.

No page in the docs names a way to cite cpmdc or OpenCPMD, and `README.md` names none either.
The tree has no `CITATION.cff`.
Two pages carry a digital object identifier (DOI), `howto/opencpmd-archive.org` and `explanation/simd-paths.org`, and both cite third-party work.

The contributing page in `docs/orgmode/contributing/index.org` lists the checks to run and the files to touch for a new keyword.
It never says that contributions are welcome or where to send them.
`README.md` has no link to it.
The start table of `docs/orgmode/index.org` links it under "Edit these docs".
The Contributing toctree on that page lists it.
`docs/orgmode/reference/feature-inventory.org` links it from the "Add a row" section.

The pages under `docs/orgmode/reference/` cover the C ABI, environment variables, `CPMDParams` fields, the CPMD option mapping, the feature inventory, the surface map, and a glossary.
The command `pixi run -e docs doxybuild` generates more pages under `api/` from the public headers.

No pages that the docs build generates name the license, and neither does `README.md`.
`LICENSE` is the Massachusetts Institute of Technology (MIT) License, with "cpmdc developers" as the copyright holder.
`docs/source/conf.py` sets the copyright line in the footer.

## Diataxis split

Each row names the pages in one quadrant and says how well they fit it.

| Quadrant | Pages | Fit |
| --- | --- | --- |
| Tutorial | `first-energy`, `eon-rgpot`, `quickstart` | `first-energy` fits the quadrant, `eon-rgpot` fits it in part, and `quickstart` does not. |
| How-to | `install`, `opencpmd-archive`, `write-cpmdparams`, `embedding`, `mpi`, `wavefunction-optimiser`, `debug-deck`, `troubleshooting` | Each page starts from a goal. |
| Reference | `c-abi`, `environment`, `cpmdparams-schema`, `cpmd-options`, `feature-inventory`, `surface-map`, `glossary`, generated `api/` | These pages fit the quadrant. |
| Explanation | `architecture`, `routes`, `wavefunction-state`, `mpi-model`, `simd-paths` | Each page explains one part of the design and fits the quadrant. |

In the Tutorial quadrant, `first-energy` is a lesson with a result at each step.
`eon-rgpot` reads like a lesson, but it needs an `eonclient` built with rgpot.
It also shows no result to compare against.
`quickstart` mixes how-to and reference material.
Its work loop, test selection, and failure tables serve readers who already know the project.

In the How-to quadrant, the `embedding` page runs long, and it repeats call contracts that `reference/c-abi` also covers.

`README.md` contains the public ABI surface, a four-row option-mapping sample, the feature identifier (ID) namespaces, and the inventory tests.
`docs/orgmode/index.org` repeats the option-mapping sample with small wording changes and has a four-row runtime-contract table.
The namespace table also appears in `docs/orgmode/reference/cpmd-options.org`, with one more row for `params.inputSections.system.*`.
The inventory tests also appear in `docs/orgmode/contributing/index.org`, and 11 of the 15 appear in `docs/orgmode/reference/feature-inventory.org`.
The ABI function list also appears, by section, in `docs/orgmode/reference/c-abi.org`.
The "Version" section in `README.md` describes 0.1.0, while `meson.build` declares 0.3.0.

## Ranked gaps

This section lists the gaps from most to least urgent.
The tutorial comes first, because One Good Tutorial treats it as the item new users need most.
A reader who follows it at `cf75065` stops at the link step of the first program.
Missing items come next, then the other partial items, then the remaining problems.
Within each group, items that a reader needs while using the code come before items needed only afterwards.

1. **The tutorial fails on a fresh run.**
   Link the host in `tutorials/first-energy.org` with the MPI compiler wrappers from the Pixi environment.
   Replace both transcripts with output from a fresh run of the default and OpenCPMD builds.
2. **Contact information is missing.**
   Name the GitHub issue tracker in the README and the docs index.
   Questions and bug reports go there.
3. **Citation instructions are missing.**
   Add a `CITATION.cff`, and add a short section on citing cpmdc and OpenCPMD to the README and the docs index.
   The project has to supply the citation record (authors, version, and any DOI), because the tree has none.
4. **The license goes unnamed.**
   Name the MIT License in the README and the docs index, and link `LICENSE`.
5. **The contribution statement is thin.**
   Start `contributing/index.org` with a sentence on which contributions are welcome and how to send them.
   Link that page from the README.
6. **Quickstart sits among the tutorials.**
   The page stays at `tutorials/quickstart.org`, so its published address, the toctree, and the links keep working.
   Its move into `howto/` is deferred.
7. **The eOn tutorial shows no output.**
   Add the `results.dat` and final energy from a real run to `tutorials/eon-rgpot.org`.
8. **Stale version text and duplicated tables.**
   Update the "Version" section in `README.md`.
   Replace the README tables that repeat the reference pages with links to those pages.

The Pixi `lychee` task passes `--config .lychee.toml`, but that file is not in the tree, so the task stops before it checks any link.

## Re-check

This section grades the checklist again on top of commit `c631ed9`.
It counts the changes that add `CITATION.cff`, the generated version
line in `README.md`, and `.lychee.toml`.

| Item | At `cf75065` | Now | Closed by |
| --- | --- | --- | --- |
| Synopsis | pass | pass | |
| Tutorial | partial | pass | the MPI link step of [#9](https://github.com/OmniPotentRPC/cpmdc/pull/9) |
| Contact information | missing | pass | the issue tracker sections of [#10](https://github.com/OmniPotentRPC/cpmdc/pull/10) |
| Install instructions | pass | pass | |
| Citation instructions | missing | pass | `CITATION.cff` and the Citation sections |
| Contribution statement | partial | pass | the opening of the contributing guide in [#10](https://github.com/OmniPotentRPC/cpmdc/pull/10) |
| Reference material | pass | pass | |
| Licensing statement | partial | pass | the License sections of [#10](https://github.com/OmniPotentRPC/cpmdc/pull/10) |
| Acknowledgments | not applicable | not applicable | |

The default half of the first tutorial, `tutorials/first-energy.org`,
ran again under GCC 13.3, Open MPI 4.1.6, and Cap'n Proto 1.0.1.
It prints the transcripts on the page: `0.958703877646`, `816` bytes,
and `cpmdc/0.3.0 available=0`.
The OpenCPMD half was not run again.

`CITATION.cff` names the author, the version, the license, and the
repository, and it carries no DOI.
`scripts/release_assert.py` regenerates its version and the README
version line from `version.txt`.
The Citation sections ask a paper with OpenCPMD energies to cite
OpenCPMD by its repository, because OpenCPMD ships no citation file.

Ranked gaps one to five are closed, and so is the README Version section
of gap eight.
Quickstart still sits among the tutorials, and its move stays deferred.
`tutorials/eon-rgpot.org` still shows no `results.dat` or final energy,
because no run of eOn through rgpot is on record.
The README tables stay, because
`tools/check_feature_inventory.py` and `tests/test_inventory_docs.py`
require the README to list each ABI symbol, the feature ID namespaces,
and the inventory tests.
A change that replaces them with links has to move those checks to the
reference pages first.

`pixi run lychee` now reads `.lychee.toml` and checks 71 links.
The two links to `CITATION.cff` on the `main` branch fail until that
file lands there.

### Gaps seven and eight

This part grades gaps seven and eight again on top of commit `c6cccb2`.
Both are closed.

Gap seven is closed. `tutorials/eon-rgpot.org` now shows the
`results.dat`, the final energy of -464.1301944953 eV, and the versions
of a real run: rgpot `5fec533f` (3.7.0), eOn `579738a93` (3.6.0),
OpenCPMD `062582b` with the archive patches, GCC 14.3.0, Open MPI
5.0.10, and Cap'n Proto 1.4.0.
The run fixed four steps that did not work as written:

- The page did not say how to build `eonclient`.
  It now has a section that builds rgpot and eOn into one prefix.
- The `pixi.toml` of eOn does not parse under Pixi 0.84.2, so the page
  builds eOn in the cpmdc Pixi environment, with Eigen 3.4 from
  `pixi global`.
- eOn builds the MPI rgpot potential only when it is configured with
  `-Dwith_mpi=enabled`, because `RGPOT_HAS_MPI` is in `rgpot-mpi.pc`
  and not in `rgpot.pc`.
- On four ranks, every rank stopped with
  `ranks_per_image=0 does not divide the MPI world`.
  rgpot opens the cpmdc session before eOn binds the calculator
  communicator, and `cpmdc_adopt_calculator_comm` refused any live
  session.
  It now refuses only a session that has evaluated, and
  the `adopt-comm` tests cover both cases.
  The four-rank run ends 1.2e-9 eV from the one-rank energy.

The first call of that minimisation reaches -17.05583351 Hartree, and the
OpenCPMD half of `tutorials/first-energy.org` ran again with the same
toolchain and matched the energy and forces on its page.
`results.dat` reports `rgpot_version` 3.2.0, which is a constant compiled
into eOn and not the version of the linked rgpot.

Gap eight is closed. The README links to `reference/c-abi.org`,
`reference/cpmd-options.org`, and `reference/feature-inventory.org`
instead of copying their tables.
`feature-inventory` now requires each ABI symbol, as a whole word, in the
C ABI reference, and each namespace term in the option reference.
`inventory-docs` requires each inventory and docs test in the feature
inventory reference, which now lists all of them.
`cpmd-public-abi-inventory` removes `cpmdc_session_create` from a copy of
the C ABI reference and expects the check to fail.

Five newsfragments linked to pages by paths that only resolve from the
repository root, so lychee failed on them.
They now give the paths as code, and lychee finds no broken link inside
the tree.

## Diataxis follow-up

`tutorials/quickstart.org` is now the how-to "Build and test cpmdc" and sits
in the How-To toctree. Its file stays under `tutorials/`, so its published
address and the links to it keep working. Its OpenCPMD section points at the
install how-to, its common failures moved to the troubleshooting how-to, and
its three test tables are one. The Meson options and Pixi tasks of the
install how-to moved to the reference page `reference/build-options.org`.
The README Documentation section links every page under `docs/orgmode`, and
`readme-navigation` fails when one is missing. `tutorials/eon-rgpot.org`
shows the `results.dat` of a real run.

The GROMACS and LAMMPS how-tos keep their steps. Their interface conventions
are in `reference/gromacs.org` and `reference/lammps.org`, their route choice
and measured runs in `explanation/md-engines.org`, and their test
descriptions in the contributing guide. The patch internals of the OpenCPMD
archive how-to, and the force-export account of the embedding how-to, are in
`explanation/opencpmd-patches.org`. The embedding how-to links to the C ABI
reference for the contracts it used to repeat. The recentring physics left
the environment reference for `explanation/isolated-cells.org`. The landing
page links to the reference tables instead of copying them.
`reference/surface-map.org` stays in the reference: it lists what the
library exports and parses.

## On develop

The audit commits now sit on develop at `c200389`, which adds the ASE
calculator `cpmdc.ase` and its how-to, `howto/ase.org`. The rename to
`CPMDCError` reaches the calculator: `CalculationFailed` is a `CPMDCError`,
and the how-to and the calculator tests use the new name. `docs-examples`
runs the blocks of the how-to.

The checks now cover the calculator. `reference-contract` skipped
`cpmdc.ase`, which had no `__all__`, and read no class members. It now
requires an `__all__` in every module of `python/cpmdc`, and each public
member of an exported class in the Python reference, which gained
`CPMDCCalculator.implemented_properties`, `ignored_changes`, and
`discard_results_on_any_change`. `readme-navigation` searched the whole
README, so it passed while the Documentation section left out
`howto/ase.org`. It now reads that section alone, and the section lists the
page. The ASE newsfragment gives its page path as code, so lychee still
finds no broken link inside the tree.

On OpenCPMD the stress of the 40 Ry periodic water deck is almost all
Pulay stress, about 16 eV over the volume, so the cell-filter example of
the ASE how-to shrinks the 10 Angstrom cube below 3 Angstrom^3 in 50
steps, and `test_bfgs_moves_the_cell_through_a_cell_filter` failed there,
on develop as well. The constant-cutoff option of
`docs/orgmode/howto/constant-cutoff.org` now removes that Pulay stress.
The how-to example and the test run with it, and the test again requires
the relaxation to end below its starting energy, within 10% of the
starting volume.

`prek run --all-files` passes, because the whitespace hooks skip the
OpenCPMD and GROMACS patches, whose blank context lines they broke. Every
Org page that the audit touches passes `snapper --check` at
`max_width = 80`.
