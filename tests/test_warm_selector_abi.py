#!/usr/bin/env python3
"""Require matching versioned selector objects, not just module interfaces."""
from pathlib import Path
import re
import subprocess
import sys

selector = Path(sys.argv[1]).resolve().read_text()
root = Path(sys.argv[2]).resolve()
compiler = sys.argv[3:]
root.mkdir(parents=True, exist_ok=True)
versioned = "embed_choose_warm_start_metric_v2"
unversioned = "embed_choose_warm_start"
consumer = """PROGRAM selector_abi
  USE warm_selector, ONLY: real_8, embed_key_kind, embed_choose_warm_start_metric_v2
  IMPLICIT NONE
  REAL(real_8) :: geo(3), box(9), dist
  INTEGER :: source
  INTEGER(embed_key_kind) :: donor
  geo = 0
  box = [4D0,0D0,0D0,0D0,5D0,0D0,0D0,0D0,6D0]
  CALL embed_choose_warm_start_metric_v2(geo,box,.TRUE.,.5D0,dist,source,donor,metric)
  IF (source /= -2) ERROR STOP 1
CONTAINS
  FUNCTION metric(a,b,cell,periodic) RESULT(distance)
    REAL(real_8), INTENT(IN) :: a(:),b(:),cell(9)
    LOGICAL, INTENT(IN) :: periodic
    REAL(real_8) :: distance
    distance = SQRT(SUM((a-b)**2)/(SIZE(a)/3))
  END FUNCTION
END PROGRAM
"""

def run(arguments, directory, expected=0):
    result = subprocess.run(arguments, cwd=directory, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if result.returncode != expected:
        raise RuntimeError("Command failed: " + " ".join(arguments) + "\n" + result.stdout)
    return result

for label, definition, main in [
    ("matched", selector, consumer),
    ("unversioned", selector.replace(versioned, unversioned),
     consumer.replace(versioned, unversioned)),
]:
    directory = root / label
    directory.mkdir(exist_ok=True)
    (directory / "selector.F90").write_text(definition)
    (directory / "consumer.F90").write_text(main)
    run(compiler + ["-c", "selector.F90", "-o", "selector.o"], directory)
    run(compiler + ["-c", "consumer.F90", "-o", "consumer.o"], directory)
    run(compiler + ["consumer.o", "selector.o", "-o", "consumer"], directory)
    run([str(directory / "consumer")], directory)

# A valid versioned module interface cannot make an unversioned object safe.
matched = root / "matched"
legacy = root / "unversioned"
result = subprocess.run(compiler + [str(matched / "consumer.o"),
                                   str(legacy / "selector.o"), "-o", "incompatible"],
                        cwd=root, text=True, stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT)
(root / "incompatible-link.txt").write_text(result.stdout)
diagnostic = result.stdout.lower()
undefined_versioned = any(
    versioned in line and ("undefined reference" in line or "undefined symbol" in line)
    for line in diagnostic.splitlines()
)
undefined_versioned = undefined_versioned or re.search(
    r"undefined symbols for architecture[^\n]*:\n\s*[\"']?[^\s]*" + versioned,
    diagnostic,
)
if result.returncode == 0 or not undefined_versioned:
    raise RuntimeError("Incompatible selector object did not fail on the versioned symbol\n" + result.stdout)
print("Matched and unversioned controls execute; mixed selector ABI is rejected at link")
