The OpenCPMD archive patches follow
[pull request 9](https://github.com/OpenCPMD/CPMD/pull/9). `tistopgm` prints
the call stack only when `trace_depth` lies inside `trace_names`, so a stop
before `tistart` writes `LocalError` instead of faulting.
