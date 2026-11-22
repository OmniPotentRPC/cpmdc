The `shared-dlopen-abi` test now pins `libibverbs.so.1` with
`RTLD_NOLOAD|RTLD_NODELETE` before `dlclose` so LeakSanitizer can
symbolize allocations without a suppression file.
