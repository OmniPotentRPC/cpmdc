The cold deck and the reference evaluator treat a deck as a cluster only when a
line of `&SYSTEM` starts with `CLUSTER`, in any case and at any indent.
`CHECK CLUSTER`, the comment text around the sections, and other sections no
longer make a periodic deck isolated or keep `CLUSTER` from going in before
`CELL VECTORS`.
