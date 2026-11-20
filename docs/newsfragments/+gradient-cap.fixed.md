Evaluations above the legacy snapshot capacity now fail before evaluation.
They no longer return a silently truncated gradient. The ABI-compatible
`hessian[4096]` field currently limits an evaluation to 1365 atoms.
On a session that already holds another system, the call reports the topology change instead of the limit.
