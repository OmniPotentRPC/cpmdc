The stopgm hook patch only adds lines. The hook still runs before the call stack is printed, and a handled stop closes the log and returns.
