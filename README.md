# qpQOASES embedded

[![CircleCI](https://circleci.com/gh/mvukov/qpoases_embedded.svg?style=svg)](https://circleci.com/gh/mvukov/qpoases_embedded)

This is an unofficial fork of [qpOASES](https://github.com/coin-or/qpOASES)
QP solver. The starting point was the embedded version of the solver that
can be found in [ACADO toolkit](https://github.com/acado/acado).

Results of refactoring:
- Zero static and global variables: all solver core memory is allocated at
  object creation. This makes possible to have multiple instances running in
  different threads in the same process.
- The solver operates only on its own memory. There are zero copies of external
  data during solver operation.
- Copying of internal memory during solver execution has been dramatically
  reduced.
- Printing is pulled out of the solver. User can provide callbacks to get info
  during solver execution.
- Dead and rarely used code has been removed.
