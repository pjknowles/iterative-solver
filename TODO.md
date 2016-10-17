- IterativeSolver is allowed to create arbitrary numbers of ParameterVector objects, but should advise for “most” of them that they do not need to be memory resident. The ParameterVector class can do what it likes with this advice.  This is not yet implemented; not sure whether it needs to be.
- Consideration of non-hermitian problems: eigensolution switches to linear equations once the eigenvalue is known to sufficient precision.
- std::string status(int level) if calling iterate directly. At the moment we have report() which uses m_verbosity.
- examples
- convergence criterion should normally be g.c but could optionally provide a metric tensor different to the kernel being solved. Convergence criterion still needs to be thought about.
- handle errors, including lack of convergency, with try/throw/catch. At the moment, lack of convergence generates a return value from solve()
- extrapolate vectors other than solution and residual. Only partly implemented.
- RSPT is unstable at high order because the perturbed wavefunctions are used directly as expansion vectors. Should rewrite in terms of orthogonal vectors.

