- Consideration of non-hermitian problems: eigensolution switches to linear equations once the eigenvalue is known to sufficient precision.
- std::string status(int level) if calling iterate directly. At the moment we have report() which uses m_verbosity.
- handle errors, including lack of convergency, with try/throw/catch. At the moment, lack of convergence generates a return value from solve()
- RSPT is unstable at high order because the perturbed wavefunctions are used directly as expansion vectors. Should rewrite in terms of orthogonal vectors.

