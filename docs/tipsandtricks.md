---
author:
  - LemurPwned
date: March 2023
title: Tips and tricks
---

# Tips and tricks

This is a loose collection of observations and tips that may help you in your work with the library.

- While scanning with a parameter use the previous value as the starting point. This will speed up the scan. However, remember to perturb the parameter by a small amount, otherwise the simulation may converge to a local minimum. Usually, for obtaining spectra this has no effect, but for obtaining the magnetisation profile it is important.
- Use `utils.Filters` for postprocessing the data. Not only logarithm, but detrending the spectra may help in obtaining a clearer picture. Using a `uniform_filter` from scipy may also help in smoothing the data.
- Try out integration times no lower than $10^{-12}$. For large IEC coupling values (in the ballpark of $10^{-4}$ or larger than that) you may need to go even much lower. You can always start up higher and then reduce step size to confirm that it has no effect on the results and convergence.
- Use `junction.clearLog()` and `stack.clearLogs()` to clear the log of the junction and stack. This will save you a lot of memory if you're doing a lot of scans and will vastly speed up the processing.
- You can define your own drivers! See the [API documentation](api/drivers.md) for more information.
- In models from `cmtj.models`, in the [`Solver`](api/models/sb-general-reference.md) class you can pass `prefer_numerical_roots=False` to the solver can speed up your computation, depending on the complexity of your model!
