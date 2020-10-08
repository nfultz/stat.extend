## Cran Notes 2020 Oct 8

This release is primarily to update documentation. There are some performance improvements as well, especially
for searching for discrete HDR. 

## Cran Notes 2020 Aug 5

```
Please reduce the length of the title to less than 65 characters.
```

Done.

```
If there are references describing the methods in your package, please
add these in the description field of your DESCRIPTION file in the form
authors (year) <doi:...>
authors (year) <arXiv:...>
authors (year, ISBN:...)
or if those are not available: <https:...>
with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for
auto-linking.
(If you want to add a title as well please put it in quotes: "Title")
```

There is no reference paper at this time. If and when a paper is published
in the future, we could submit a new version of the package with an updated
DESCRIPTION file.


```
Please write TRUE and FALSE instead of T and F. (Please don't use 'T' or
'F' as vector names.)
```

We believe this is a false positive in your search. 

base::T and base::F are never used.

Here are the matches:

`T` is used for Student's T distribution in a string:

```
[nfultz@nicco stat.extend]$ grep '\bT\b' R/*.R
R/HDR.R:                 paste0('Student\'s T distribution with ', df, ' degrees-of-freedom'),
R/HDR.R:                 paste0('Student\'s T distribution with ', df,
```

`F` is used  in two places, one for the F distribution, and another for the cumulative distribution function
as is standard in probability notation. F always is used as a function internally to this package.

```
[nfultz@nicco stat.extend]$ grep '\bF\b' R/*.R
R/custom.R:#' @param F a CDF of a distribution
R/custom.R:#' HDR.discrete.unimodal(.95, Q=qpois, F=ppois, lambda=1)
R/custom.R:HDR.discrete.unimodal <- function(cover.prob, Q, F, f = NULL, u = NULL, distribution = UNSPECIFIED_LABEL,
R/custom.R:  hdr(cover.prob = cover.prob, modality=discrete.unimodal, Q=Q, F=F, f=f, u=u, distribution=distribution,
R/HDR.R:                 paste0('F distribution with ', df1,
R/HDR.R:                 paste0('F distribution with ', df1,
R/HDR.R:  hdr(cover.prob, modality=discrete.unimodal, Q = qhyper, F = phyper, distribution = DIST,
R/HDR.R:  hdr(cover.prob, discrete.unimodal, Q = qgeom, F = pgeom, distribution = DIST,
R/HDR.R:  hdr(cover.prob, discrete.unimodal, Q = qbinom, F = pbinom, distribution = DIST,
R/HDR.R:  hdr(cover.prob, discrete.unimodal, Q = qpois, F = ppois, distribution = DIST,
R/HDR.R:           #Q = qnbinom, F = pbinom,
R/HDR.R:           Q = QQ, F = FF,
R/internals.R:discrete.unimodal <- function(cover.prob, Q, F, f = NULL, s = NULL, ...,
R/internals.R:  F <- partial(F, ...);
R/internals.R:  TT  <- F(MIN:MAX);
R/internals.R:  for (L in MIN:MAX) { LP     <- ifelse(L > MIN, F(L-1), 0);
R/internals.R:  P[L-MIN+1] <- F(U) - LP; }
R/internals.R:                   probability = F(U) - ifelse(L > 0, F(L-1), 0),
```

```
You also seem to be a copyright holder [cph].
Please add this information to the Authors@R field.
```
Done. Ben is the copyright holder, Neal is not.


## Test environments
* local R installation, R 3.6.0
* rhub::check_for_cran()

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
