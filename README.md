# fishgrowth

Fit growth models to otoliths and/or tagging data, using `RTMB` and maximum
likelihood.

The otoliths (or similar measurements of age) provide direct observed
coordinates of age and length. The tagging data provide information about the
observed length at release and length at recapture at a later time, where the
age at release is unknown and estimated as a vector of parameters.

The growth models provided by this package can be fitted to otoliths only (or
other direct measurements of age), tagging data only, or a combination of the
two. Growth variability can be modelled as constant or increasing with length.
