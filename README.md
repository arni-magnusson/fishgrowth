[![CRAN Status](https://r-pkg.org/badges/version/fishgrowth)](https://cran.r-project.org/package=fishgrowth)
[![CRAN Monthly](https://cranlogs.r-pkg.org/badges/fishgrowth)](https://cran.r-project.org/package=fishgrowth)
[![CRAN Total](https://cranlogs.r-pkg.org/badges/grand-total/fishgrowth)](https://cran.r-project.org/package=fishgrowth)

# fishgrowth

Fit growth models to otoliths and/or tagging data, using the `RTMB` package and
maximum likelihood.

The otoliths (or similar measurements of age) provide direct observed
coordinates of age and length. The tagging data provide information about the
observed length at release and length at recapture at a later time, where the
age at release is unknown and estimated as a vector of parameters.

The growth models provided by this package can be fitted to otoliths only (or
other direct measurements of age), tagging data only, or a combination of the
two. Growth variability can be modelled as constant or increasing with length.

## Installation

The package can be installed from
[CRAN](https://cran.r-project.org/package=fishgrowth) using the
`install.packages` command:

```R
install.packages("fishgrowth")
```

## Usage

For a summary of the package:

```R
library(fishgrowth)
?fishgrowth
```

## Development

The package is developed openly on
[GitHub](https://github.com/arni-magnusson/fishgrowth).

Feel free to open an
[issue](https://github.com/arni-magnusson/fishgrowth/issues) there if you
encounter problems or have suggestions for future versions.

The current development version can be installed using:

```R
library(remotes)
install_github("arni-magnusson/fishgrowth")
```
