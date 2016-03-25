# wbgt

This R package can be used to estimate wet bulb globe temperature (WBGT) from datasets of standard meterological measurements using models developed by Liljegren et al (1).  

### What is WBGT?

Wet bulb globe temperature is an established heat index for determining workplace heat stress levels in order to reduce the risk of heat-related illness in workers (2,3).  It is calculated from measured natural wet-bulb temperature, globe thermometer temperature, and dry bulb temperature (air temperature).  WBGT can also be estimated from standard meteorolgical measurements, allowing for the analysis of WBGT levels over a larger spatiotemporal scope (3).  

Several models for estimating WBGT from meteorological data have been developed, but models developed by Liljegren et al (1) are recommended for calculating outdoor WBGT (3).  Manual calculation of each WBGT value using these models can be performed using an excecutable program available from Liljegren et al (1).  

### What is in this package?

This package is currently a simple R wrapper around Liljegren et al's C code to allow for batch estimation of WBGTs from a meteorological input dataset.

### Future plans

As discussed in (3), Bernard's method for computing the WBGT appears to be more accurate indoors. We hope to include an implementation of that method in this package in the next version.

### How to install

#### Requirements

1. If you are using Windows, make sure you have installed [Rtools](https://cran.r-project.org/bin/windows/Rtools/) or installation will fail.
2. On any platform, you need a suitable C compiler installed. Here's the easiest way to make sure you have a compiler:
    1. Windows: install Rtools.
    2. Mac: install Xcode or the Command Line Tools for Xcode.
    3. Linux: different distributions have their own ways of accomplishing this (`apt`, `yum`, etc.). By default, R uses            `g++` to compile C/C++ code. You can also use `clang` if you like (see below).
3. Hadley Wickham's [excellent book](http://r-pkgs.had.co.nz/intro.html) has more information that applies to all platforms. Read (enough of) it (to make things work)!

#### Installing from GitHub, assuming the above requirements are met

Until this is on CRAN (which may never happen), the best way to install this package is using the `devtools` package. Here's an example of code you can run in the R console to install the package.
```
# run the following line if devtools is not installed
install.packages("devtools")

# this line uses devtools to install the package from GitHub
devtools::install_github("mdljts/wbgt")
```

### Acknowledgments

This product includes software produced by UChicago Argonne, LLC under Contract No. DE-AC02-06CH11357 with the Department of Energy.


### References

1. Liljegren JC, Carhart RA, Lawday P, Tschopp S, Sharp R. Modeling the wet bulb globe temperature using standard meteorological measurements. J Occup Environ Hyg. 2008;5(10):645-55. 
2. Yaglou CP, Minard D. Control of heat casualties at military training centers. AMA Arch Ind Health. 1957;16(4):302-16. 
3. Lemke B, Kjellstrom T. Calculating workplace WBGT from meteorological data: a tool for climate change assessment. Ind Health. 2012;50(4):267-78. 
