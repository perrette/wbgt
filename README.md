# wbgt

This code is a python wrapper around Liljegren's original C code.
Builds on [work by Max Lieblich, University of Washington, 2016](https://github.com/mdljts/wbgt), who modernized the original code so that it interfaces nicely with R (or python's ctypes).

This python package can be used to estimate wet bulb globe temperature (WBGT) from datasets of standard meterological measurements using models developed by Liljegren et al (1).
See documentation in [the original repository](https://github.com/mdljts/wbgt).

This readme focuses on modifications made since it was forked. I only keep the references at the end.

### References

1. Liljegren JC, Carhart RA, Lawday P, Tschopp S, Sharp R. Modeling the wet bulb globe temperature using standard meteorological measurements. J Occup Environ Hyg. 2008;5(10):645-55.
2. Yaglou CP, Minard D. Control of heat casualties at military training centers. AMA Arch Ind Health. 1957;16(4):302-16.
3. Lemke B, Kjellstrom T. Calculating workplace WBGT from meteorological data: a tool for climate change assessment. Ind Health. 2012;50(4):267-78.
