# BP2 (best bounce solver in the universe)

Early prototype - just the core bounce algorithm without much in the way of wrapping / abstraction + an example script. 

Current dependencies:
- CasADi
- Gnuplot-iostream (https://github.com/dstahlke/gnuplot-iostream)
- Boost (iostreams, regex, system)

State of code: 
- CasadiBounceSolver.cpp needs significant refactoring
- The algorithm isn't totally finalised; need to at least look @ dynamic grid scaling

Note on performance:
- Creating a CasadiBounceSolver is "slow"
- Calling CasadiBounceSolver->solve is *fast*

IE, the optimal way to use the solver is using parameterised potentials (i.e. by temperature, other model parameters). Create a single CasadiBounceSolver, then call solve many times while varying the parameters. 