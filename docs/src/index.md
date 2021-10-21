# Air.jl

An **A**tomic **I**mmutable **R**esource for Julia that is as light as air.

## Author 

Noah C. Benson &lt;<nben@uw.edu>&gt;


## Description 

Air is an [immutable data structure](https://en.wikipedia.org/wiki/Persistent_data_structure)
and [software transactional memory](https://en.wikipedia.org/wiki/Software_transactional_memory)
tool for Julia.  It provides several persistent data structures including
dictionaries, sets, N-dimensional arrays, heaps, weighted dictionaries, and
weighted sets, and it includes a transaction system that allows one to compose
atomic operations over such data.

Air is currently under development but includes substantial testing and is
generally stable.  Inspiration for Air's design is derived largely from
paradigms in [Clojure](https://en.wikipedia.org/wiki/Clojure).

The full API reference for Air can be found [here](API.md).

