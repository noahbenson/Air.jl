# Air.jl

An **A**tomic **I**mmutables **R**esource for Julia that's as light as air.

## Author 

Noah C. Benson <nben@uw.edu>


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

    
## License

MIT License

Copyright (c) 2019-2021 Noah C. Benson

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
