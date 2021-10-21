# Persistent Lazy Dictionaries

One of the core pieces of `Air` is the persistent dictionary type,
[`PDict`](pdict.md). A nearly identical type is the `LazyDict`. Lazy
dictionaries differ from normal persistent dictionaries in that a
`LazyDict{K,V}` may be given a key and values of types `K` and `Delay{V}`,
respectively, and in this case, only calculates and yields the actual value (of
type `V`) if and when requested. Otherwise, the `LazyDict` and `LazyIdDict` types
are broadly similar to the `PDict` and `PIdDict` types.

## Examples

In progress.
