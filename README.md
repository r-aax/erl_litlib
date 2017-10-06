# erl_litlib

Erlang little library for various purposes.

## inflists

Module for infinite lists support - lists containing infinite number of elements. An infinite list is represented as record containing head element, accumulator values and function for genereating next (the second) element and next accumulator.

```Erlang
-record(inflist,
{
    h :: term(),
    acc :: term(),
    f :: fun((term(), term()) -> {term(), term()})
}).
```

Main functions:

- Constructors:

```Erlang
iterate/3, iterate/2
```

- Getting infinite lists parts (head, tail, sublist):

```Erlang
hd/1, tl/1, ht/1, take/2,
is_begin/2,
nth/2,
drop/2, drop_less/2, nthtail/2,
sublist/2, sublist/3, split/2
```

- Basic simple infinite lists:

```Erlang
repeat/1, cycle/1
```

- Arithmetic series:

```Erlang
seq/2, odds/0, evens/0, seq/1, naturals/0, naturals/1
```

- Geometric series:

```Erlang
geometric_series/2, power_series/1
```

- Zip/unzip:

```Erlang
zip/2, zip_3/3, zipwith/3, unzip/1, unzip_3/1
```

- High order functions:

```Erlang
map/2, filter/2, adj_pairs_map/2, fold/3,
is_all/3, is_any/3
```

- Mathematical functions:

```Erlang
add/2, inc/1, sub/2, dec/1, neg/1,
mul/2, twice/1, dvs/2, half1/1, inv/1, ndvs/2, nhalf/1, nrem/2,
square/1, sqrt/1, cube/1, pow/2, npow/2
```

- More complex mathematical functions:

```Erlang
pows/2, npows/2,
partial_sums/1, partial_products/1, partial_avgs/1,
dirichlet_series/1, dirichlet_series/2, sign_alternate/1
```

- Some usefull infinite lists:

```Erlang
fib/0, trib/0,
harmonic_series/0, anharmonic_series/0,
grundy_series/0,
facts/0, inv_facts/0,
squares/0, sqrts/0, cubes/0, triangulars/0
```

- Prime numbers:

```Erlang
primes/0
```

- Concatenation:

```Erlang
concat/2
```

- Sparse infinite lists:

```Erlang
sparse/2, odds/1, evens/1
```

- Merge/umerge:

```Erlang
merge/2, unmerge/1
```

- Teylor's series:

```Erlang
taylor_exp/1, taylor_lnxp1/1, taylor_sin/1, taylor_cos/1, taylor_arctg/1
```

- Monotonous infinite lists actions:

```Erlang
mono_merge/2, mono_unique/1, mono_union/2, mono_intersection/2, mono_complement/2
```

# Build

```
$ rebar3 compile
```
