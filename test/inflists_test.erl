%% @doc
%% Tests for jdlib_inflists.
%%
%% @author Alexey Rybakov

% Module name.
-module(inflists_test).

% We don't need hd and tl functions for lists because we use the same names for infinite lists.
-compile({no_auto_import, [hd/1, tl/1]}).

% Unit-test using.
-include_lib("eunit/include/eunit.hrl").

% Functions import.
-import(inflists,
        [% Constructors.
         iterate/3, iterate/2,
         % Get infinite lists parts.
         hd/1, tl/1, ht/1,
         take/2, is_begin/2, nth/2, drop/2, drop_less/2, nthtail/2, sublist/2, sublist/3, split/2,
         find/2,
         % Basic simple infinite lists.
         repeat/1, cycle/1,
         % Arithmetic series.
         seq/2, odds/0, evens/0, seq/1, naturals/0, naturals/1,
         % Geometric series.
         geometric_series/2, power_series/1,
         % Zip/unzip.
         zip/2, zip_3/3, zipwith/3, unzip/1, unzip_3/1,
         % High order functions.
         map/2, filter/2, adj_pairs_map/2, fold/3, is_all/3, is_any/3,
         % Mathematical functions.
         add/2, inc/1, sub/2, dec/1, neg/1,
         mul/2, twice/1, dvs/2, half/1, inv/1, ndvs/2, nhalf/1, nrem/2,
         square/1, sqrt/1, cube/1, pow/2, npow/2,
         % More complex mathematical functions.
         partial_sums/1, partial_products/1, partial_avgs/1,
         dirichlet_series/1, dirichlet_series/2, sign_alternate/1,
         % Some usefull infinite lists.
         fib/0, trib/0, harmonic_series/0, anharmonic_series/0, grundy_series/0,
         facts/0, inv_facts/0,
         squares/0, sqrts/0, cubes/0, triangulars/0,
         % Prime numbers.
         primes/0,
         % Concatenation.
         concat/2,
         % Sparse infinite lists.         
         sparse/2, odds/1, evens/1,
         % Merge/umerge.
         merge/2, unmerge/1,
         % Teylor's series.
         taylor_exp/1, taylor_lnxp1/1, taylor_sin/1, taylor_cos/1, taylor_arctg/1,
         % Monotonic infinite lists actions.
         mono_merge/2, mono_unique/1, mono_union/2, mono_intersection/2, mono_complement/2]).

%---------------------------------------------------------------------------------------------------
% Tests.
%---------------------------------------------------------------------------------------------------

-spec take_drop_sublist_split_test() -> ok.
%% @doc
%% Functions take/drop/sublist/split test.
take_drop_sublist_split_test() ->
    IL1 = repeat(a),
    IL2 = cycle([1, 2, 3]),
    ?assertEqual([a, a, a, a, a], take(drop(IL1, 5), 5)),
    ?assertEqual([3, 1, 2, 3, 1], take(drop(IL2, 5), 5)),
    ?assert(is_begin(IL1, [a, a, a])),
    ?assert(is_begin(IL1, [a, a, a, a, a])),
    ?assert(is_begin(IL2, [1, 2, 3, 1, 2, 3, 1, 2])),
    ?assert(is_begin(IL2, [])),
    ?assert(not is_begin(IL2, [3, 2, 1])),
    ?assertEqual(a, nth(IL1, 5)),
    ?assertEqual(2, nth(IL2, 5)),
    ?assertMatch({[1, 2, 3, 1], _}, split(IL2, 4)),
    ?assertEqual([3, 1, 2], sublist(IL2, 3, 3)),
    ok.

%---------------------------------------------------------------------------------------------------

-spec repeat_test() -> ok.
%% @doc
%% Function repeat test.
repeat_test() ->
    ?assertEqual([a, a, a, a, a], take(repeat(a), 5)),
    ok.

%---------------------------------------------------------------------------------------------------

-spec cycle_test() -> ok.
%% @doc
%% Function cycle test.
cycle_test() ->
    ?assertThrow({badarg, []}, cycle([])),
    ?assertEqual([a, b, c, a, b, c, a], take(cycle([a, b, c]), 7)),
    ok.

%---------------------------------------------------------------------------------------------------

-spec arithmetic_series_test() -> ok.
%% @doc
%% Function seq test.
arithmetic_series_test() ->
    ?assertEqual([5, 7, 9], take(drop(seq(1, 2), 2), 3)),
    ?assertEqual([1, 3, 5, 7, 9], take(odds(), 5)),
    ?assertEqual([2, 4, 6, 8, 10], take(evens(), 5)),
    ?assertEqual([10, 11, 12, 13, 14], take(drop(seq(1), 9), 5)),
    ?assertEqual([1, 2, 3, 4, 5], take(naturals(), 5)),
    ok.

%---------------------------------------------------------------------------------------------------

-spec geometric_series_test() -> ok.
%% @doc
%% Function geometric_series test.
geometric_series_test() ->
    ?assertEqual([12, 24, 48], take(drop(geometric_series(3, 2), 2), 3)),
    ok.

%---------------------------------------------------------------------------------------------------

-spec zip_unzip_test() -> ok.
%% @doc
%% Zip two infinite lists.
zip_unzip_test() ->
    IL1 = seq(3, 4),
    IL2 = geometric_series(2, 2),
    IL3 = zip(IL1, IL2),
    ?assertEqual([{3, 2}, {7, 4}, {11, 8}, {15, 16}, {19, 32}], take(IL3, 5)),
    IL4 = cycle([a, b, c]),
    IL5 = zip_3(IL4, IL2, IL1),
    ?assertEqual([{c, 8, 11}, {a, 16, 15}], sublist(IL5, 3, 2)),
    IL6 = zipwith(IL1, IL2, fun(E1, E2) -> E1 + E2 end),
    ?assertEqual([5, 11, 19], take(IL6, 3)),
    {UL1, UL2} = unzip(IL3),
    ?assertEqual([3, 7, 11], take(UL1, 3)),
    ?assertEqual([2, 4, 8], take(UL2, 3)),
    {UL3, UL4, UL5} = unzip_3(IL5),
    ?assertEqual([b, c, a], sublist(UL3, 2, 3)),
    ?assertEqual([4, 8, 16], sublist(UL4, 2, 3)),
    ?assertEqual([7, 11, 15], sublist(UL5, 2, 3)),
    ok.

%---------------------------------------------------------------------------------------------------

-spec map_test() -> ok.
%% @doc
%% Function map test.
map_test() ->
    IL1 = seq(1, 2),
    IL2 = map(IL1, fun(E) -> E + 10 end),
    ?assertEqual([11, 13, 15, 17, 19], take(IL2, 5)),
    IL3 = geometric_series(1, 2),
    IL4 = map(IL3, fun(E) -> E - 1 end),
    ?assertEqual([3, 7, 15], sublist(IL4, 3, 3)),
    IL5 = seq(1),
    IL6 = adj_pairs_map(IL5, fun(X, Y) -> X + Y end),
    ?assertEqual([7, 9, 11], sublist(IL6, 3, 3)),
    IL7 = adj_pairs_map(IL5, fun(X, Y) -> Y - X end),
    ?assertEqual([1, 1, 1, 1, 1], sublist(IL7, 5, 5)),
    ?assertEqual([1, 3, 6, 10, 15], take(partial_sums(seq(1)), 5)),
    ?assertEqual([1, 2, 6, 24, 120], take(partial_products(seq(1)), 5)),
    ok.

%---------------------------------------------------------------------------------------------------

-spec fold_test() -> ok.
%% @doc
%% Function mapfold test.
fold_test() ->
    IL1 = seq(1),
    IL2 = fold(IL1, fun(X, Y) -> X + Y end, 0),
    ?assertEqual([1, 3, 6, 10, 15], take(IL2, 5)),
    ok.

%---------------------------------------------------------------------------------------------------

-spec math_test() -> ok.
%% @doc
%% Function add, sub, mul, dvs, square, sqrt, pow.
math_test() ->
    IL1 = seq(2),
    IL2 = geometric_series(2, 2),
    ?assertEqual([4, 7, 12], take(add(IL1, IL2), 3)),
    ?assertEqual([3, 4, 5], take(add(IL1, 1), 3)),
    ?assertEqual([3, 5, 9], take(add(1, IL2), 3)),
    ?assertThrow({badarg, _}, add(1, 2)),
    ?assertEqual([0, -1, -4], take(sub(IL1, IL2), 3)),
    ?assertEqual([-3, -2, -1], take(sub(IL1, 5), 3)),
    ?assertEqual([1, -1, -5], take(sub(3, IL2), 3)),
    ?assertThrow({badarg, _}, sub(1, 2)),
    ?assertEqual([4, 12, 32], take(mul(IL1, IL2), 3)),
    ?assertEqual([10, 15, 20], take(mul(IL1, 5), 3)),
    ?assertEqual([6, 12, 24], take(mul(3, IL2), 3)),
    ?assertThrow({badarg, _}, mul(1, 2)),
    ?assertEqual([2.0, 3.0, 4.0], take(sqrt(square(IL1)), 3)),
    ?assertEqual([1.0, 8.0, 27.0, 64.0, 125.0], take(pow(seq(1), 3), 5)),
    ok.

%---------------------------------------------------------------------------------------------------

-spec primes_test() -> ok.
%% @doc
%% Test for primes.
primes_test() ->
    P =
        [
              2,   3,   5,   7,  11,  13,  17,  19,  23,  29,
             31,  37,  41,  43,  47,  53,  59,  61,  67,  71,
             73,  79,  83,  89,  97, 101, 103, 107, 109, 113,
            127, 131, 137, 139, 149, 151, 157, 163, 167, 173
        ],
    ?assertEqual(P, take(primes(), 40)),
    ok.

%---------------------------------------------------------------------------------------------------

-spec concat_test() -> ok.
%% @doc
%% Function attach test.
concat_test() ->
    IL = naturals(),
    ?assertEqual([a, b, c, 1, 2, 3], take(concat([a, b, c], IL), 6)),
    ?assertEqual([a, 1, 2, 3], take(concat(a, IL), 4)),
    ok.

%---------------------------------------------------------------------------------------------------

-spec sparse_test() -> ok.
%% @doc
%% Function sparse test.
sparse_test() ->
    IL1 = seq(1),
    IL2 = sparse(IL1, 0),
    IL3 = sparse(IL1, 1),
    IL4 = sparse(IL1, 2),
    ?assertThrow({badarg, _}, sparse(IL1, -1)),
    ?assertEqual([1, 2, 3], take(IL2, 3)),
    ?assertEqual([1, 3, 5], take(IL3, 3)),
    ?assertEqual([1, 4, 7], take(IL4, 3)),
    ?assertEqual([a, a, a, a, a], take(odds(cycle([a, b])), 5)),
    ?assertEqual([b, b, b, b, b], take(evens(cycle([a, b])), 5)),
    ok.

%---------------------------------------------------------------------------------------------------

-spec merge_test() -> ok.
%% @doc
%% Function merge test.
merge_test() ->
    IL1 = seq(1),
    IL2 = cycle([a, b, c]),
    IL3 = merge(IL1, IL2),
    ?assertEqual([1, a, 2, b, 3, c, 4, a], take(IL3, 8)),
    {IL4, IL5} = unmerge(IL1),
    ?assertEqual([1, 3, 5], take(IL4, 3)),
    ?assertEqual([2, 4, 6], take(IL5, 3)),
    ok.

%---------------------------------------------------------------------------------------------------

-spec taylor_test() -> ok.
%% @doc
%% Some taylor series tests.
taylor_test() ->
    V = 0.12345,
    ?assert((lists:sum(take(taylor_exp(V), 5)) - math:exp(V)) < 0.001),
    ?assertThrow({badarg, _}, taylor_lnxp1(2)),
    ?assert((lists:sum(take(taylor_lnxp1(V), 10)) - math:log(V + 1)) < 0.001),
    ?assert((lists:sum(take(taylor_sin(V), 5)) - math:sin(V)) < 0.001),
    ?assert((lists:sum(take(taylor_cos(V), 5)) - math:cos(V)) < 0.001),
    ?assertThrow({badarg, _}, taylor_arctg(2)),
    ?assert((lists:sum(take(taylor_arctg(V), 10)) - math:atan(V)) < 0.001),
    ok.

%---------------------------------------------------------------------------------------------------

-spec mono_test() -> ok.
%% @doc
%% Tests for monotonous infinite lists.
mono_test() ->
    A = odds(),
    B = squares(),
    ?assert(is_begin(mono_merge(A, B), [1, 1, 3, 4, 5, 7, 9, 9, 11, 13, 15, 16, 17])),
    ?assert(is_begin(mono_union(A, B), [1, 3, 4, 5, 7, 9, 11, 13, 15, 16, 17])),
    ?assert(is_begin(mono_intersection(A, B), [1, 9, 25, 49, 81])),
    ?assert(is_begin(mono_complement(A, B), [3, 5, 7, 11, 13, 15, 17])),
    ok.
%---------------------------------------------------------------------------------------------------

-spec other_test() -> ok.
%% @doc
%% Other functions test.
other_test() ->
    ?assertEqual([1, -2, 3, -4, 5], take(sign_alternate(naturals()), 5)),
    ?assert(is_begin(fib(), [1, 1, 2, 3, 5, 8, 13])),
    ?assert(is_begin(trib(), [0, 1, 1, 2, 4, 7, 13, 24, 44])),
    ?assertEqual([1, -1, 1, -1], take(grundy_series(), 4)),
    ok.

%---------------------------------------------------------------------------------------------------

