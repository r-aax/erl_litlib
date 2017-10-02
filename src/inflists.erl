%% @doc
%% Infinite lists realization.
%%
%% @author Alexey Rybakov

% Module name.
-module(inflists).

% We don't need hd and tl functions for lists because we use the same names for infinite lists.
-compile({no_auto_import, [hd/1, tl/1]}).

% Export.
-export([% Constructors.
         iterate/3, iterate/2,
         % Get infinite lists parts.
         hd/1, tl/1, ht/1,
         take/2, is_begin/2, nth/2, drop/2, drop_less/2, nthtail/2, sublist/2, sublist/3, split/2,
         % Basic simple infinite lists.
         repeat/1, cycle/1,
         % Arithmetic series.
         seq/2, odds/0, evens/0, seq/1, naturals/0, naturals/1,
         % Geometric series.
         geometric_series/2, power_series/1,
         % Zip/unzip.
         zip/2, zip_3/3, zipwith/3, unzip/1, unzip_3/1,
         % High order functions.
         map/2, filter/2, adj_pairs_map/2, mapfold/3, is_all/3, is_any/3,
         % Mathematical functions.
         add/2, inc/1, sub/2, dec/1, neg/1, mul/2, dvs/2, inv/1, ndvs/2, nrem/2,
         square/1, sqrt/1, cube/1, pow/2, npow/2,
         % More complex mathematical functions.
         pows/2, npows/2, partial_sums/1, partial_products/1, partial_avgs/1,
         dirichlet_series/1, dirichlet_series/2, sign_alternate/1,
         % Some usefull infinite lists.
         fib/0, harmonic_series/0, anharmonic_series/0, grundy_series/0, facts/0, inv_facts/0,
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
         mono_union/2, mono_intersection/2, mono_complement/2]).

%---------------------------------------------------------------------------------------------------
% Types.
%---------------------------------------------------------------------------------------------------

% Types export.
-export_type([inflist/0]).

% Define infinite list as record.
-record(inflist,
{
    h :: term(),
    acc :: term(),
    f :: fun((term(), term()) -> {term(), term()})
}).

% Inifinite list.
-type inflist() :: #inflist{}.

%---------------------------------------------------------------------------------------------------
% Infinite lists constructors.
%---------------------------------------------------------------------------------------------------

-spec iterate(H, Acc, F) -> inflist()
      when F :: fun((H, Acc) -> {H, Acc}),
           H :: term(),
           Acc :: term().
%% @doc
%% Create infinite list with head, accumulator and iterate function.
iterate(H, Acc, F) when is_function(F, 2) ->
    #inflist
    {
        h = H,
        acc = Acc,
        f = F
    };
iterate(_, _, F) ->
    throw({badarg, {F, wrong_arity}}).

%---------------------------------------------------------------------------------------------------

-spec iterate(H, F :: fun((H) -> H)) -> inflist()
      when H :: term().
%% @doc
%% Create infinite list with head and iterate function.
iterate(H, F) when is_function(F, 1) ->
    iterate
    (
        H,
        0,
        fun(Cur_H, _) ->
            {F(Cur_H), 0}
        end
    );
iterate(_, F) ->
    throw({badarg, {F, wrong_arity}}).

%---------------------------------------------------------------------------------------------------
% Get infinite lists parts (head, tail, sublists).
%---------------------------------------------------------------------------------------------------

-spec hd(IL :: inflist()) -> term().
%% @doc
%% Head of infinite list.
hd(IL) when is_record(IL, inflist) ->
    IL#inflist.h.

%---------------------------------------------------------------------------------------------------

-spec tl(IL :: inflist()) -> inflist().
%% @doc
%% Tail of infinite list.
tl(#inflist{h = H, acc = Acc, f = F} = IL) ->
    {New_H, New_Acc} = F(H, Acc),
    IL#inflist{h = New_H, acc = New_Acc}.

%---------------------------------------------------------------------------------------------------

-spec ht(IL :: inflist()) -> {term(), inflist()}.
%% @doc
%% Take head and tail simultaneously.
ht(IL) ->
    {hd(IL), tl(IL)}.

%---------------------------------------------------------------------------------------------------

-spec take(IL :: inflist(), N :: integer()) -> list().
%% @doc
%% Take first N elements of infinite list (synonym for sublist).
%% @see sublist/2
take(_, N) when (N < 0) ->
    throw({badarg, N});
take(IL, N) ->
    take(IL, N, []).

-spec take(IL :: inflist(), N :: integer(), [E]) -> [E]
      when E :: term().
%% @private
%% @doc
%% Take first elements of infinite list.
take(_, 0, R) ->
    lists:reverse(R);
take(IL, N, R) ->
    take(tl(IL), N - 1, [hd(IL) | R]).

%---------------------------------------------------------------------------------------------------

-spec is_begin(IL :: inflist(), B :: term() | list()) -> boolean().
%% @doc
%% Check if infinite list begins with given term or list.
is_begin(IL, T) when not is_list(T) ->
    hd(IL) =:= T;
is_begin(IL, []) when is_record(IL, inflist) ->
    true;
is_begin(IL, [H | T]) ->
    {IH, IT} = ht(IL),
    if
        (IH =:= H) ->
            is_begin(IT, T);
        true ->
            false
    end.

%---------------------------------------------------------------------------------------------------

-spec nth(IL :: inflist(), N :: integer()) -> term().
%% @doc
%% Take N-th element for N > 0.
nth(_, N) when (N < 1) ->
    throw({badarg, N});
nth(IL, 1) ->
    hd(IL);
nth(IL, N) ->
    nth(tl(IL), N - 1).

%---------------------------------------------------------------------------------------------------

-spec drop(IL :: inflist(), N :: integer()) -> inflist().
%% @doc
%% Drop first N elements of infinite list (synonym for nthtail).
%% @see nthtail/2
drop(_, N) when (N < 0) ->
    throw({badarg, N});
drop(IL, 0) ->
    IL;
drop(IL, N) ->
    drop(tl(IL), N - 1).

%---------------------------------------------------------------------------------------------------

-spec drop_less(IL :: inflist(), N :: number()) -> inflist().
%% @doc
%% Drop first elements less than N.
drop_less(IL, N) ->
    {H, T} = ht(IL),
    if
        H < N ->
            drop_less(T, N);
        true ->
            IL
    end.

%---------------------------------------------------------------------------------------------------

-spec nthtail(IL :: inflist(), N :: integer()) -> inflist().
%% @doc
%% Tail of infinite list without N elements (synonym for drop).
%% @see drop/2
nthtail(IL, N) ->
    drop(IL, N).

%---------------------------------------------------------------------------------------------------

-spec sublist(IL :: inflist(), N :: integer()) -> list().
%% @doc
%% Sublist from first position (synonym for take).
%% @see take/2
sublist(IL, N) ->
    take(IL, N).

%---------------------------------------------------------------------------------------------------

-spec sublist(IL :: inflist(), Start :: integer(), N :: integer()) -> list().
%% @doc
%% Sublist from given position.
sublist(IL, Start, N) ->
    take(drop(IL, Start - 1), N).

%---------------------------------------------------------------------------------------------------

-spec split(IL :: inflist(), N :: integer()) -> {list(), inflist()}.
%% @doc
%% Split infinite list by position.
split(_, N) when (N < 0) ->
    throw({badarg, N});
split(IL, N) ->
    split(IL, N, []).

-spec split(IL :: inflist(), N :: integer(), [E]) -> {[E], inflist()}
      when E :: term().
%% @doc
%% Split infinite list by position.
split(IL, 0, R) ->
    {lists:reverse(R), IL};
split(IL, N, R) ->
    split(tl(IL), N - 1, [hd(IL) | R]).

%---------------------------------------------------------------------------------------------------
% Basic simple infinite lists.
%---------------------------------------------------------------------------------------------------

-spec repeat(T :: term()) -> inflist().
%% @doc
%% Construct infinite list, containing one repeating element (Haskell analogue).
%% T -> [T, T, ..]
repeat(T) ->
    iterate
    (
        T,
        fun(_) ->
            T
        end
    ).

%---------------------------------------------------------------------------------------------------

-spec cycle(L :: list()) -> inflist().
%% @doc
%% Construct infinite list, containing infinite number of list L copies (Haskell analogue).
%% [A, B, C] -> [A, B, C, A, B, C, ..]
cycle([]) ->
    throw({badarg, []});
cycle([H | T]) ->
    iterate
    (
        H,
        T,
        fun
            (_, []) ->
                {H, T};
            (_, [Cur_H | Cur_T]) ->
                {Cur_H, Cur_T}
        end
    ).

%---------------------------------------------------------------------------------------------------
% Arithmetic series.
%---------------------------------------------------------------------------------------------------

-spec seq(From :: number(), Step :: number()) -> inflist().
%% @doc
%% Construct inflinite list [From, From + Step, From + 2 * Step, ..].
seq(From, Step) ->
    iterate
    (
        From,
        fun(H) ->
            H + Step
        end
    ).

%---------------------------------------------------------------------------------------------------

-spec odds() -> inflist().
%% @doc
%% Odd natural numbers.
odds() ->
    seq(1, 2).

%---------------------------------------------------------------------------------------------------

-spec evens() -> inflist().
%% @doc
%% Even natural numbers.
evens() ->
    seq(2, 2).

%---------------------------------------------------------------------------------------------------

-spec seq(From :: number()) -> inflist().
%% @doc
%% Construct infinite list [From, From + 1, From + 2, ..] (synonym for naturals).
%% @see naturals/1
seq(From) ->
    seq(From, 1).

%---------------------------------------------------------------------------------------------------

-spec naturals() -> inflist().
%% @doc
%% Infinite list of natural numbers.
naturals() ->
    seq(1).

%---------------------------------------------------------------------------------------------------

-spec naturals(From :: integer()) -> inflist().
%% @doc
%% Naturals from given number (synonym for seq).
%% @see seq/1
naturals(From) ->
    seq(From).

%---------------------------------------------------------------------------------------------------
% Geometric series.
%---------------------------------------------------------------------------------------------------

-spec geometric_series(Base :: number(), K :: number()) -> inflist().
%% @doc
%% Construct infinite list [Base, Base * K, Base * K^2, ..].
geometric_series(Base, K) ->
    iterate
    (
        Base,
        fun(H) ->
            H * K
        end
    ).

%---------------------------------------------------------------------------------------------------

-spec power_series(X :: number()) -> inflist().
%% @doc
%% Series of number powers.
power_series(X) ->
    geometric_series(1, X).

%---------------------------------------------------------------------------------------------------
% Zip/unzip functions and functors.
%---------------------------------------------------------------------------------------------------

-spec zip(IL1 :: inflist(), IL2 :: inflist()) -> inflist().
%% @doc
%% Zip two infinite lists.
zip(#inflist{h = H1, acc = Acc1, f = F1}, #inflist{h = H2, acc = Acc2, f = F2}) ->
    iterate
    (
        {H1, H2},
        {Acc1, Acc2},
        fun({Cur_H1, Cur_H2}, {Cur_Acc1, Cur_Acc2}) ->
            {New_H1, New_Acc1} = F1(Cur_H1, Cur_Acc1),
            {New_H2, New_Acc2} = F2(Cur_H2, Cur_Acc2),
            {{New_H1, New_H2}, {New_Acc1, New_Acc2}}
        end
    );
zip(IL1, IL2) ->
    throw({badarg, {IL1, IL2}}).

%---------------------------------------------------------------------------------------------------

-spec zip_3(IL1 :: inflist(), IL2 :: inflist(), IL3 :: inflist()) -> inflist().
%% @doc
%% Zip three infinite lists.
zip_3(#inflist{h = H1, acc = Acc1, f = F1},
      #inflist{h = H2, acc = Acc2, f = F2},
      #inflist{h = H3, acc = Acc3, f = F3}) ->
    iterate
    (
        {H1, H2, H3},
        {Acc1, Acc2, Acc3},
        fun({Cur_H1, Cur_H2, Cur_H3}, {Cur_Acc1, Cur_Acc2, Cur_Acc3}) ->
            {New_H1, New_Acc1} = F1(Cur_H1, Cur_Acc1),
            {New_H2, New_Acc2} = F2(Cur_H2, Cur_Acc2),
            {New_H3, New_Acc3} = F3(Cur_H3, Cur_Acc3),
            {{New_H1, New_H2, New_H3}, {New_Acc1, New_Acc2, New_Acc3}}
        end
    );
zip_3(IL1, IL2, IL3) ->
    throw({badarg, {IL1, IL2, IL3}}).

%---------------------------------------------------------------------------------------------------

-spec zipwith(IL1 :: inflist(), IL2 :: inflist(), Zip_F) -> inflist()
      when Zip_F :: fun((T1, T2) -> {T1, T2}),
      T1 :: term(),
      T2 :: term().
%% @doc
%% Zip two infinite lists with given function.
zipwith(#inflist{h = H1, acc = Acc1, f = F1},
        #inflist{h = H2, acc = Acc2, f = F2},
        Zip_F) when is_function(Zip_F, 2) ->
    iterate
    (
        Zip_F(H1, H2),
        {{H1, Acc1}, {H2, Acc2}},
        fun(_, {{Cur_H1, Cur_Acc1}, {Cur_H2, Cur_Acc2}}) ->
            {New_H1, New_Acc1} = F1(Cur_H1, Cur_Acc1),
            {New_H2, New_Acc2} = F2(Cur_H2, Cur_Acc2),
            {Zip_F(New_H1, New_H2), {{New_H1, New_Acc1}, {New_H2, New_Acc2}}}
        end
    );
zipwith(IL1, IL2, Zip_F) ->
    throw({badarg, {IL1, IL2, Zip_F}}).

%---------------------------------------------------------------------------------------------------

-spec unzip(IL :: inflist()) -> {inflist(), inflist()}.
%% @doc
%% Unzip infinite list into two lists.
unzip(#inflist{h = {H1, H2}, acc = {Acc1, Acc2}, f = F}) ->
    {
        iterate
        (
            H1,
            Acc1,
            fun(Cur_H1, Cur_Acc1) ->
                {{New_H1, _}, {New_Acc1, _}} = F({Cur_H1, H2}, {Cur_Acc1, Acc2}),
                {New_H1, New_Acc1}
            end
        ),
        iterate
        (
            H2,
            Acc2,
            fun(Cur_H2, Cur_Acc2) ->
                {{_, New_H2}, {_, New_Acc2}} = F({H1, Cur_H2}, {Acc1, Cur_Acc2}),
                {New_H2, New_Acc2}
            end
        )
    };
unzip(IL) ->
    throw({badarg, IL}).

%---------------------------------------------------------------------------------------------------

-spec unzip_3(IL :: inflist()) -> {inflist(), inflist(), inflist()}.
%% @doc
%% Unzip infinite list into three lists.
unzip_3(#inflist{h = {H1, H2, H3}, acc = {Acc1, Acc2, Acc3}, f = F}) ->
    {
        iterate
        (
            H1,
            Acc1,
            fun(Cur_H1, Cur_Acc1) ->
                {{New_H1, _, _}, {New_Acc1, _, _}} = F({Cur_H1, H2, H3}, {Cur_Acc1, Acc2, Acc3}),
                {New_H1, New_Acc1}
            end
        ),
        iterate
        (
            H2,
            Acc2,
            fun(Cur_H2, Cur_Acc2) ->
                {{_, New_H2, _}, {_, New_Acc2, _}} = F({H1, Cur_H2, H3}, {Acc1, Cur_Acc2, Acc3}),
                {New_H2, New_Acc2}
            end
        ),
        iterate
        (
            H3,
            Acc3,
            fun(Cur_H3, Cur_Acc3) ->
                {{_, _, New_H3}, {_, _, New_Acc3}} = F({H1, H2, Cur_H3}, {Acc1, Acc2, Cur_Acc3}),
                {New_H3, New_Acc3}
            end
        )
    };
unzip_3(IL) ->
    throw({badarg, IL}).

%---------------------------------------------------------------------------------------------------
% High order functions.
%---------------------------------------------------------------------------------------------------

-spec map(IL :: inflist(), Map_F :: fun((term()) -> term())) -> inflist().
%% @doc
%% Apply function to every element of infinite list.
map(#inflist{h = H, acc = Acc, f = F}, Map_F) when is_function(Map_F, 1) ->
    iterate
    (
        Map_F(H),
        {H, Acc},
        fun(_, {Cur_H, Cur_Acc}) ->
            {New_H, New_Acc} = F(Cur_H, Cur_Acc),
            {Map_F(New_H), {New_H, New_Acc}}
        end
    );
map(IL, Map_F) ->
    throw({badarg, {IL, Map_F}}).

%---------------------------------------------------------------------------------------------------

-spec filter(IL :: inflist(), Filter_F :: fun((term()) -> boolean())) -> inflist().
%% @doc
%% Filter infinite list.
%% Warning! Dangerous function.
%% It can cause infinite recursion if result list is not finite. 
filter(IL, Filter_F) ->
    New_IL =
        iterate
        (
            0,
            IL,
            fun
                F_(_, L) ->
                    {H, T} = ht(L),
                    F_Res = Filter_F(H),
                    if
                        F_Res ->
                            {H, T};
                        true ->
                            F_(none, T)
                    end
            end
        ),
    tl(New_IL).

%---------------------------------------------------------------------------------------------------

-spec adj_pairs_map(IL :: inflist(), Map_F :: fun((term(), term()) -> term())) -> inflist().
%% @doc
%% Apply map function to every pair of adjacent elements.
adj_pairs_map(IL, Map_F) ->
    zipwith(IL, tl(IL), Map_F).

%---------------------------------------------------------------------------------------------------

-spec mapfold(IL :: inflist(), Fold_F :: fun((term(), Fold_Acc) -> Fold_Acc), Fold_Acc) -> inflist()
      when Fold_Acc :: term().
%% @doc
%% Fold infinite list with saving accumulators.
%% Note. We will never reach fold result because list is infinite.
mapfold(#inflist{h = H, acc = Acc, f = F}, Fold_F, Fold_Acc) when is_function(Fold_F, 2) ->
    iterate
    (
        Fold_F(H, Fold_Acc),
        {H, Acc},
        fun(Cur_Fold_Acc, {Cur_H, Cur_Acc}) ->
            {New_H, New_Acc} = F(Cur_H, Cur_Acc),
            {Fold_F(New_H, Cur_Fold_Acc), {New_H, New_Acc}}
        end
    );
mapfold(IL, Fold_F, Fold_Acc) ->
    throw({badarg, {IL, Fold_F, Fold_Acc}}).

%---------------------------------------------------------------------------------------------------

-spec is_all(IL :: inflist(), Pred :: fun((term()) -> boolean()), N :: integer()) -> boolean().
%% @doc
%% Check predicate for all of N first infinite list members.
is_all(_, _, 0) ->
    true;
is_all(IL, Pred, N) when (N > 0) ->
    Is = Pred(hd(IL)),
    if
        not Is ->
            false;
        true ->
            is_all(tl(IL), Pred, N - 1)
    end.

%---------------------------------------------------------------------------------------------------

-spec is_any(IL :: inflist(), Pred :: fun((term()) -> boolean()), N :: integer()) -> boolean().
%% @doc
%% Check predicate for any of N first infinite list members.
is_any(_, _, 0) ->
         false;
is_any(IL, Pred, N) when (N > 0) ->
    Is = Pred(hd(IL)),
    if
        Is ->
            true;
        true ->
            is_any(tl(IL), Pred, N - 1)
    end.

%---------------------------------------------------------------------------------------------------
% Mathematical functions.
%---------------------------------------------------------------------------------------------------

-spec add(Arg, Arg) -> inflist()
      when Arg :: inflist() | term().
%% @doc
%% Add function.
add(A, B) ->
    Is_A = is_record(A, inflist),
    Is_B = is_record(B, inflist),
    if
        Is_A andalso Is_B ->
            zipwith(A, B, fun(X, Y) -> X + Y end);
        Is_A ->
            map(A, fun(X) -> X + B end);
        Is_B ->
            map(B, fun(X) -> X + A end);
        true ->
            throw({badarg, {A, B}})
    end.

%---------------------------------------------------------------------------------------------------

-spec inc(IL :: inflist()) -> inflist().
%% @doc
%% Increment.
inc(IL) ->
    add(IL, 1).

%---------------------------------------------------------------------------------------------------

-spec sub(Arg, Arg) -> inflist()
      when Arg :: inflist() | term().
%% @doc
%% Sub function.
sub(A, B) ->
    Is_A = is_record(A, inflist),
    Is_B = is_record(B, inflist),
    if
        Is_A andalso Is_B ->
            zipwith(A, B, fun(X, Y) -> X - Y end);
        Is_A ->
            map(A, fun(X) -> X - B end);
        Is_B ->
            map(B, fun(X) -> A - X end);
        true ->
            throw({badarg, {A, B}})
    end.

%---------------------------------------------------------------------------------------------------

-spec dec(IL :: inflist()) -> inflist().
%% @doc
%% Decrement.
dec(IL) ->
    sub(IL, 1).

%---------------------------------------------------------------------------------------------------

-spec neg(IL :: inflist()) -> inflist().
%% @doc
%% Negate infinite list.
neg(IL) ->
    map(IL, fun(X) -> -X end).

%---------------------------------------------------------------------------------------------------

-spec mul(Arg, Arg) -> inflist()
      when Arg :: inflist() | term().
%% @doc
%% Multiplication.
mul(A, B) ->
    Is_A = is_record(A, inflist),
    Is_B = is_record(B, inflist),
    if
        Is_A andalso Is_B ->
            zipwith(A, B, fun(X, Y) -> X * Y end);
        Is_A ->
            map(A, fun(X) -> X * B end);
        Is_B ->
            map(B, fun(X) -> X * A end);
        true ->
            throw({badarg, {A, B}})
    end.

%---------------------------------------------------------------------------------------------------

-spec dvs(Arg, Arg) -> inflist()
      when Arg :: inflist() | term().
%% @doc
%% Division.
dvs(A, B) ->
    Is_A = is_record(A, inflist),
    Is_B = is_record(B, inflist),
    if
        Is_A andalso Is_B ->
            zipwith(A, B, fun(X, Y) -> X / Y end);
        Is_A ->
            map(A, fun(X) -> X / B end);
        Is_B ->
            map(B, fun(X) -> A / X end);
        true ->
            throw({badarg, {A, B}})
    end.

%---------------------------------------------------------------------------------------------------

-spec inv(IL :: inflist()) -> inflist().
%% @doc
%% Invert list.
inv(IL) ->
    map(IL, fun(X) -> 1 / X end).

%---------------------------------------------------------------------------------------------------

-spec ndvs(Arg, Arg) -> inflist()
      when Arg :: inflist() | term().
%% @doc
%% Division without remainder (for natural numbers).
ndvs(A, B) ->
    Is_A = is_record(A, inflist),
    Is_B = is_record(B, inflist),
    if
        Is_A andalso Is_B ->
            zipwith(A, B, fun(X, Y) -> X div Y end);
        Is_A ->
            map(A, fun(X) -> X div B end);
        Is_B ->
            map(B, fun(X) -> A div X end);
        true ->
            throw({badarg, {A, B}})
    end.

%---------------------------------------------------------------------------------------------------

-spec nrem(Arg, Arg) -> inflist()
      when Arg :: inflist() | term().
%% @doc
%% Get remainder (infinite list of remainders).
nrem(A, B) ->
    Is_A = is_record(A, inflist),
    Is_B = is_record(B, inflist),
    if
        Is_A andalso Is_B ->
            zipwith(A, B, fun(X, Y) -> X rem Y end);
        Is_A ->
            map(A, fun(X) -> X rem B end);
        Is_B ->
            map(B, fun(X) -> A rem X end);
        true ->
            throw({badarg, {A, B}})
    end.

%---------------------------------------------------------------------------------------------------

-spec square(IL :: inflist()) -> inflist().
%% @doc
%% Infinite list of squares.
square(IL) ->
    mul(IL, IL).

%---------------------------------------------------------------------------------------------------

-spec sqrt(IL :: inflist()) -> inflist().
%% @doc
%% Square root of infinite list.
sqrt(IL) ->
    map(IL, fun(X) -> math:sqrt(X) end).

%---------------------------------------------------------------------------------------------------

-spec cube(IL :: inflist()) -> inflist().
%% @doc
%% Cobe of infinite list.
cube(IL) ->
    mul(square(IL), IL).

%---------------------------------------------------------------------------------------------------

-spec pow(IL :: inflist(), P :: number()) -> inflist().
%% @doc
%% Power of infinite list.
pow(IL, P) ->
    map(IL, fun(X) -> math:pow(X, P) end).

%---------------------------------------------------------------------------------------------------

-spec npow(IL :: inflist(), P :: number()) -> inflist().
%% @doc
%% Power of infinite list.
npow(IL, 1) ->
    IL;
npow(IL, P) when (P > 1) ->
    mul(IL, npow(IL, P - 1)).

%---------------------------------------------------------------------------------------------------
% More complex mathematical functions.
%---------------------------------------------------------------------------------------------------

-spec pows(N :: number(), IL :: inflist()) -> inflist().
%% @doc
%% Return powers of given number (powers are taken from IL).
pows(N, IL) ->
    map(IL, fun(X) -> math:pow(N, X) end).

%---------------------------------------------------------------------------------------------------

-spec npows(N :: number(), IL :: inflist()) -> inflist().
%% @doc
%% Return natural powers of given number (powers are taken from IL).
npows(N, IL) ->
    map
    (
        IL,
        fun(P) ->
            N_Pow =
                fun
                    NP_(N_, 1) ->
                        N_;
                    NP_(N_, P_) when (P_ > 1) ->
                        N_ * NP_(N_, P_ - 1)
                end,
            N_Pow(N, P)
        end
    ).

%---------------------------------------------------------------------------------------------------
% Partial sums and products of infinite list.
%---------------------------------------------------------------------------------------------------

-spec partial_sums(IL :: inflist()) -> inflist().
%% @doc
%% Partial sums.
partial_sums(IL) ->
    mapfold(IL, fun(X, Y) -> X + Y end, 0).

%---------------------------------------------------------------------------------------------------

-spec partial_products(IL :: inflist()) -> inflist().
%% @doc
%% Partial products.
partial_products(IL) ->
    mapfold(IL, fun(X, Y) -> X * Y end, 1).

%---------------------------------------------------------------------------------------------------

-spec partial_avgs(IL :: inflist()) -> inflist().
%% @doc
%% Average values of infinite list.
partial_avgs(IL) ->
    dvs(partial_sums(IL), naturals()).

%---------------------------------------------------------------------------------------------------

-spec dirichlet_series(S :: number()) -> inflist().
%% @doc
%% Dirichlet series sum 1 / n^s.
dirichlet_series(S) ->
    pow(harmonic_series(), S).

%---------------------------------------------------------------------------------------------------

-spec dirichlet_series(IL :: inflist(), S :: number()) -> inflist().
%% @doc
%% Dirichlet series of base inflist IL and pow degree S.
dirichlet_series(IL, S) ->
    mul(IL, dirichlet_series(S)).

%---------------------------------------------------------------------------------------------------

-spec sign_alternate(IL :: inflist()) -> inflist().
%% @doc
%% Alternate sign of infinite list.
%% Odd position elements are unchanged, even position elements are negated.
sign_alternate(IL) ->
    mul(IL, grundy_series()).

%---------------------------------------------------------------------------------------------------
% Some usefull sequences.
%---------------------------------------------------------------------------------------------------

-spec fib() -> inflist().
%% @doc
%% Fibonacci numbers.
fib() ->
    iterate
    (
        1,
        0,
        fun(H, Acc) ->
            {H + Acc, H}
        end
    ).

%---------------------------------------------------------------------------------------------------

-spec harmonic_series() -> inflist().
%% @doc
%% Harmonic series.
harmonic_series() ->
    inv(naturals()).

%---------------------------------------------------------------------------------------------------

-spec anharmonic_series() -> inflist().
%% @doc
%% Anharmonic series (Leibniz series).
anharmonic_series() ->
    sign_alternate(inv(odds())).

%---------------------------------------------------------------------------------------------------

-spec grundy_series() -> inflist().
%% @doc
%% Grundy series (see Patrick Carmelo Grundy).
grundy_series() ->
    cycle([1, -1]).

%---------------------------------------------------------------------------------------------------

-spec facts() -> inflist().
%% @doc
%% Factorials series (starts with 0! = 1).
facts() ->
    iterate
    (
        1,
        1,
        fun(H, Acc) ->
            M = H * Acc,
            {M, Acc + 1}
        end
    ).

%---------------------------------------------------------------------------------------------------

-spec inv_facts() -> inflist().
%% @doc
%% Series of inverted factorials.
inv_facts() ->
    inv(facts()).

%---------------------------------------------------------------------------------------------------

-spec squares() -> inflist().
%% @doc
%% Natural squares.
squares() ->
    square(naturals()).

%---------------------------------------------------------------------------------------------------

-spec sqrts() -> inflist().
%% @doc
%% Square roots of naturals.
sqrts() ->
    sqrt(naturals()).

%---------------------------------------------------------------------------------------------------

-spec cubes() -> inflist().
%% @doc
%% Natural cubes.
cubes() ->
    cube(naturals()).

%---------------------------------------------------------------------------------------------------

-spec triangulars() -> inflist().
%% @doc
%% Triangulars numbers.
triangulars() ->
    iterate
    (
        1,
        1,
        fun(H, Acc) ->
            {H + Acc + 1, Acc + 1}
        end
    ).

%---------------------------------------------------------------------------------------------------
% Prime numbers.
%---------------------------------------------------------------------------------------------------

-spec primes() -> inflist().
%% @doc
%% Prime numbers.
primes() ->
    iterate
    (
        2,   % first prime number
        [2], % list of factors

        fun(H, Acc) ->

            % Get next prime number function.
            Next_P_Fun =
                fun
                    NPF_(N, [], _) ->
                        N;
                    NPF_(N, [P | _], _) when (P * P > N) ->
                        N;
                    NPF_(N, [P | T], Ps) ->
                        if
                            (N rem P) =:= 0 ->
                                NPF_(N + 1, Ps, Ps);
                            true ->
                                NPF_(N, T, Ps)
                        end
                end,

            P = Next_P_Fun(H + 1, Acc, Acc),
            {P, Acc ++ [P]}
        end
    ).

%---------------------------------------------------------------------------------------------------
% Concatenate.
%---------------------------------------------------------------------------------------------------

-spec concat(L :: list(), IL :: inflist()) -> inflist().
%% @doc
%% Attach list to infinite list from the beginning.
concat([], IL) when is_record(IL, inflist) ->
    IL;
concat(T, IL) when not is_list(T) ->
    concat([T], IL);
concat([LH | LT], #inflist{h = H, acc = Acc, f = F}) ->
    iterate
    (
        LH,
        {false, LT},
        fun
            (_, {false, Cur_L}) ->
                case Cur_L of
                    [Cur_LH | Cur_LT] ->
                        {Cur_LH, {false, Cur_LT}};
                    [] ->
                        {H, {true, Acc}}
                end;
            (Cur_H, {true, Cur_Acc}) ->
                {New_H, New_Acc} = F(Cur_H, Cur_Acc),
                {New_H, {true, New_Acc}}
        end
    );
concat(L, IL) ->
    throw({badarg, {L, IL}}).

%---------------------------------------------------------------------------------------------------
% Sparse.
%---------------------------------------------------------------------------------------------------

-spec sparse(IL :: inflist(), N :: integer()) -> inflist().
%% @doc
%% Take sparse infinite list (first element, and then every (N + 1)-th).
%% For example:
%% sparse(IL, 0) = IL
%% sparse([E1, E2, E3, ...], 1) = [E1, E3, E5, ...]
%% sparse([E1, E2, E3, ...], 2) = [E1, E4, E7, ...]
sparse(IL, 0) when is_record(IL, inflist) ->
    IL;
sparse(#inflist{h = H, acc = Acc, f = F}, N) when (is_integer(N) andalso (N > 0)) ->
    iterate
    (
        H,
        Acc,
        fun(Cur_H, Cur_Acc) ->
            FN =
                fun
                    Loc_FN(Loc_H, Loc_Acc, 0) ->
                        {Loc_H, Loc_Acc};
                    Loc_FN(Loc_H, Loc_Acc, Loc_N) ->
                        {New_H, New_Acc} = F(Loc_H, Loc_Acc),
                        Loc_FN(New_H, New_Acc, Loc_N - 1)
                end,
            FN(Cur_H, Cur_Acc, N + 1)
        end
    );
sparse(IL, N) ->
    throw({badarg, {IL, N}}).

%---------------------------------------------------------------------------------------------------

-spec odds(IL :: inflist()) -> inflist().
%% @doc
%% Odd elements of list.
odds(IL) ->
    sparse(IL, 1).

%---------------------------------------------------------------------------------------------------

-spec evens(IL :: inflist()) -> inflist().
%% @doc
%% Even elements of list.
evens(IL) ->
    sparse(tl(IL), 1).

%---------------------------------------------------------------------------------------------------
% Merge/unmerge.
%---------------------------------------------------------------------------------------------------

-spec merge(IL1 :: inflist(), IL2 :: inflist()) -> inflist().
%% @doc
%% Merge two infinite lists.
merge(#inflist{h = H1, acc = Acc1, f = F1}, #inflist{h = H2, acc = Acc2, f = F2}) ->
    iterate
    (
        H1,
        {{H2, Acc2}, F1(H1, Acc1), false},
        fun(_, {{Cur_H, Cur_Acc}, Next, Is_F1}) ->
            {
                Cur_H,
                {
                    Next,
                    (if Is_F1 -> F1; true -> F2 end)(Cur_H, Cur_Acc),
                    not Is_F1
                }
            }
        end
    );
merge(IL1, IL2) ->
    throw({badarg, {IL1, IL2}}).

%---------------------------------------------------------------------------------------------------

-spec unmerge(IL :: inflist()) -> {inflist(), inflist()}.
%% @doc
%% Split infinite list to odd and even elements infinite lists.
unmerge(IL) ->
    {odds(IL), evens(IL)}.

%---------------------------------------------------------------------------------------------------
% Taylor series.
%---------------------------------------------------------------------------------------------------

-spec taylor_exp(X :: number()) -> inflist().
%% @doc
%% Taylor series of e^x for (-inf, inf).
%% <pre>
%%           x    x^2   x^3
%% e^x = 1 + -  + --- + --- + ...
%%           1!    2!    3!
%% </pre>
taylor_exp(X) ->
    dvs(power_series(X), facts()).

%---------------------------------------------------------------------------------------------------

-spec taylor_lnxp1(X :: number()) -> inflist().
%% @doc
%% Taylor series of ln(1 + x) for (-1, 1].
%% <pre>
%%                 x^2   x^3
%% ln(x + 1) = x - --- + --- - ...
%%                  2     3
%% </pre>
taylor_lnxp1(X) when ((X =< -1) orelse (X > 1)) ->
    throw({badarg, X});
taylor_lnxp1(X) ->
    sign_alternate(dvs(tl(power_series(X)), naturals())).

%---------------------------------------------------------------------------------------------------

-spec taylor_sin(X :: number()) -> inflist().
%% @doc
%% Taylor series of sin(x) for (-inf, inf).
%% <pre>
%%              x^3   x^5
%% sin(x) = x - --- + --- - ...
%%               3!    5!
%% </pre>
taylor_sin(X) ->
    sign_alternate(evens(taylor_exp(X))).

%---------------------------------------------------------------------------------------------------

-spec taylor_cos(X :: number()) -> inflist().
%% @doc
%% Taylor series of cos(x) for (-inf, inf).
%% <pre>
%%              x^2   x^4
%% cos(x) = 1 - --- + --- - ...
%%               2!    4!
%% </pre>
taylor_cos(X) ->
    sign_alternate(odds(taylor_exp(X))).

%---------------------------------------------------------------------------------------------------

-spec taylor_arctg(X :: number()) -> inflist().
%% @doc
%% Taylor series of arctg(x) for (-1, 1).
%% <pre>
%%                x^3   x^5
%% arctg(x) = x - --- + --- - ...
%%                 3     5
%% </pre>
taylor_arctg(X) when (abs(X) >= 1) ->
    throw({badarg, X});
taylor_arctg(X) ->
    sign_alternate(dvs(evens(power_series(X)), naturals())).

%---------------------------------------------------------------------------------------------------
% Sets operations.
%---------------------------------------------------------------------------------------------------

-spec mono_union(IL1 :: inflist(), IL2 :: inflist()) -> inflist().
%% @doc
%% Union of two monotonous lists.
mono_union(IL1, IL2) ->
    IL =
        iterate
        (
            0,
            {IL1, IL2},
            fun(_, {L1, L2}) ->
                {H1, T1} = ht(L1),
                {H2, T2} = ht(L2),
                if
                    H1 < H2 ->
                        {H1, {T1, L2}};
                    true ->
                        {H2, {L1, T2}}
                end
            end
        ),
    tl(IL).

%---------------------------------------------------------------------------------------------------

-spec mono_intersection(IL1 :: inflist(), IL2 :: inflist()) -> inflist().
%% @doc
%% Union of two monotonous lists.
mono_intersection(IL1, IL2) ->
    IL =
        iterate
        (
            0,
            {IL1, IL2},
            fun
                F_(_, {L1, L2}) ->
                    {H1, T1} = ht(L1),
                    {H2, T2} = ht(L2),
                    if
                        H1 < H2 ->
                            F_(none, {drop_less(L1, H2), L2});
                        H1 > H2 ->
                            F_(none, {L1, drop_less(L2, H1)});
                        true ->
                            {H1, {T1, T2}}
                    end
            end
        ),
    tl(IL).

%---------------------------------------------------------------------------------------------------

-spec mono_complement(IL1 :: inflist(), IL2 :: inflist()) -> inflist().
%% @doc
%% Union of two monotonous lists.
mono_complement(IL1, IL2) ->
    IL =
        iterate
        (
            0,
            {IL1, IL2},
            fun
                F_(_, {L1, L2}) ->
                    {H1, T1} = ht(L1),
                    {H2, T2} = ht(L2),
                    if
                        H1 < H2 ->
                            {H1, {T1, L2}};
                        H1 > H2 ->
                            F_(none, {L1, drop_less(L2, H1)});
                        true ->
                            F_(none, {T1, T2})
                    end
            end
        ),
    tl(IL).

%---------------------------------------------------------------------------------------------------
