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
         map/2, filter/2, adj_pairs_map/2, fold/3, is_all/3, is_any/3,
         % Mathematical functions.
         add/2, inc/1, sub/2, dec/1, neg/1, mul/2, dvs/2, inv/1, ndvs/2, nrem/2,
         square/1, sqrt/1, cube/1, pow/2, npow/2,
         % More complex mathematical functions.
         pows/2, npows/2, partial_sums/1, partial_products/1, partial_avgs/1,
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

% Inifinite list - list containing infinite number of element.
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
%% Iterate function produces the second infinite list element from its head and accumulator.
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
%% Iterate function produces the second infinite list element from its head.
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
%%
%% Example:
%% <pre>
%% repeat(T) -> [T, T, T, ..]
%% </pre>
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
%%
%% Example:
%% <pre>
%% repeat([A, B, C]) -> [A, B, C, A, B, C, ..]
%% </pre>
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
%% Construct arithmetic series.
%%
%% Example:
%% <pre>
%% seq(From, Step) -> [From, From + Step, From + 2 * Step, ..]
%% </pre>
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
%%
%% Example:
%% <pre>
%% odds() -> [1, 3, 5, 7, 9, ..]
%% </pre>
odds() ->
    seq(1, 2).

%---------------------------------------------------------------------------------------------------

-spec evens() -> inflist().
%% @doc
%% Even natural numbers.
%%
%% Example:
%% <pre>
%% evens() -> [2, 4, 6, 8, 10, ..]
%% </pre>
evens() ->
    seq(2, 2).

%---------------------------------------------------------------------------------------------------

-spec seq(From :: number()) -> inflist().
%% @doc
%% Construct infinite list of naturals from given number.
%%
%% Example:
%% <pre>
%% seq(From) -> [From, From + 1, From + 2, ..]
%% </pre>
%%
%% @see naturals/1
seq(From) ->
    seq(From, 1).

%---------------------------------------------------------------------------------------------------

-spec naturals() -> inflist().
%% @doc
%% Infinite list of natural numbers.
%%
%% Example:
%% <pre>
%% naturals() -> [1, 2, 3, 4, 5, ..]
%% </pre>
naturals() ->
    seq(1).

%---------------------------------------------------------------------------------------------------

-spec naturals(From :: integer()) -> inflist().
%% @doc
%% Naturals from given number (synonym for seq).
%%
%% Example:
%% <pre>
%% seq(From) -> [From, From + 1, From + 2, ..]
%% </pre>
%%
%% @see seq/1
naturals(From) ->
    seq(From).

%---------------------------------------------------------------------------------------------------
% Geometric series.
%---------------------------------------------------------------------------------------------------

-spec geometric_series(Base :: number(), K :: number()) -> inflist().
%% @doc
%% Construct geometric series.
%%
%% Example:
%% <pre>
%% geometric_series(Base, K) -> [Base, Base * K, Base * K^2, Base * K^3, ..]
%% </pre>
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
%%
%% Example:
%% <pre>
%% power_series(X) -> [1, X, X^2, X^3, ..]
%% </pre>
power_series(X) ->
    geometric_series(1, X).

%---------------------------------------------------------------------------------------------------
% Zip/unzip functions and functors.
%---------------------------------------------------------------------------------------------------

-spec zip(IL1 :: inflist(), IL2 :: inflist()) -> inflist().
%% @doc
%% Zip two infinite lists.
%%
%% Example:
%% <pre>
%% A = [a1, a2, a3, a4, a5, ..]
%% B = [b1, b2, b3, b4, b5, ..]
%% 
%% zip(A, B) -> [{a1, b1}, {a2, b2}, {a3, b3}, {a4, b4}, {a5, b5}, ..]
%% </pre>
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
%%
%% Example:
%% <pre>
%% A = [a1, a2, a3, a4, a5, ..]
%% B = [b1, b2, b3, b4, b5, ..]
%% C = [c1, c2, c3, c4, c5, ..]
%%
%% zip_3(A, B, C) -> [{a1, b1, c1}, {a2, b2, c2}, {a3, b3, c3}, ..]
%% </pre>
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
%%
%% Example:
%% <pre>
%% A = [a1, a2, a3, a4, a5, ..]
%% B = [b1, b2, b3, b4, b5, ..]
%%
%% zipwith(A, B, Zip_F) -> [Zip_F(a1, b1), Zip_F(a2, b2), Zip_F(a3, b3), ..]
%% </pre>
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
%%
%% Example:
%% <pre>
%% IL = [{a1, b1}, {a2, b2}, {a3, b3}, {a4, b4}, {a5, b5}, ..]
%%
%% unzip(IL) -> {[a1, a2, a3, a4, a5, ..],
%%               [b1, b2, b3, b4, b5, ..]}
%% </pre>
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
%%
%% Example:
%% <pre>
%% IL = [{a1, b1, c1}, {a2, b2, c2}, {a3, b3, c3}, {a4, b4, c4}, ..]
%%
%% unzip_3(IL) -> {[a1, a2, a3, a4, ..],
%%                 [b1, b2, b3, b4, ..],
%%                 [c1, c2, c3, c4, ..]}
%% </pre>
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
%%
%% Example:
%% <pre>
%% IL = [a1, a2, a4, a4, a5, ..]
%%
%% map(IL, Map_F) -> [Map_F(a1), Map_F(a2), Map_F(a3), Map_F(a4), Map_F(a5), ..]
%% </pre>
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
%% Warning!
%% Dangerous function.
%% It can cause infinite recursion if result list is not finite. 
%%
%% Example:
%% <pre>
%% IL = [1, 2, 3, 4, 5, 6, ..]
%%
%% filter(IL, fun(X) -> X rem 2 =:= 1 end) -> [1, 3, 5, 7, ..]
%% filter(IL, fun(X) -> X =:= 0 end) -> infinite loop
%% </pre>
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
%%
%% Example:
%% <pre>
%% IL = [a1, a2, a3, a4, a5, ..]
%%
%% adj_pairs_map(IL, Map_F) -> [Map_F(a1, a2), Map_F(a2, a3), Map_F(a3, a4), ..]
%% </pre>
adj_pairs_map(IL, Map_F) ->
    zipwith(IL, tl(IL), Map_F).

%---------------------------------------------------------------------------------------------------

-spec fold(IL :: inflist(), Fold_F :: fun((term(), Fold_Acc) -> Fold_Acc), Fold_Acc) -> inflist()
      when Fold_Acc :: term().
%% @doc
%% Partial folds of sublists.
%% Note. We will never reach fold result because list is infinite.
%%
%% Example:
%% <pre>
%% IL = [a1, a2, a3, a4, a5, ..]
%%
%% fold(IL, Fold_F) -> [lists:foldl(Fold_F, Fold_Acc, [a1]),
%%                      lists:foldl(Fold_F, Fold_Acc, [a1, a2]),
%%                      lists:foldl(Fold_F, Fold_Acc, [a1, a2, a3]),
%%                      lists:foldl(Fold_F, Fold_Acc, [a1, a2, a3, a4]),
%%                      lists:foldl(Fold_F, Fold_Acc, [a1, a2, a3, a4, a5]),
%%                      ..]
%% </pre>
fold(#inflist{h = H, acc = Acc, f = F}, Fold_F, Fold_Acc) when is_function(Fold_F, 2) ->
    iterate
    (
        Fold_F(H, Fold_Acc),
        {H, Acc},
        fun(Cur_Fold_Acc, {Cur_H, Cur_Acc}) ->
            {New_H, New_Acc} = F(Cur_H, Cur_Acc),
            {Fold_F(New_H, Cur_Fold_Acc), {New_H, New_Acc}}
        end
    );
fold(IL, Fold_F, Fold_Acc) ->
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
%%
%% Example:
%% <pre>
%% A = [A1, A2, A3, A4, A5, ..]
%% B = [B1, B2, B3, B4, B5, ..]
%% V - not inflist
%%
%% add(A, B) -> [A1 + B1, A2 + B2, A3 + B3, A4 + B4, A5 + B5, ..]
%% add(A, V) -> [A1 + V, A2 + V, A3 + V, A4 + V, A5 + V, ..]
%% add(V, A) -> [V + A1, V + A2, V + A3, V + A4, V + A5, ..]
%% </pre>
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
%%
%% Example:
%% <pre>
%% A = [A1, A2, A3, A4, A5, ..]
%%
%% inc(A) -> [A1 + 1, A2 + 1, A3 + 1, A4 + 1, A5 + 1, ..]
%% </pre>
inc(IL) ->
    add(IL, 1).

%---------------------------------------------------------------------------------------------------

-spec sub(Arg, Arg) -> inflist()
      when Arg :: inflist() | term().
%% @doc
%% Sub function.
%%
%% Example:
%% <pre>
%% A = [A1, A2, A3, A4, A5, ..]
%% B = [B1, B2, B3, B4, B5, ..]
%% V - not inflist
%%
%% sub(A, B) -> [A1 - B1, A2 - B2, A3 - B3, A4 - B4, A5 - B5, ..]
%% sub(A, V) -> [A1 - V, A2 - V, A3 - V, A4 - V, A5 - V, ..]
%% sub(V, A) -> [V - A1, V - A2, V - A3, V - A4, V - A5, ..]
%% </pre>
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
%%
%% Example:
%% <pre>
%% A = [A1, A2, A3, A4, A5, ..]
%%
%% dec(A) -> [A1 - 1, A2 - 1, A3 - 1, A4 - 1, A5 - 1, ..]
%% </pre>
dec(IL) ->
    sub(IL, 1).

%---------------------------------------------------------------------------------------------------

-spec neg(IL :: inflist()) -> inflist().
%% @doc
%% Negate infinite list.
%%
%% Example:
%% <pre>
%% A = [A1, A2, A3, A4, A5, ..]
%%
%% neg(A) -> [-A1, -A2, -A3, -A4, -A5, ..]
%% </pre>
neg(IL) ->
    map(IL, fun(X) -> -X end).

%---------------------------------------------------------------------------------------------------

-spec mul(Arg, Arg) -> inflist()
      when Arg :: inflist() | term().
%% @doc
%% Multiplication.
%%
%% Example:
%% <pre>
%% A = [A1, A2, A3, A4, A5, ..]
%% B = [B1, B2, B3, B4, B5, ..]
%% V - not inflist
%%
%% mul(A, B) -> [A1 * B1, A2 * B2, A3 * B3, A4 * B4, A5 * B5, ..]
%% mul(A, V) -> [A1 * V, A2 * V, A3 * V, A4 * V, A5 * V, ..]
%% mul(V, A) -> [V * A1, V * A2, V * A3, V * A4, V * A5, ..]
%% </pre>
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
%%
%% Example:
%% <pre>
%% A = [A1, A2, A3, A4, A5, ..]
%% B = [B1, B2, B3, B4, B5, ..]
%% V - not inflist
%%
%% dvs(A, B) -> [A1 / B1, A2 / B2, A3 / B3, A4 / B4, A5 / B5, ..]
%% dvs(A, V) -> [A1 / V, A2 / V, A3 / V, A4 / V, A5 / V, ..]
%% dvs(V, A) -> [V / A1, V / A2, V / A3, V / A4, V / A5, ..]
%% </pre>
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
%% Infinite list of inverted values.
%%
%% Example:
%% <pre>
%% A = [A1, A2, A3, A4, A5, ..]
%%
%% inv(A) -> [1 / A1, 1 / A2, 1 / A3, 1 / A4, 1 / A5, ..]
%% </pre>
inv(IL) ->
    map(IL, fun(X) -> 1 / X end).

%---------------------------------------------------------------------------------------------------

-spec ndvs(Arg, Arg) -> inflist()
      when Arg :: inflist() | term().
%% @doc
%% Division without remainder (for natural numbers).
%%
%% Example:
%% <pre>
%% A = [A1, A2, A3, A4, A5, ..]
%% B = [B1, B2, B3, B4, B5, ..]
%% V - not inflist
%%
%% ndvs(A, B) -> [A1 div B1, A2 div B2, A3 div B3, A4 div B4, A5 div B5, ..]
%% ndvs(A, V) -> [A1 div V, A2 div V, A3 div V, A4 div V, A5 div V, ..]
%% ndvs(V, A) -> [V div A1, V div A2, V div A3, V div A4, V div A5, ..]
%% </pre>
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
%%
%% Example:
%% <pre>
%% A = [A1, A2, A3, A4, A5, ..]
%% B = [B1, B2, B3, B4, B5, ..]
%% V - not inflist
%%
%% nrem(A, B) -> [A1 rem B1, A2 rem B2, A3 rem B3, A4 rem B4, A5 rem B5, ..]
%% nrem(A, V) -> [A1 rem V, A2 rem V, A3 rem V, A4 rem V, A5 rem V, ..]
%% nrem(V, A) -> [V rem A1, V rem A2, V rem A3, V rem A4, V rem A5, ..]
%% </pre>
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
%%
%% Example:
%% <pre>
%% A = [A1, A2, A3, A4, A5, ..]
%%
%% square(A) -> [A1^2, A2^2, A3^2, A4^2, A5^2, ..]
%% </pre>
square(IL) ->
    mul(IL, IL).

%---------------------------------------------------------------------------------------------------

-spec sqrt(IL :: inflist()) -> inflist().
%% @doc
%% Square root of infinite list.
%%
%% Example:
%% <pre>
%% A = [A1, A2, A3, A4, A5, ..]
%%
%% sqrt(A) -> [sqrt(A1), sqrt(A2), sqrt(A3), sqrt(A4), sqrt(A5), ..]
%% </pre>
sqrt(IL) ->
    map(IL, fun(X) -> math:sqrt(X) end).

%---------------------------------------------------------------------------------------------------

-spec cube(IL :: inflist()) -> inflist().
%% @doc
%% Cube of infinite list.
%%
%% Example:
%% <pre>
%% A = [A1, A2, A3, A4, A5, ..]
%%
%% cube(A) -> [A1^3, A2^3, A3^3, A4^3, A5^3, ..]
%% </pre>
cube(IL) ->
    mul(square(IL), IL).

%---------------------------------------------------------------------------------------------------

-spec pow(IL :: inflist(), P :: number()) -> inflist().
%% @doc
%% Power of infinite list.
%%
%% Example:
%% <pre>
%% A = [A1, A2, A3, A4, A5, ..]
%%
%% pow(A, P) -> [A1^P, A2^P, A3^P, A4^P, A5^P, ..]
%% </pre>
pow(IL, P) ->
    map(IL, fun(X) -> math:pow(X, P) end).

%---------------------------------------------------------------------------------------------------

-spec npow(IL :: inflist(), P :: number()) -> inflist().
%% @doc
%% Power of infinite list (for natural numbers).
%%
%% Example:
%% <pre>
%% A = [A1, A2, A3, A4, A5, ..]
%%
%% npow(A, P) -> [A1^P, A2^P, A3^P, A4^P, A5^P, ..]
%% </pre>
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
%%
%% Example:
%% <pre>
%% A = [A1, A2, A3, A4, A5, ..]
%%
%% pows = [N^A1, N^A2, N^A3, N^A4, N^A5, ..]
%% </pre>
pows(N, IL) ->
    map(IL, fun(X) -> math:pow(N, X) end).

%---------------------------------------------------------------------------------------------------

-spec npows(N :: number(), IL :: inflist()) -> inflist().
%% @doc
%% Return natural powers of given number (powers are taken from IL).
%%
%% Example:
%% <pre>
%% A = [A1, A2, A3, A4, A5, ..]
%%
%% pows = [N^A1, N^A2, N^A3, N^A4, N^A5, ..]
%% </pre>
npows(N, IL) ->
    map
    (
        IL,
        fun(P) ->
            N_Pow =
                fun
                    % We need x^0 for sequences 1, x, x^2, ..
                    NP_(_, 0) ->
                        1;
                    NP_(N_, P_) when (P_ > 0) ->
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
%%
%% Example:
%% <pre>
%% A = [A1, A2, A3, A4, A5, ..]
%%
%% partial_sums(A) -> [A1, 
%%                     A1 + A2,
%%                     A1 + A2 + A3,
%%                     A1 + A2 + A3 + A4,
%%                     A1 + A2 + A3 + A4 + A5,
%%                     ..]
%% </pre>
partial_sums(IL) ->
    fold(IL, fun(X, Y) -> X + Y end, 0).

%---------------------------------------------------------------------------------------------------

-spec partial_products(IL :: inflist()) -> inflist().
%% @doc
%% Partial products.
%%
%% Example:
%% <pre>
%% A = [A1, A2, A3, A4, A5, ..]
%%
%% partial_products(A) -> [A1, 
%%                         A1 * A2,
%%                         A1 * A2 * A3,
%%                         A1 * A2 * A3 * A4,
%%                         A1 * A2 * A3 * A4 * A5,
%%                         ..]
%% </pre>
partial_products(IL) ->
    fold(IL, fun(X, Y) -> X * Y end, 1).

%---------------------------------------------------------------------------------------------------

-spec partial_avgs(IL :: inflist()) -> inflist().
%% @doc
%% Average values of infinite list.
%%
%% Example:
%% <pre>
%% A = [A1, A2, A3, A4, A5, ..]
%%
%% partial_avgs(A) -> [A1, 
%%                     (A1 + A2) / 2,
%%                     (A1 + A2 + A3) / 3,
%%                     (A1 + A2 + A3 + A4) / 4,
%%                     (A1 + A2 + A3 + A4 + A5) / 5,
%%                     ..]
%% </pre>
partial_avgs(IL) ->
    dvs(partial_sums(IL), naturals()).

%---------------------------------------------------------------------------------------------------

-spec dirichlet_series(S :: number()) -> inflist().
%% @doc
%% Dirichlet series.
%%
%% Example:
%% <pre>
%%                              1      1      1      1 
%% dirichlet_series(S) -> [1, -----, -----, -----, -----, ..]
%%                             2^S    3^S    4^S    5^S
%% </pre>
dirichlet_series(S) ->
    pow(harmonic_series(), S).

%---------------------------------------------------------------------------------------------------

-spec dirichlet_series(IL :: inflist(), S :: number()) -> inflist().
%% @doc
%% Dirichlet series of base inflist IL and pow degree S.
%%
%% Example:
%% <pre>
%% A = [A1, A2, A3, A4, A5, ..]
%%
%%                               A2     A3     A4     A5 
%% dirichlet_series(S) -> [A1, -----, -----, -----, -----, ..]
%%                              2^S    3^S    4^S    5^S
%% </pre>
dirichlet_series(IL, S) ->
    mul(IL, dirichlet_series(S)).

%---------------------------------------------------------------------------------------------------

-spec sign_alternate(IL :: inflist()) -> inflist().
%% @doc
%% Alternate sign of infinite list.
%% Odd position elements are unchanged, even position elements are negated.
%%
%% Example:
%% <pre>
%% A = [A1, A2, A3, A4, A5, ..]
%%
%% sign_alternate(A) -> [A1, -A2, A3, -A4, A5, ..]
%% </pre>
sign_alternate(IL) ->
    mul(IL, grundy_series()).

%---------------------------------------------------------------------------------------------------
% Some usefull sequences.
%---------------------------------------------------------------------------------------------------

-spec fib() -> inflist().
%% @doc
%% Fibonacci numbers.
%%
%% Example:
%% <pre>
%% fib() -> [1, 1, 2, 3, 5, 8, 13, 21, ..]
%% </pre>
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

-spec trib() -> inflist().
%% @doc
%% Tribonacci numbers.
%%
%% Example:
%% <pre>
%% trib() -> [0, 1, 1, 2, 4, 7, 13, 24, 44, ..]
%% </pre>
trib() ->
    iterate
    (
        0,
        {1, 1},
        fun(H, {A1, A2}) ->
            {A1, {A2, H + A1 + A2}}
        end
    ).

%---------------------------------------------------------------------------------------------------

-spec harmonic_series() -> inflist().
%% @doc
%% Harmonic series.
%%
%% Example:
%% <pre>
%%                           1    1    1    1    1    1  
%% harmonic_series() -> [1, ---, ---, ---, ---, ---, ---, ..]
%%                           2    3    4    5    6    7
%% </pre>
harmonic_series() ->
    inv(naturals()).

%---------------------------------------------------------------------------------------------------

-spec anharmonic_series() -> inflist().
%% @doc
%% Anharmonic series (Leibniz series).
%%
%% Example:
%% <pre>
%%                               1    1      1    1   
%% anharmonic_series() -> [1, - ---, ---, - ---, ---, ..]
%%                               3    5      7    9
%% </pre>
anharmonic_series() ->
    sign_alternate(inv(odds())).

%---------------------------------------------------------------------------------------------------

-spec grundy_series() -> inflist().
%% @doc
%% Grundy series (see Patrick Carmelo Grundy).
%%
%% Example:
%% <pre>
%% grundy_series() -> [1, -1, 1, -1, 1, -1, ..]
%% </pre>
grundy_series() ->
    cycle([1, -1]).

%---------------------------------------------------------------------------------------------------

-spec facts() -> inflist().
%% @doc
%% Factorials series (starts with 0! = 1).
%%
%% Example:
%% <pre>
%% facts() -> [1, 1, 2!, 3!, 4!, 5!, ..]
%% </pre>
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
%%
%% Example:
%% <pre>
%%                        1    1    1
%% inv_facts() -> [1, 1, ---, ---, ---, ..]
%%                        2!   3!   4!
%% </pre>
inv_facts() ->
    inv(facts()).

%---------------------------------------------------------------------------------------------------

-spec squares() -> inflist().
%% @doc
%% Natural squares.
%%
%% Example:
%% <pre>
%% squares() -> [1, 4, 9, 16, 25, ..]
%% </pre>
squares() ->
    square(naturals()).

%---------------------------------------------------------------------------------------------------

-spec sqrts() -> inflist().
%% @doc
%% Square roots of naturals.
%%
%% Example:
%% <pre>
%% sqrts() -> [1, sqrt(2), sqrt(3), 2, sqrt(5), ..]
%% </pre>
sqrts() ->
    sqrt(naturals()).

%---------------------------------------------------------------------------------------------------

-spec cubes() -> inflist().
%% @doc
%% Natural cubes.
%%
%% Example:
%% <pre>
%% cubes() -> [1, 8, 27, 64, 125, ..]
%% </pre>
cubes() ->
    cube(naturals()).

%---------------------------------------------------------------------------------------------------

-spec triangulars() -> inflist().
%% @doc
%% Triangulars numbers.
%%
%% Example:
%% <pre>
%% triangulars() -> [1, 3, 6, 10, 15, 21, ..]
%% </pre>
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
%%
%% Example:
%% <pre>
%% primes() -> [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, ..]
%% </pre>
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
%%
%% Example:
%% <pre>
%% IL = [a1, a2, a3, a4, a5, ..]
%% L = [b1, b2, b3]
%% T - not inflist
%% 
%% concat(L, IL) -> [b1, b2, b3, a1, a2, a3, a4, a5, ..]
%% concat(T, IL) -> [T, a1, a2, a3, a4, a5, ..]
%% </pre>
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
%%
%% Example:
%% <pre>
%% A = [a1, a2, a3, a4, a5, ..]
%%
%% sparse(A, 0) -> A
%% sparse(A, 1) -> [a1, a3, a5, ...]
%% sparse(A, 2) -> [a1, a4, a7, ...]
%% </pre>
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
%%
%% Example:
%% <pre>
%% A = [a1, a2, a4, a4, a5, ..]
%%
%% odds(A) -> [a1, a3, a5, ..]
%% </pre>
odds(IL) ->
    sparse(IL, 1).

%---------------------------------------------------------------------------------------------------

-spec evens(IL :: inflist()) -> inflist().
%% @doc
%% Even elements of list.
%%
%% Example:
%% <pre>
%% A = [a1, a2, a3, a4, a5, ..]
%%
%% evens(A) -> [a2, a4, a6, ..]
%% </pre>
evens(IL) ->
    sparse(tl(IL), 1).

%---------------------------------------------------------------------------------------------------
% Merge/unmerge.
%---------------------------------------------------------------------------------------------------

-spec merge(IL1 :: inflist(), IL2 :: inflist()) -> inflist().
%% @doc
%% Merge two infinite lists.
%%
%% Example:
%% <pre>
%% A = [a1, a2, a3, a4, a5, ..]
%% B = [b1, b2, b3, b4, b5, ..]
%%
%% merge(A, B) -> [a1, b1, a2, b2, a3, b3, a4, b4, ..]
%% </pre>
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
%%
%% Example:
%% <pre>
%% A = [a1, a2, a3, a4, a4, ..]
%%
%% unmerge(A) -> {[a1, a3, a5, ..], [a2, a4, a6, ..]}
%% </pre>
unmerge(IL) ->
    {odds(IL), evens(IL)}.

%---------------------------------------------------------------------------------------------------
% Taylor series.
%---------------------------------------------------------------------------------------------------

-spec taylor_exp(X :: number()) -> inflist().
%% @doc
%% Taylor series of e^x for (-inf, inf).
%% <pre>
%%                       X    X^2    X^3
%% taylor_exp(X) -> [1, ---, -----, -----, ..]
%%                       1!    2!     3!
%% </pre>
taylor_exp(X) ->
    dvs(power_series(X), facts()).

%---------------------------------------------------------------------------------------------------

-spec taylor_lnxp1(X :: number()) -> inflist().
%% @doc
%% Taylor series of ln(1 + x) for (-1, 1].
%% <pre>
%%                                       X^2    X^3      X^4
%% taylor_lnxp1(X) -> ln(x + 1) = [X, - -----, -----, - -----, ..]
%%                                        2      3        4
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
%%                        X^3    X^5      X^7
%% taylor_sin(X) = [X, - -----, -----, - -----, ..]
%%                         3!     5!       7!
%% </pre>
taylor_sin(X) ->
    sign_alternate(evens(taylor_exp(X))).

%---------------------------------------------------------------------------------------------------

-spec taylor_cos(X :: number()) -> inflist().
%% @doc
%% Taylor series of cos(x) for (-inf, inf).
%% <pre>
%%                        X^2    X^4      X^6
%% taylor_cos(X) = [1, - -----, -----, - -----, ..]
%%                         2!     4!       6!
%% </pre>
taylor_cos(X) ->
    sign_alternate(odds(taylor_exp(X))).

%---------------------------------------------------------------------------------------------------

-spec taylor_arctg(X :: number()) -> inflist().
%% @doc
%% Taylor series of arctg(x) for (-1, 1).
%% <pre>
%%                          X^3    X^5      X^7
%% taylor_arctg(X) = [X, - -----, -----, - -----, ..]
%%                           3      5        7
%% </pre>
taylor_arctg(X) when (abs(X) >= 1) ->
    throw({badarg, X});
taylor_arctg(X) ->
    sign_alternate(dvs(evens(power_series(X)), naturals())).

%---------------------------------------------------------------------------------------------------
% Sets operations.
%---------------------------------------------------------------------------------------------------

-spec mono_merge(IL1 :: inflist(), IL2 :: inflist()) -> inflist().
%% @doc
%% Merge of two monotonous lists.
%%
%% Example:
%% <pre>
%% A = [1, 3, 5, 7, 9, ..]
%% B = [1, 4, 9, 16, 25, ..]
%%
%% mono_merge(A, B) -> [1, 1, 3, 4, 5, 7, 9, 9, 11, 13, ..]
%% </pre>
mono_merge(IL1, IL2) ->
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

-spec mono_unique(IL :: inflist()) -> inflist().
%% @doc
%% Take only unique elements of infinite list.
%% Warning.
%% This is dangerous function, because it can lead to infinite loop
%% when list has tail consisting of the same element.
%%
%% Example:
%% <pre>
%% mono_unique([1, 1, 2, 3, 3, 3, 4, 5, ..]) -> [1, 2, 3, 4, 5, ..]
%% mono_unique([a, a, a, a, a, ..]) -> infinite loop
%% </pre>
mono_unique(#inflist{h = H, acc = Acc, f = F}) ->
    iterate
    (
        H,
        Acc,
        fun
            _F(Cur_H, Cur_Acc) ->
                {Next_H, Next_Acc} = F(Cur_H, Cur_Acc),
                if
                    Next_H =:= Cur_H ->
                        _F(Next_H, Next_Acc);
                    true ->
                        {Next_H, Next_Acc}
                end
        end
    ).
  
%---------------------------------------------------------------------------------------------------

-spec mono_union(IL1 :: inflist(), IL2 :: inflist()) -> inflist().
%% @doc
%% Union of two monotolnous lists.
%%
%% Example:
%% <pre>
%% A = [1, 3, 5, 7, 9, ..]
%% B = [1, 4, 9, 16, 25, ..]
%%
%% mono_union(A, B) -> [1, 3, 4, 5, 7, 9, 11, 13, ..]
%% </pre>
mono_union(IL1, IL2) ->
    mono_unique(mono_merge(IL1, IL2)).

%---------------------------------------------------------------------------------------------------

-spec mono_intersection(IL1 :: inflist(), IL2 :: inflist()) -> inflist().
%% @doc
%% Intersection of two monotonous lists.
%% Warning.
%% Dangerous function, because it can lead to infinite loop.
%%
%% Example:
%% <pre>
%% A = [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, ..]
%% B = [1, 4, 9, 16, 25, ..]
%%
%% mono_intersection(A, B) -> [1, 9, 25, ..]
%% mono_intersection(odds(), evens()) -> infinite loop
%% </pre>
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
%% Complement of monotonous infinite lists.
%% All element of the first list which are not elements of the second list.
%%
%% Example:
%% <pre>
%% A = [1, 3, 5, 7, 9, 11, 13, ..]
%% B = [1, 4, 9, 16, 25, ..]
%%
%% mono_complement(A, B) -> [3, 5, 7, 11, 13, 15, 17, 19, 21, 23, 27, 29, ..]
%% </pre>
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
