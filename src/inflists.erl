%% @doc
%% Infinite lists realization.
%%
%% @author Alexey Rybakov

% Module name.
-module(inflists).

% Export.
-export([iterate/3, iterate/2,
         hd/1, tl/1, ht/1,
         take/2,
         repeat/1, cycle/1, seq/2]).

% We don't need hd and tl functions for lists because we use the same names for infinite lists.
-compile({no_auto_import, [hd/1, tl/1]}).

%---------------------------------------------------------------------------------------------------
% Types.
%---------------------------------------------------------------------------------------------------

% Types export.
-export_type([inflist/0]).

% Define infinite list as a record.
-record(inflist,
{
    h :: term(),
    acc :: term(),
    f :: fun((term(), term()) -> {term(), term()})
}).

% Infinite list.
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

-spec iterate(H, F) -> inflist()
      when H :: term(),
           F :: fun((H) -> H).
%% @doc
%% Create infinite list with head and iterate function (without accumulator).
iterate(H, F) when is_function(F, 1) ->
    iterate
    (
        H,
        none,
        fun(Cur_H, _) ->
            {F(Cur_H), none}
        end
    );
iterate(_, F) ->
    throw({badarg, {F, wrong_arity}}).
 
%---------------------------------------------------------------------------------------------------
% Take list elements and sublists.
%---------------------------------------------------------------------------------------------------

-spec hd(IL :: inflist()) -> term().
%% @doc
%% Infinite list head.
hd(IL) when is_record(IL, inflist) ->
    IL#inflist.h.

%---------------------------------------------------------------------------------------------------

-spec tl(IL :: inflist()) -> inflist().
%% @doc
%% Infinite list tail.
tl(#inflist{h = H, acc = Acc, f = F} = IL) ->
    {New_H, New_Acc} = F(H, Acc),
    IL#inflist{h = New_H, acc = New_Acc}.

%---------------------------------------------------------------------------------------------------

-spec ht(IL :: inflist()) -> {term(), inflist()}.
%% @doc
%% Infinite list head and tail.
ht(IL) ->
    {hd(IL), tl(IL)}.

%---------------------------------------------------------------------------------------------------

-spec take(IL :: inflist(), N :: integer()) -> list().
%% @doc
%% Take first elements of the infinite list.
take(_, N) when (N < 0) ->
    throw({badarg, N});
take(IL, N) ->
    take (IL, N, []).

-spec take(IL :: inflist(), N :: integer(), [E]) -> [E]
      when E :: term().
%% @private
%% @doc
%% Take first elements of the infinite list.
take(_, 0, R) ->
    lists:reverse(R);
take(IL, N, R) ->
    take(tl(IL), N - 1, [hd(IL) | R]).

%---------------------------------------------------------------------------------------------------
% Simple infinite lists.
%---------------------------------------------------------------------------------------------------

-spec repeat(T :: term()) -> inflist().
%% @doc
%% Haskell like repeat function.
%% T -> [T, T, ..]
repeat(T) ->
    iterate(T, fun(_) -> T end).

%---------------------------------------------------------------------------------------------------

-spec cycle(L :: list()) -> inflist().
%% @doc
%% Haskell like cycle function.
%% Constructs infinite list, containing infinite number of list L copies.
%% [A, B, C] -> [A, B, C, A, B, C, ..]
cycle([]) ->
    throw({badarg, []});
cycle([H | T]) ->
    iterate
    (
        H, T,
        fun
            (_, []) ->
                {H, T};
            (_, [Cur_H | Cur_T]) ->
                {Cur_H, Cur_T}
        end
    ).

%---------------------------------------------------------------------------------------------------

-spec seq(From :: number(), Step :: number()) -> inflist().
%% @doc
%% Arithmetic series.
seq(From, Step) ->
    iterate(From, fun(H) -> H + Step end).

%---------------------------------------------------------------------------------------------------
