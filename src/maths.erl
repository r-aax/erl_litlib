%% @doc
%% Mathematical functions realization.
%%
%% @author Alexey Rybakov

% Module name.
-module(maths).

% Export.
-export([randf01/0, randf/1, randf/2, randi/1, randi/2, rand_dev/2,
         factorization/1]).

%---------------------------------------------------------------------------------------------------
% Random numbers generation.
%---------------------------------------------------------------------------------------------------

-spec randf01() -> float().
%% @doc
%% Generate random float value in [0.0, 1.0] segment.
randf01() ->
    rand:uniform().

%---------------------------------------------------------------------------------------------------

-spec randf(H :: float()) -> float().
%% @doc
%% Generate random float value in [0.0, H] segment.
randf(H) ->
    H * randf01().

%---------------------------------------------------------------------------------------------------

-spec randf(L :: float(), H :: float()) -> float().
%% @doc
%% Generate random float value in [L, H] segment.
randf(L, H) when (H >= L) ->
    L + randf(H - L).

%---------------------------------------------------------------------------------------------------

-spec randi(H :: integer()) -> integer().
%% @doc
%% Generate random integer value in [0, H] segment.
randi(H) ->
    rand:uniform(H + 1) - 1.

%---------------------------------------------------------------------------------------------------

-spec randi(L :: integer(), H :: integer()) -> integer().
%% @doc
%% Generate random integer value in [L, H] segment.
randi(L, H) when (H >= L) ->
    L + randi(H - L).

%---------------------------------------------------------------------------------------------------

-spec rand_dev(V :: float(), D :: float()) -> float().
%% @doc
%% Generate random value from [V - abs(V * D), V + abs(V* D)] segment.
rand_dev(V, D) ->
    Dev = abs(V * D),
    randf(V - Dev, V + Dev).

%---------------------------------------------------------------------------------------------------
% Factorization.
%---------------------------------------------------------------------------------------------------

-spec factorization(N :: integer()) -> [integer()].
%% @doc
%% Factorization of natural number.
factorization(N) ->
    factorization(N, 2, []).

-spec factorization(N :: integer(), F :: integer(), Fs :: [integer()]) -> [integer()].
%% @doc
%% Factorization of natural number.
factorization(1, _, []) ->
    [1];
factorization(1, _, Fs) ->
    Fs;
factorization(N, F, Fs) ->

    if
        % If F > sqrt(N) then there is no other factors.
        F * F > N ->
            [N | Fs];

        % If F divides N.
        N rem F =:= 0 ->
            factorization(N div F, F, [F | Fs]);

        % Try next factor.
        % Actually we have to try here next prime number.
        true ->
            factorization(N, F + 1, Fs)
    end.

%---------------------------------------------------------------------------------------------------
