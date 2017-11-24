%% @doc
%% Mathematical functions realization.
%%
%% @author Alexey Rybakov

% Module name.
-module(maths).

% Export.
-export([factorization/1]).

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
