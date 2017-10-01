%% @doc
%% Tests for infinite lists.
%%
%% @author Alexey Rybakov

% Module name.
-module(inflists_test).

% Unit-test using.
-include_lib("eunit/include/eunit.hrl").

% We don't need hd and tl functions for lists because we use the same names for infinite lists.
-compile({no_auto_import, [hd/1, tl/1]}).

% Functions import.
-import(inflists,
	    [take/2,
	     repeat/1, cycle/1, seq/2]).

%---------------------------------------------------------------------------------------------------
% Simple infinite lists.
%---------------------------------------------------------------------------------------------------

-spec repeat_test() -> ok.
%% @doc
%% Function repeat test.
repeat_test() ->
	?assertEqual([x, x, x, x, x], take(repeat(x), 5)).

%---------------------------------------------------------------------------------------------------

-spec cycle_test() -> ok.
%% @doc
%% Function cycle test.
cycle_test() ->
	?assertEqual([a, b, c, a, b, c, a, b, c, a], take(cycle([a, b, c]), 10)).

%---------------------------------------------------------------------------------------------------

-spec seq_test() -> ok.
%% @doc
%% Function seq test.
seq_test() ->
	?assertEqual([5, 6, 7, 8, 9], take(seq(5, 1), 5)).

%---------------------------------------------------------------------------------------------------