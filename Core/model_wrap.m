function [J,sol,rout] = model_wrap(pars,data)
%function model_wrap.m
%Input: model parameters and data structure
%Output: least squares cost, solution and model residual
%Uses: model_sol.m solving the ODEs
%Description: Wrapper integrates estimated parameters into the total list
%of parameters.

ALLPARS = data.gpars.ALLPARS;
INDMAP  = data.gpars.INDMAP;

tpars = ALLPARS;
tpars(INDMAP') = pars;

[sol,rout,J] = model_sol(tpars,data);