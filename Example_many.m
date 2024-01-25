 % ----------------------------------------------------------------------- %
% Example of use of the Benchmark Functions                          %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
clear all;
clc;

% Objective function

%ObjFnc = 'F1';
ObjFnc = 'F2';
%ObjFnc = 'F3';
%ObjFnc = 'F4';
%ObjFnc = 'F5';
%ObjFnc = 'F6';

switch ObjFnc
    case 'F1'   
        Obj.nVar = 3;
        Obj.var_min = -5.12;
        Obj.var_max = 5.12;
        Obj.fun = @Sphere;
    case 'F2'
        Obj.nVar = 3;
        Obj.var_min = -2.48;
        Obj.var_max = 2.48;
        Obj.fun = @rosenbrock;
    case 'F3'
        Obj.nVar = 3;
        Obj.var_min = -5.12;
        Obj.var_max = 5.12;
        Obj.fun = @rastrigin;
    case 'F4'
        Obj.nVar = 10;
        Obj.var_min = -600;
        Obj.var_max = 600;
        Obj.fun = @griewank;
    case 'F5'
        Obj.nVar = 5;
        Obj.var_min = -32.768;
        Obj.var_max = 32.768;
        Obj.fun = @ackley;
    case 'F6'
        Obj.nVar = 10;
        Obj.var_min = -2.5;
        Obj.var_max = 2.5;
        Obj.fun = @beale;
end
% LEAO
[sol_final,media,desvio_padrao]=AOL1(Obj);
% REP;
% % Display info
% display('Repository fitness values are stored in REP.pos_fit');
% display('Repository particles positions are store in REP.pos');
