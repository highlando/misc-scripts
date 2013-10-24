function [fvt] = ftWav(ft,t)
%evaluates the time component of the control corresponding to the specified
%basis function no. ft
%ftLin uses as basis the (hirarchical) ordered subspace of piecewise linear
%functions
%ATTENTION: ftLin uses te, which has to be defined accordingly to runIoMap


te = 0.1;
tif = sqrt(1/te);
fvt = 0;            %local basis functions - in most cases f =0
sr2 = sqrt(2);
lf0 =    0.5*te;
lf1 =   0.25*te; 
lf2=   0.125*te; 
lf3=  0.0625*te;    % factors for the definition of the local basis functions
lf4= 0.03125*te;
lf5=0.015625*te;

if ft > 64
    fprintf(1,'Warning: only less than 64 wavelet basis functions implemented \n')
    return
end

switch ft
%% level 1
    case 1
        fvt = 1;
%% level 2
    case 2
       if t >= 0*lf0 && t <= 1*lf0 
           fvt = 1;
       elseif t >= 1*lf0 && t <= 2*lf0 
           fvt = -1;
       end
%% level 3
    case 3 
       if t >= 0*lf1 && t <= 1*lf1 
           fvt = sr2;
       elseif t >= 1*lf1 && t <= 2*lf1 
           fvt = -sr2;
       end
    case 4 
       if t >= 2*lf1 && t <= 3*lf1 
           fvt = sr2;
       elseif t >= 3*lf1 && t <= 4*lf1 
           fvt = -sr2;
       end  
%% level 4
    case 5
       if t >= 0*lf2 && t <= 1*lf2 
           fvt = 2;
       elseif t >= 1*lf2 && t <= 2*lf2 
           fvt = -2;
       end
    case 6
       if t >= 2*lf2 && t <= 3*lf2 
           fvt = 2;
       elseif t >= 3*lf2 && t <= 4*lf2 
           fvt = -2;
       end 
    case 7
       if t >= 4*lf2 && t <= 5*lf2 
           fvt = 2;
       elseif t >= 5*lf2 && t <= 6*lf2 
           fvt = -2;
       end
    case 8
       if t >= 6*lf2 && t <= 7*lf2 
           fvt = 2;
       elseif t >= 7*lf2 && t <= 8*lf2 
           fvt = -2;
       end
%% level 5
    case 9
       if t >= 0*lf3 && t <= 1*lf3 
           fvt = 2*sr2;
       elseif t >= 1*lf3 && t <= 2*lf3 
           fvt = -2*sr2;
       end
    case 10
       if t >= 2*lf3 && t <= 3*lf3 
           fvt = 2*sr2;
       elseif t >= 3*lf3 && t <= 4*lf3 
           fvt = -2*sr2;
       end 
    case 11
       if t >= 4*lf3 && t <= 5*lf3 
           fvt = 2*sr2;
       elseif t >= 5*lf3 && t <= 6*lf3 
           fvt = -2*sr2;
       end
    case 12
       if t >= 6*lf3 && t <= 7*lf3 
           fvt = 2*sr2;
       elseif t >= 7*lf3 && t <= 8*lf3 
           fvt = -2*sr2;
       end
    case 13
       if t >= 8*lf3 && t <= 9*lf3 
           fvt = 2*sr2;
       elseif t >= 9*lf3 && t <= 10*lf3 
           fvt = -2*sr2;
       end
    case 14
       if t >= 10*lf3 && t <= 11*lf3 
           fvt = 2*sr2;
       elseif t >= 11*lf3 && t <= 12*lf3 
           fvt = -2*sr2;
       end 
    case 15   
       if t >= 12*lf3 && t <= 13*lf3 
           fvt = 2*sr2;
       elseif t >= 13*lf3 && t <= 14*lf3 
           fvt = -2*sr2;
       end
    case 16
       if t >= 14*lf3 && t <= 15*lf3 
           fvt = 2*sr2;
       elseif t >= 15*lf3 && t <= 16*lf3 
           fvt = -2*sr2;
       end
%% level 6
    case 17
       if t >= 0*lf4 && t <= 1*lf4 
           fvt = 4;
       elseif t >= 1*lf4 && t <= 2*lf4 
           fvt = -4;
       end
    case 18
       if t >= 2*lf4 && t <= 3*lf4 
           fvt = 4;
       elseif t >= 3*lf4 && t <= 4*lf4 
           fvt = -4;
       end 
    case 19
       if t >= 4*lf4 && t <= 5*lf4 
           fvt = 4;
       elseif t >= 5*lf4 && t <= 6*lf4 
           fvt = -4;
       end
    case 20
       if t >= 6*lf4 && t <= 7*lf4 
           fvt = 4;
       elseif t >= 7*lf4 && t <= 8*lf4 
           fvt = -4;
       end
    case 21
       if t >= 8*lf4 && t <= 9*lf4 
           fvt = 4;
       elseif t >= 9*lf4 && t <= 10*lf4 
           fvt = -4;
       end
    case 22
       if t >= 10*lf4 && t <= 11*lf4 
           fvt = 4;
       elseif t >= 11*lf4 && t <= 12*lf4 
           fvt = -4;
       end 
    case 23   
       if t >= 12*lf4 && t <= 13*lf4 
           fvt = 4;
       elseif t >= 13*lf4 && t <= 14*lf4 
           fvt = -4;
       end
    case 24
       if t >= 14*lf4 && t <= 15*lf4 
           fvt = 4;
       elseif t >= 15*lf4 && t <= 16*lf4 
           fvt = -4;
       end
    case 25
       if t >= 16*lf4 && t <= 17*lf4 
           fvt = 4;
       elseif t >= 17*lf4 && t <= 18*lf4 
           fvt = -4;
       end
    case 26
       if t >= 18*lf4 && t <= 19*lf4 
           fvt = 4;
       elseif t >= 19*lf4 && t <= 20*lf4 
           fvt = -4;
       end 
    case 27
       if t >= 20*lf4 && t <= 21*lf4 
           fvt = 4;
       elseif t >= 21*lf4 && t <= 22*lf4 
           fvt = -4;
       end
    case 28
       if t >= 22*lf4 && t <= 23*lf4 
           fvt = 4;
       elseif t >= 23*lf4 && t <= 24*lf4 
           fvt = -4;
       end
    case 29
       if t >= 24*lf4 && t <= 25*lf4 
           fvt = 4;
       elseif t >= 25*lf4 && t <= 26*lf4 
           fvt = -4;
       end
    case 30
       if t >= 26*lf4 && t <= 27*lf4 
           fvt = 4;
       elseif t >= 27*lf4 && t <= 28*lf4 
           fvt = -4;
       end 
    case 31   
       if t >= 28*lf4 && t <= 29*lf4 
           fvt = 4;
       elseif t >= 29*lf4 && t <= 30*lf4 
           fvt = -4;
       end
    case 32
       if t >= 30*lf4 && t <= 31*lf4 
           fvt = 4;
       elseif t >= 31*lf4 && t <= 32*lf4 
           fvt = 4;                    %% changed to reduce aliasing error
       end
%% level 8
    case 33
       if t >= 0*lf5 && t <= 1*lf5 
           fvt = sr2*4;
       elseif t >= 1*lf5 && t <= 2*lf5 
           fvt = -sr2*4;
       end
    case 34
       if t >= 2*lf5 && t <= 3*lf5 
           fvt = sr2*4;
       elseif t >= 3*lf5 && t <= 4*lf5 
           fvt = -sr2*4;
       end 
    case 35
       if t >= 4*lf5 && t <= 5*lf5 
           fvt = sr2*4;
       elseif t >= 5*lf5 && t <= 6*lf5 
           fvt = -sr2*4;
       end
    case 36
       if t >= 6*lf5 && t <= 7*lf5 
           fvt = sr2*4;
       elseif t >= 7*lf5 && t <= 8*lf5 
           fvt = -sr2*4;
       end
    case 37
       if t >= 8*lf5 && t <= 9*lf5 
           fvt = sr2*4;
       elseif t >= 9*lf5 && t <= 10*lf5 
           fvt = -sr2*4;
       end
    case 38
       if t >= 10*lf5 && t <= 11*lf5 
           fvt = sr2*4;
       elseif t >= 11*lf5 && t <= 12*lf5 
           fvt = -sr2*4;
       end 
    case 39   
       if t >= 12*lf5 && t <= 13*lf5 
           fvt = sr2*4;
       elseif t >= 13*lf5 && t <= 14*lf5 
           fvt = -sr2*4;
       end
    case 40
       if t >= 14*lf5 && t <= 15*lf5 
           fvt = sr2*4;
       elseif t >= 15*lf5 && t <= 16*lf5 
           fvt = -sr2*4;
       end
    case 41
       if t >= 16*lf5 && t <= 17*lf5 
           fvt = sr2*4;
       elseif t >= 17*lf5 && t <= 18*lf5 
           fvt = -sr2*4;
       end
    case 42
       if t >= 18*lf5 && t <= 19*lf5 
           fvt = sr2*4;
       elseif t >= 19*lf5 && t <= 20*lf5 
           fvt = -sr2*4;
       end 
    case 43
       if t >= 20*lf5 && t <= 21*lf5 
           fvt = sr2*4;
       elseif t >= 21*lf5 && t <= 22*lf5 
           fvt = -sr2*4;
       end
    case 44
       if t >= 22*lf5 && t <= 23*lf5 
           fvt = sr2*4;
       elseif t >= 23*lf5 && t <= 24*lf5 
           fvt = -sr2*4;
       end
    case 45
       if t >= 24*lf5 && t <= 25*lf5 
           fvt = sr2*4;
       elseif t >= 25*lf5 && t <= 26*lf5 
           fvt = -sr2*4;
       end
    case 46
       if t >= 26*lf5 && t <= 27*lf5 
           fvt = sr2*4;
       elseif t >= 27*lf5 && t <= 28*lf5 
           fvt = -sr2*4;
       end 
    case 47   
       if t >= 28*lf5 && t <= 29*lf5 
           fvt = sr2*4;
       elseif t >= 29*lf5 && t <= 30*lf5 
           fvt = -sr2*4;
       end
    case 48
       if t >= 30*lf5 && t <= 31*lf5 
           fvt = sr2*4;
       elseif t >= 31*lf5 && t <= 32*lf5 
           fvt = -sr2*4;
       end
     case 49
       if t >= 0*lf5+0.5*te && t <= 1*lf5+0.5*te 
           fvt = sr2*4;
       elseif t >= 1*lf5+0.5*te && t <= 2*lf5+0.5*te 
           fvt = -sr2*4;
       end
    case 50
       if t >= 2*lf5+0.5*te && t <= 3*lf5+0.5*te 
           fvt = sr2*4;
       elseif t >= 3*lf5+0.5*te && t <= 4*lf5+0.5*te 
           fvt = -sr2*4;
       end 
    case 51
       if t >= 4*lf5+0.5*te && t <= 5*lf5+0.5*te 
           fvt = sr2*4;
       elseif t >= 5*lf5+0.5*te && t <= 6*lf5+0.5*te 
           fvt = -sr2*4;
       end
    case 52
       if t >= 6*lf5+0.5*te && t <= 7*lf5+0.5*te 
           fvt = sr2*4;
       elseif t >= 7*lf5+0.5*te && t <= 8*lf5+0.5*te 
           fvt = -sr2*4;
       end
    case 53
       if t >= 8*lf5+0.5*te && t <= 9*lf5+0.5*te 
           fvt = sr2*4;
       elseif t >= 9*lf5+0.5*te && t <= 10*lf5+0.5*te 
           fvt = -sr2*4;
       end
    case 54
       if t >= 10*lf5+0.5*te && t <= 11*lf5+0.5*te 
           fvt = sr2*4;
       elseif t >= 11*lf5+0.5*te && t <= 12*lf5+0.5*te 
           fvt = -sr2*4;
       end 
    case 55   
       if t >= 12*lf5+0.5*te && t <= 13*lf5+0.5*te 
           fvt = sr2*4;
       elseif t >= 13*lf5+0.5*te && t <= 14*lf5+0.5*te 
           fvt = -sr2*4;
       end
    case 56
       if t >= 14*lf5+0.5*te && t <= 15*lf5+0.5*te 
           fvt = sr2*4;
       elseif t >= 15*lf5+0.5*te && t <= 16*lf5+0.5*te 
           fvt = -sr2*4;
       end
    case 57
       if t >= 16*lf5+0.5*te && t <= 17*lf5+0.5*te 
           fvt = sr2*4;
       elseif t >= 17*lf5+0.5*te && t <= 18*lf5+0.5*te 
           fvt = -sr2*4;
       end
    case 58
       if t >= 18*lf5+0.5*te && t <= 19*lf5+0.5*te 
           fvt = sr2*4;
       elseif t >= 19*lf5+0.5*te && t <= 20*lf5+0.5*te 
           fvt = -sr2*4;
       end 
    case 59
       if t >= 20*lf5+0.5*te && t <= 21*lf5+0.5*te 
           fvt = sr2*4;
       elseif t >= 21*lf5+0.5*te && t <= 22*lf5+0.5*te 
           fvt = -sr2*4;
       end
    case 60
       if t >= 22*lf5+0.5*te && t <= 23*lf5+0.5*te 
           fvt = sr2*4;
       elseif t >= 23*lf5+0.5*te && t <= 24*lf5+0.5*te 
           fvt = -sr2*4;
       end
    case 61
       if t >= 24*lf5+0.5*te && t <= 25*lf5+0.5*te 
           fvt = sr2*4;
       elseif t >= 25*lf5+0.5*te && t <= 26*lf5+0.5*te 
           fvt = -sr2*4;
       end
    case 62
       if t >= 26*lf5+0.5*te && t <= 27*lf5+0.5*te 
           fvt = sr2*4;
       elseif t >= 27*lf5+0.5*te && t <= 28*lf5+0.5*te 
           fvt = -sr2*4;
       end 
    case 63   
       if t >= 28*lf5+0.5*te && t <= 29*lf5+0.5*te 
           fvt = sr2*4;
       elseif t >= 29*lf5+0.5*te && t <= 30*lf5+0.5*te 
           fvt = -sr2*4;
       end
    case 64
       if t >= 30*lf5+0.5*te && t <= 31*lf5+0.5*te 
           fvt = sr2*4;
       elseif t >= 31*lf5+0.5*te && t <= 32*lf5+0.5*te 
           fvt = -sr2*4;
       end
end
fvt = tif*fvt;                       % adjusting to the time intervall
        
