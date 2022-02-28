clear all
syms x1 x2
fprintf('\nNewton-verfahren für 2 gleichungen\n')
% Eingabedaten 
start1=[1;1];start2=[-1;2];start=[start1,start2];l=length(start);
g1=@(x1,x2) sin(x1)*cos(x2);
g2=@(x1,x2) x1^2+x2^2-2;
% um die Ergebnisse zu sammeln
xx1=[];xx2=[];kk=[];
% wenden Sie die Newton-Methode auf die Startvektoren an, 
% die am Anfang des Skripts angegeben sind
for i=1:1:l
    [xxk,kki]=newtoneq(g1,g2,start(:,i));
    xx1(end+1)=xxk(1);
    xx2(end+1)=xxk(2);
    kk(end+1)=kki;
end
fprintf('\nmit dem Startvektor [%d %d] das Newton-Verfahren die erforderliche Toleranz von 0,0001 in %d Iterationen erreicht hat und eine numerische Lösung von [%d %0.3f] \n',start1(1),start1(2),kk(1),xx1(1),xx2(1))
fprintf('\nmit dem Startvektor [%d %d] das Newton-Verfahren die erforderliche Toleranz von 0,0001 in %d Iterationen erreicht hat und eine numerische Lösung von [%d %0.3f] \n',start2(1),start2(2),kk(2),xx1(2),xx2(2));

function [xk,k]=newtoneq(g1,g2,x0)
    syms x1 x2
    ex1=1;ex2=1;tol=1e-4;kmax=100;k=1;
    df11=inline(diff(g1,x1));df12=inline(diff(g1,x2));
    df21=inline(diff(g2,x1));df22=inline(diff(g2,x2));
    fprintf('\n')
    disp('k       Xo1         Xk1         ex1         Xo2         Xk2         ex2 ')
    disp('¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯')
    while ((ex1>=tol)||(ex2>=tol))
        x1=x0(1);x2=x0(2);
        j11=df11(x1,x2);j12=df12(x1,x2);j21=df21(x1);j22=df22(x2);
        J=[j11,j12;j21,j22];
        F=[g1(x1,x2);g2(x1,x2)];
        xk=x0-double(J^-1*F);
        ex1=abs(xk(1)-x0(1));
        ex2=abs(xk(2)-x0(2));
        fprintf('%d\t  %0.4f\t  %0.4f\t  %0.4f\t  %0.4f\t  %0.4f\t  %0.4f\n',k,x0(1),xk(1),ex1,x0(2),xk(2),ex2)
        x0=xk;
        k=k+1;
        if k>=kmax, break, end
    end
    if xk(1)~=0
    xk(1)=abs(round(xk(1)));    
    end
end
