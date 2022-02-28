clear
syms x
% Eingabedaten
f=x^3-2*x^2+x;start=[-0.9 2];l=length(start);
d=polynomialDegree(f);rd=double(root(f,x));
fprintf('\nDie Gleichung hat %d Wurzeln deren Werte sind: ',d)
disp(transpose(rd))
fprintf('\n')
% Die folgende "Schleife für" berechnet die Näherungen der Wurzeln des 
% Polynoms mit den für das einfache und das modifizierte Newton-Verfahren 
% definierten Funktionen für die beiden angegebenen Startvektoren
for i=1:1:l
%     newton(f,start(i),rd(i));
    newtonmod(f,start(i),rd(i));
end

% In den folgenden Zeilen werden die Funktionen für 
% das einfache und das modifizierte Newton-Verfahren definiert

function newton(f,x0,r)
fprintf('\nSimple Newton-verfahren\n')
syms x
tol=1e-6;kmax=100;err1=1;a=x0;x10=x0;k1=1;
fx=inline(f);dfx=inline(diff(f));
disp("k        Xo          Xk       │Xk-Xk-1│  │X-Xk│/│X-Xk-1│²  ")
disp("¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯")
while (err1>=tol) 
    xk1=x10-fx(x10)/dfx(x10);
    err1=abs(xk1-x10);
    c1=abs(r-xk1)/abs(r-x10)^2;
    fprintf('%d\t  %f\t  %f\t   %f\t\t  %0.2f\n',k1,x10,xk1,err1,c1)
    x10=xk1;
    k1=k1+1;
    if k1>=kmax, break, end
end
fprintf('\ndie Wurzel für den Startwerte X0=%d mit hilfe auf Simple Newton-verfahren ist Xk=%f\t\n',a,xk1)
end 

function newtonmod(f,x0,r)
fprintf('\nModifizierte Newton-verfahren\n')
syms x
tol=1e-6;kmax=100;err2=10;a=x0;x20=x0;k2=1;
fx=inline(f);dfx=inline(diff(f));d2fx=inline(diff(diff(f)));
disp("k        Xo          Xk       │Xk-Xk-1│  │X-Xk│/│X-Xk-1│²  ")
disp("¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯")
while (err2>=tol) 
    xk2=x20-((fx(x20)*dfx(x20))/((dfx(x20)^2)-(fx(x20)*d2fx(x20))));
    err2=abs(xk2-x20);
    c2=abs(r-xk2)/abs(r-x20)^2;
    fprintf('%d\t  %f\t  %f\t   %f\t\t  %0.2f\n',k2,x20,xk2,err2,c2)
    x20=xk2;
    k2=k2+1;
    if k2>=kmax, break, end
end
fprintf('\ndie Wurzel für den Startwerte X0=%d mit hilfe auf Modifizierte Newton-verfahren ist Xk=%f\t\n',a,xk2)
end 
