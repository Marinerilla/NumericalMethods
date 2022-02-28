clear all
syms x
format shortG
% Eingabedaten
f=6*x^5+2*x;a=1;b=2;
It=trapezint(f,a,b);% Die It-Variable erfasst den Wert des Integrals

function [Th]=trapezint(f,a,b)% die Funktion hängt nur von dem zu integrierenden Polynom
% und den Integrationsgrenzen a und b ab
    syms x
    fx=inline(f);
    % Zweite Ableitung der Funktion zur Analyse des Fehlers berechnet.
    f2=diff(diff(f)); f2x=inline(f2);
    fa=fx(a); fb=fx(b);f2a=f2x(a); f2b=f2x(b);
    % In der folgenden Zeile wurde mit Hilfe von Matlab 
    % der Wert des analytischen Integrals berechnet.
    I=double(int(fx(x),a,b));    
    h=0;
    disp('j    n            h         T(h)     │T(h)-I│/h²  (b-a)*M/12  ')
    disp('¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯')
    for j = 0:10% Um die in der Aufgabe n=2^j angegebene Bedingung zu implementieren 
        n=2^j;   
        h=(b-a)/n;
        fs=0;f2s=0;
        % Diese Schleife startet die Summierung für die Trapezregel
        for i=1:n-1
            y=a+h*i;
            fs=fs+fx(y);% Für Integrand-Funktion
            f2s=f2s+f2x(y);% Für die zweite Ableitung der Funktion
        end
        Th=0.5*h*(fa+fb+2*fs);
        M=(b-a)*(0.5*h*(f2a+f2b+2*f2s))/12; 
        e=abs(Th-I)/h^2;
        fprintf('%d\t %d      \t%0.3f\t  %0.4f\t    %0.4f\t     %0.4f\n',j,n,h,Th,e,M) 
    end
    fprintf('\nDie numerische Lösung des Integrals nach der Trapezregel ist %0.2f',Th)
end