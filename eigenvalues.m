
A = [6 -2 5;-2 6 -5; 5 -5 13];
b=round(eigen_QR(A));
d=diag(b);
fprintf('die Eigenwerte der Matrix sind µ¹=%d, µ²=%d  µ³=%d',d(1),d(2),d(3))
function[Ak, QQ] = eigen_QR(A)
    iterations = 5;
    n = size(A,1);
    QQ = eye(n);
    Ak = A;
    for k=1:iterations
        [Q, R] = qr(Ak);
        Ak = R * Q;
        QQ = QQ * Q;
    end
end

