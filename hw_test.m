hw07_worker = hw07();

t0 = 0; tf=3; y0=1;

% test func1
y=hw07_worker.p1(@func1, y0, [t0, tf], 10000, 'euler');
hw_assert( abs (y(end) - exact1(tf)) < 1e-2 );
y=hw07_worker.p1(@func1, y0, [t0, tf], 100, 'midpoint');
hw_assert( abs (y(end) - exact1(tf)) < 1e-2 );
y=hw07_worker.p1(@func1, y0, [t0, tf], 10, 'rk4');
hw_assert( abs (y(end) - exact1(tf)) < 1e-2 );

% test func2
y=hw07_worker.p1(@func2, y0, [t0, tf], 10000, 'euler');
hw_assert( abs (y(end) - exact2(tf)) < 1e-4 );
y=hw07_worker.p1(@func2, y0, [t0, tf], 100, 'midpoint');
hw_assert( abs (y(end) - exact2(tf)) < 1e-4 );
y=hw07_worker.p1(@func2, y0, [t0, tf], 10, 'rk4');
hw_assert( abs (y(end) - exact2(tf)) < 1e-4 );

% test func3
y=hw07_worker.p1(@func3, y0, [t0, tf], 10000, 'euler');
hw_assert( abs (y(end) - exact3(tf)) < 1e-4 );
y=hw07_worker.p1(@func3, y0, [t0, tf], 100, 'midpoint');
hw_assert( abs (y(end) - exact3(tf)) < 1e-4 );
y=hw07_worker.p1(@func3, y0, [t0, tf], 10, 'rk4');
hw_assert( abs (y(end) - exact3(tf)) < 1e-4 );

% test for p2
%> hw07.p2('euler');
%> hw07.p2('midpoint');
%> hw07.p2('rk4');

function ret = func1(t, y)
ret = y;
end

function ret = func2(t, y)
ret = -y;
end

function ret = func3(t, y)
ret = -y + t;
end

function exact = exact1(t)
exact = exp(t);
end

function exact = exact2(t)
exact = exp(-t);
end

function exact = exact3(t)
exact = 2 * exp(-t) + t - 1;
end

function hw_assert(X)
if X; fprintf('\t PASS\n'); else; fprintf('\t FAIL\n'); end
end

hw07.p3()