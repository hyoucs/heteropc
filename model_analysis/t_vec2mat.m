function A = t_vec2mat(v, n)

	A = zeros(n,n);
	% map to lower triangular matrix
	A(tril(true(size(A)),-1)) = v;
	% add one to diagonal and symmetrize
	A = A + A' + eye(n);

end