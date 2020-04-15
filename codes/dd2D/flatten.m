function xvec = flatten(x)
	xvec = reshape(x', prod(size(x)), 1);
end