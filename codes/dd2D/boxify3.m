function xbox = boxify3(x, Nx, Nz)
	xbox = reshape(x, 3, Nx*Nz)';
end