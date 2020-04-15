function xbox = boxify3(x)
	xlen = length(x)/3;
	xbox = [x(1:xlen) x(xlen+1:2*xlen) x(2*xlen+1:end)];
end