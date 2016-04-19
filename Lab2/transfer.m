function f = transfer(x,str)
	if strcmpi(str,'lin')
		f = x;
	elseif strcmpi(str,'sigmoid')
		f = 1/(1+exp(-x));
	elseif strcmpi(str,'hyperbolic')
		f = (exp(x)-exp(-x))/(exp(x)+exp(-x))
	elseif strcmpi(str,'satlin')
		if x<0
			f=0;
		elseif and(x>0, x<=1)
			f=x;
		else
			f=1;
		end
	elseif strcmpi(str,'symsatlin')
		if x<0
			f=-1;
		elseif and(x>0,x<=1)
			f=x;
		else
			f=1;
		end
	elseif strcmpi(str,'poslin')
		if x<0
			f=0;
		else
			f=x;
		end
	elseif strcmpi(str,'bipolar')
		if x<0
			f=-1;
		else
			f=1;
		end
	end
end
	
	
	 
