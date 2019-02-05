#!/usr/local/bin/julia

using PyPlot

# ==============================================================================
#	       		               Plots
# ==============================================================================

function data_plot(x,y,p1,a,c)
		#PyPlot.ioff()
		loglog(x, y, marker="o", linestyle="none")
		loglog(x, x.^(c[2]).*exp(c[1]), linestyle="--")
		xlabel("h")
		ylabel("err")
		grid("on")
		PyPlot.show()
end
