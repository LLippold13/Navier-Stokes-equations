#top
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, MaxNLocator
#
def plotting(image, title='', axes=[], maxN=[0,0], ax_labels=['',''], dec=2):

	N = len(image)

	if maxN==[0,0]:
		maxN = [N,N]

	if axes==[]:
		axes = range(N)

	def format_fnx(tick_val, tick_pos):
		if int(tick_val) in range(N):
			# return np.int64(round(axes[int(tick_val)], dec))
			return round(axes[int(tick_val)], dec)
		else:
			return ''

	def format_fny(tick_val, tick_pos):
		if int(tick_val) in range(N):
			# return np.int64(round(np.flip(axes)[int(tick_val)],dec))
			return round(np.flip(axes)[int(tick_val)],dec)
		else:
			return ''



	fig, ax = plt.subplots()
	# fig = pylab.gcf()#get current figure window
	# fig.suptitle('M_xy')
	fig.canvas.set_window_title(title)

	
	im = ax.imshow(np.flipud(image), cmap='gray')
	fig.colorbar(im)
	xticks = ax.get_xticks()
	yticks = ax.get_yticks()


	ax.set_xlabel(ax_labels[0])
	ax.set_ylabel(ax_labels[1])

	ax.xaxis.set_major_formatter(FuncFormatter(format_fnx))
	ax.xaxis.set_major_locator(MaxNLocator(maxN[0]))
	ax.yaxis.set_major_formatter(FuncFormatter(format_fny))
	ax.yaxis.set_major_locator(MaxNLocator(maxN[1]))
#bottom