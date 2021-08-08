
# Seaborn splits matplotlib parameters into

    # Those that control figure style
    #   - axes_style() : returns a dictionary of parameters
    #   - set_style() : sets the matplotlib defaults. Five preset
    #                   themes: darkgrid (deflt), whitegrid, dark,
    #                   white, and ticks.
    #                   e.g. sns.set_style("whitegrid")

    # Those that scale figure elements to different contexts
    #   - plotting_context() : returns a dictionary of parameters
    #   - set_context() : sets the matplotlib defaults. Four preset
    #                     contexts, sorted by size: paper,
    #                     notebook(deflt), talk, and poster.
    #                     e.g.sns.set_context("paper")


#=======================================================================

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


# Offset sin waves, to play with
def sinplot(flip=1):
    x = np.linspace(0, 14, 100)
    for i in range(1, 7):
        plt.plot(x, np.sin(x + i * .5) * (7 - i) * flip)


#=======================================================================

# Original plot

sinplot() # <------------------------------

plt.show()
plt.cla()
plt.clf()

#=======================================================================

# Switch to seaborn defaults

sns.set_theme()

sinplot() # <------------------------------

plt.show()
plt.cla()
plt.clf()


#=======================================================================

# Play with set_syle
# BEFORE THE PLOT HAS BEEN CALLED

sns.set_style("dark")

sinplot() # <------------------------------

plt.show()
plt.cla()
plt.clf()

#=======================================================================

# (White and ticks only) Removing the top and right axes spines
# AFTER THE PLOT HAS BEEN CALLED

sns.set_style("white")

sinplot() # <------------------------------

sns.despine()

plt.show()
plt.cla()
plt.clf()

#=======================================================================

# Use the axes_style() function in a with statement
# to temporarily set plot parameters.

f = plt.figure(figsize=(6, 6))
# Width, height in inches, default 6.4, 4.8

gs = f.add_gridspec(2, 2)
# Initializes a 2x2 grid specification
# Indexing and slices to allocate multiple "cells" subplot.
# Can get funky with 3x3... but why would you do that...

with sns.axes_style("darkgrid"):
    ax = f.add_subplot(gs[0, 0])
    sinplot() # <------------------------------

with sns.axes_style("white"):
    ax = f.add_subplot(gs[0, 1])
    sinplot() # <------------------------------

with sns.axes_style("ticks"):
    ax = f.add_subplot(gs[1, 0])
    sinplot() # <------------------------------

with sns.axes_style("whitegrid"):
    ax = f.add_subplot(gs[1, 1])
    sinplot() # <------------------------------

f.tight_layout()
# Automatically adjusts subplot params so that the subplot(s) fits in to the figure area

plt.show()
plt.cla()
plt.clf()

#=======================================================================

# Legends

import matplotlib.patches as mpatches

FailingPatch = mpatches.Patch(color='red', label='Fail')
PassingPatch = mpatches.Patch(color='blue', label='Pass')
InvalidPatch = mpatches.Patch(color='black', label='Invalid')

leg = plt.legend(
    handles = [FailingPatch, PassingPatch, InvalidPatch],
    frameon = False,
    fontsize = 10,
)

for patch in leg.get_patches():
    patch.set_height(8)
    patch.set_width(8)
    patch.set_y(0)
    patch.set_x(0)

plt.show()
plt.cla()
plt.clf()
