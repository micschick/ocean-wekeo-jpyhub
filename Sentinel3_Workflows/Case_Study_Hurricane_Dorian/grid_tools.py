#!/usr/bin/env python

import numpy as np
from scipy.ndimage.measurements import label

def surrounding_mean(xx,yy,array):
   '''
     takes mean of points around given point; fails on borders
   '''
   try:
       sum_vals = np.nanmean([array[xx-1,yy-1],\
                              array[xx,yy-1],\
                              array[xx+1,yy-1],\
                              array[xx-1,yy],\
                              array[xx+1,yy],\
                              array[xx-1,yy+1],\
                              array[xx,yy+1],\
                              array[xx+1,yy+1]])
   except:   
       sum_vals = np.nan

   return sum_vals

def gap_fill(input_array):
    '''
     fills artefact regions (assumed negative) with average of perimeter
    '''
    # set up blank array
    input_blank = input_array.copy()
    input_array[input_array<0]=np.nan
    # mask boolean array, positive values assumed legitimate
    input_blank[input_blank>=0]=0
    input_blank[input_blank<0]=1
    # set structure
    structure = np.ones((3, 3), dtype=np.int)
    # isolate 'islands'
    labelled, ncomponents = label(input_blank, structure)

    # cycle through labelled islands and get perimeter value
    for the_label in np.unique(labelled):
      sum_value = []
      indices = np.where(labelled==the_label)
      Xs = indices[0]
      if len(Xs) > np.shape(input_array)[0]*np.shape(input_array)[1]*0.5:
          continue
      Ys = indices[1]
      for ii in np.arange(0,len(Xs)):
           sum_value.append(surrounding_mean(Xs[ii],Ys[ii],input_array))
      input_array[labelled == the_label] = np.nanmean(sum_value)

    return input_array

def nine_point_spread(the_mask):
    ''' spreads the nan mask '''

    m, n = np.shape(the_mask)
    # pythonise indices
    m = m - 1
    n = n - 1

    # corners
    the_mask[0, 0] = the_mask[0, 0] + the_mask[0, 1] + \
                     the_mask[1, 0] + the_mask[1, 1]

    the_mask[0, n] = the_mask[0, n] + the_mask[0, n - 1] + \
                     the_mask[1, n] + the_mask[1, n - 1]

    the_mask[m, 0] = the_mask[m, 0] + the_mask[m, 1] + \
                     the_mask[m - 1, 0] + the_mask[m - 1, 1]

    the_mask[m, n] = the_mask[m, n] + the_mask[m, n - 1] + \
                     the_mask[m - 1, n] + the_mask[m - 1, n - 1]
    # edges
    the_mask[1:m - 1, 0] = the_mask[1:m - 1, 0] + the_mask[1:m - 1, 1] + \
                           the_mask[0:m - 2, 0] + the_mask[2:m, 0] + \
                           the_mask[0:m - 2, 1] + the_mask[2:m, 1]

    the_mask[1:m - 1, n] = the_mask[1:m - 1, n] + the_mask[1:m - 1, n - 1] + \
                           the_mask[0:m - 2, n] + the_mask[2:m, n] + \
                           the_mask[0:m - 2, n - 1] + the_mask[2:m, n - 1]

    the_mask[0, 1:n - 1] = the_mask[0, 1:n - 1] + the_mask[1, 1:n - 1] + \
                           the_mask[0, 0:n - 2] + the_mask[0, 2:n] + \
                           the_mask[1, 0:n - 2] + the_mask[1, 2:n]

    the_mask[m, 1:n - 1] = the_mask[m, 1:n - 1] + the_mask[m - 1, 1:n - 1] + \
                           the_mask[m, 0:n - 2] + the_mask[m, 2:n] + \
                           the_mask[m - 1, 0:n - 2] + the_mask[m - 1, 2:n]

    # body
    the_mask[1:m - 1, 1:n - 1] = the_mask[1:m - 1, 1:n - 1] + \
                                 the_mask[0:m - 2, 1:n - 1] + \
                                 the_mask[2:m, 1:n - 1] + \
                                 the_mask[1:m - 1, 0:n - 2] + \
                                 the_mask[0:m - 2, 0:n - 2] + \
                                 the_mask[2:m, 0:n - 2] + \
                                 the_mask[1:m - 1, 2:n] + \
                                 the_mask[0:m - 2, 2:n] + \
                                 the_mask[2:m, 2:n]
    return the_mask
