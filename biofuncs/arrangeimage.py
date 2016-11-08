#! /bin/env python3

import re, os, sys
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib as mpl
import numpy as np
import pandas as pd
from PIL import Image
from functools import reduce
from multiprocessing import Pool

mpl.use('Agg')                  # When no X service run

# ------------------------------------------------------------------------------

def arrange_image(image_df):
    # {{{

    '''arrange_image is used for arranging cell photos into one big image by using matplotlib. Input a pandas DataFrame contains the figures' path arranged the way you wanted. Then you will get the plot object. Finally, save the fig by matplotlib method.'''
    # arrange images according to filename DataFrame
    fig, ax = plt.subplots(image_df.shape[0], image_df.shape[1])
    for i in range(image_df.shape[0]):
        for j in range(image_df.shape[1]):
            ax[i, j].set_frame_on(False)
            ax[i, j].tick_params(
                axis        = 'both',
                which       = 'both',
                bottom      = 'off',
                labelbottom = 'off',
                top         = 'off',
                labeltop    = 'off',
                left        = 'off',
                labelleft   = 'off',
                right       = 'off',
                labelright  = 'off'
            )
            if not os.path.exists(image_df.ix[i, j]):
                continue
            img = Image.open(image_df.ix[i,j])
            ax[i, j].imshow(img)
            if j == 0:
                ax[i, j].set_ylabel(image_df.index[i])
            if i == 0:
                ax[i, j].set_title(image_df.columns[j])
    fig.subplots_adjust(
        left=None,
        bottom=None,
        right=None,
        top=None,
        wspace=0,
        hspace=0.01
    )
    #fig.tight_layout()
    return fig
    # }}}
 
# ------------------------------------------------------------------------------
if __name__ == '__main__':
    base = '~/YD'
    ctrl = ['DMSO', 'CTRL']
    drug = ['6TG', 'BORTEZONIB', 'BI2536', 'ALISERTIB']
    cell = ['K562', 'GM12878']
    times = ['24H', '38H', '48H', '61H', '72H',
             '4days', '5days', '6days', '7days']
    dose = ['1.1ng', '11ng', '110ng', '1100ng', '11000ng']


    ctrl_df = pd.DataFrame(
        np.array(
            [[
                os.path.join(base, x.split('_')[1], cell[0], x)
                for x in reduce(
                        lambda x, y: x+y,
                        map(
                            os.listdir,
                            [os.path.join(base, actrl, cell[0])
                             for actrl in ctrl]
                        )
                ) if x.find(timepoint) > 0
            ] for timepoint in times]
        ).T,
        columns = [x.replace('days', 'D') for x in times]
    )
                

    exp_df = [
        pd.DataFrame(
            np.array(
                [[os.path.join(base, drug[0], cellone,
                               timepoint.lower(), x)
                  for x in [
                          '_'.join(
                              [cellone[0],
                               drug[0][:3].lower(),
                               dos,
                               timepoint + '.tif']
                          ) for dos in dose]
                ] for timepoint in times]
            ).T,
            columns = [x.replace('days', 'D') for x in times]
        ) for cellone in cell
    ] + [
        pd.DataFrame(
            np.array(
                [[
                    os.path.join(base, drugone, cellone,
                                 timepoint.lower(), x)
                    for x in [
                            '_'.join(
                                [cellone[0],
                                 drugone[:3].title(),
                                 dos,
                                 timepoint + '.tif']
                            ) for dos in dose]
                ] for timepoint in times]
            ).T,
            columns = [x.replace('days', 'D') for x in times]
        ) for cellone in cell for drugone in drug[1:]
    ]

    img_df = [pd.concat([exp, ctrl_df], ignore_index = True)
              for exp in exp_df]

    # for i in range(len(img_df)):
    #     img_df[i].index = dose + ctrl
    #     filename = '_'.join(
    #         os.path.basename(img_df[4].ix[0,0]).split('_')[0:2]
    #     ) + '.pdf'
    #     myimg = arrange_image(img_df[i])
    #     myimg.savefig(filename, format = 'pdf', dpi = 1800)
    #     del myimg
    def tmpplotfunc(i):                                     
        img_df[i].index = dose + ctrl                       
        filename = '_'.join(                                
            os.path.basename(img_df[i].ix[0,0]).split('_')[0:2]
        ) + 'p.pdf'                                         
        myimg = arrange_image(img_df[i])                    
        myimg.savefig(filename, format = 'pdf', dpi = 1800) 
        del myimg                                           
                                                            
    with Pool(10) as p:                                     
        p.map(tmpplotfunc, range(len(img_df)))


