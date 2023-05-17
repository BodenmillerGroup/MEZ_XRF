# -*- coding: utf-8 -*-
"""
Created on Mon Jun  7 08:42:25 2021

This module contains classes and functions for repacking and normalising 2D 
deconvoluted XRF plots. The main input is the `.h5` files output by PyMCA 
XRF deconvolution. 

@author: MerrickS
"""

import h5py
import numpy as np
import pandas as pd

import high_plex_hdf

def unpack_pymca_h5(hdf_fpath):
    """
    Unpacks XRF deconvoluted channels from the PyMCA hdf output. 
    
    Image channels [c] are imported from seperately named 2D array datasets 
    [x, y] and stack to a 3D numpy array [x, y, c].
    
    A dataframe of plots names that is indexed matched to the generated image
    stack is also generated. 
    
    The fpico mask for that datset if also collected. 
    
    :param hdf_fpath: path to hdf from which the sample name is extracted
    :type hdf_fpath: string
    :param node: hdf node holding plots of interest
    :type node: string
    
    :return image_stack: a stacked numpy array of channel images
    :return df_plots: a dataframe of plot names held in 'plot_channel' where 
        row index equates to plot_images channel index
    :return fpico_mask: a mask of beam intensity recordings
    """
    with h5py.File(hdf_fpath, 'r') as hdf:           
        # turn plots to 3D array
        plots_node = f'{hdf_fpath.stem}/plotselect'
        if plots_node in hdf:
            channels = list(hdf[plots_node].keys())
            df_plots = pd.DataFrame(channels, columns=['plot_channel'])

            for i, plot in enumerate(hdf[plots_node].keys()):
                if i == 0:
                    image_stack = hdf[f'{plots_node}/{plot}']
                    image_stack = np.expand_dims(image_stack, -1)
                else:
                    image_add = hdf[f'{plots_node}/{plot}']
                    image_add = np.expand_dims(image_add, -1)
                    image_stack = np.concatenate([image_stack, image_add], -1) 
        assert image_stack.shape[-1] == len(channels)
        
        # add fpico mask to masks
        if 'fpico_mask' in hdf:
            fpico_mask = hdf['fpico_mask'][()]
        else:
            fpico_mask = np.zeros(image_stack.shape[:1]) # empty fpico ask if none       
        
    return image_stack, df_plots, fpico_mask

class XrfImageMaskHDF:
    """
    This class stores image stacks (3D numpy arrays [x, y, c], masks and 
    associated metadata. It has functions for normalising 2D plots held as 
    [c]hannels in image stacks. Normalised image stacks are stored under a 
    name related to the normalisation method in an `images` dictionary, 
    alongside raw image stacks. 
    
    Multiple image stacks, as well as related masks and metadata can be 
    exported from the `DeconvolutedXrfHdfRepack` class to an hdf file with a 
    comparable structure using the classes `export_hdf` function. 
    
    The exported hdf is useful for keeping related image and mask datasets 
    together for archival or sharing. The exported hdf structure also 
    streamlines downstream analysis of images using the related mask (e.g. mask
    label measurement of a chosen image stack).
    """
    def __init__(self, images, channel_metadata, masks, sample_metadata):
        self.images = images
        self.channel_metadata = channel_metadata
        self.sample_metadata = sample_metadata
        self.masks = masks
        
    def scatter_normalise(self, scatter_set_max):
        scatter_factor = scatter_set_max[self.sample_metadata['scatter_set']]
        #idx = self.channel_metadata[self.channel_metadata['plot_channel'] == 'Scatter_Compton000'].index[0]
        idx = self.channel_metadata[self.channel_metadata['plot_channel'] == 'Scatter_Peak000'].index[0]
        scatter_image = self.images['raw'][:, :, idx]
        scatter_factor = np.mean(scatter_image) / scatter_factor
        
        image_norm = self.images['raw'].copy()
        for i in range(image_norm.shape[-1]):
            image_norm[:, :, i] = image_norm[:, :, i]/scatter_factor
            
        self.images['1_compton_norm'] = image_norm
        
    def beam_intensity_normalise(self):
        image_norm = self.images['1_compton_norm'].copy()
               
        for i in range(image_norm.shape[-1]):
            image_norm[:, :, i] = np.divide(
                image_norm[:, :, i], self.masks['fpico_mask'], where=self.masks['fpico_mask']!=0)
            
        self.images['2_beam_intensity_norm'] = image_norm
              
    
    def export_hdf(self, output_fpath):
        with h5py.File(output_fpath, 'a') as hdf:
            # images exported as 3D stack [x, y, c]
            images_node = 'images'
            hdf.require_group(images_node)
            for image in self.images:
                image_sub_node = f'{images_node}/{image}'
                if image_sub_node in hdf:
                    hdf[image_sub_node][...] = self.images[image]
                if image_sub_node not in hdf:
                    hdf.create_dataset(image_sub_node, data=self.images[image], compression='gzip')

            # masks exported as 2D array [x, y]
            masks_node = 'masks'                
            hdf.require_group(masks_node)    
            if self.masks == None:
                pass
            if isinstance(self.masks, dict):
                for mask in self.masks:
                    mask_sub_node = f'{masks_node}/{mask}'
                    if mask_sub_node in hdf:
                        hdf[mask_sub_node][...] = self.masks[mask]
                    if mask_sub_node not in hdf:
                        hdf.create_dataset(mask_sub_node, data=self.masks[mask], compression='gzip')

            # sample_metadata exported as attributes
            if self.sample_metadata == None:
                pass
            else:
                for sample_data in self.sample_metadata:
                    if sample_data in hdf:
                        hdf[sample_data][...] = self.sample_metadata[sample_data]                
                    if sample_data not in hdf:
                        hdf[sample_data] = self.sample_metadata[sample_data]                
                                
        # dataframe of channel information exported with pandas standard hdf export function
        if f'{images_node}/df_channel_metadata' not in h5py.File(output_fpath, 'r'):
            self.channel_metadata.to_hdf(output_fpath, f'{images_node}/df_channel_metadata', mode='a')