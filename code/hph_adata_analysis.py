# -*- coding: utf-8 -*-
"""
Created on Mon Jun  7 08:42:25 2021

This module contains functions for analysing an anndata object. 

It also has functions for plotting images related to anndata objects. 

@author: MerrickS
"""

import pathlib 
import h5py
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.colors

from matplotlib.colors import FuncNorm, PowerNorm

from matplotlib_scalebar.scalebar import ScaleBar
from mpl_toolkits.axes_grid1 import ImageGrid

from skimage import io, filters, morphology

import high_plex_hdf as hph

def column_percentile_normalise(array, percentile, crop=False):
    for i, col in enumerate(array.T):
        norm = np.percentile(col, percentile)
        
        if crop == True:
            col[col > norm] = norm
            array[:,i] = col/norm
        else:
            array[:,i] = col/norm
            #print(i, col.max())
    
    return array

def column_percentile_normalise_test(array, low_percentile, upper_percentile, crop=False):
    for i, col in enumerate(array.T):
        low_crop = np.percentile(col, low_percentile)
        high_crop = np.percentile(col, upper_percentile)
       
        col[col < low_crop] = low_crop
        col = col - low_crop
        col[col > high_crop] = high_crop
        array[:,i] = col/high_crop
           
    return array

def plot_cell_size_distribution(adata, bin_w=20):
    df_temp = pd.DataFrame(adata.obs.area, columns=['area'])
    df_temp['sample'] = adata.obs['sample']

    sns.displot(        df_temp, 
        x="area", 
        col="sample", 
        col_wrap=2,
        binwidth=bin_w,
        # ax=ax
        )
        
def cell_size_percentage_filter(adata, low_percentile, upper_percentile):
    lower_crop = np.percentile(a_merge.obs['area'].to_list(), low_percentile)
    upper_crop = np.percentile(a_merge.obs['area'].to_list(), upper_percentile)

    adata_crop = adata[(adata.obs['area'] > lower_crop) & (adata.obs['area'] < upper_crop)]
    
    return adata_crop

def cell_size_absolute_filter(adata, lower_crop, upper_crop):
    adata_crop = adata[(adata.obs['area'] > lower_crop) & (adata.obs['area'] < upper_crop)]
    
    return adata_crop

def mask_erode(mask):
    eroded_mask = np.zeros(shape=mask.shape).astype(int)
    
    for label in np.unique(mask[mask!=0]):
        label_mask = mask==label
        eroded_label_mask = morphology.binary_erosion(label_mask)
        eroded_mask[eroded_label_mask] = label
        
    return eroded_mask

def adata_uns_hdf_fpaths(adata, hph_keys=None):
    """
    Generates filepaths pointing to high_plex_hdf structured `hdf` files stored
    as dictionary values in `adata.uns['high_plex_hdf_paths']`. Dictionary keys
    are high_plex_hdf ids that are named in adata.obs['hph_id'], so linking 
    adata events to their source image and mask, stored in 
    adata.obs['hph_image_key'] and adata.obs['hph_mask_key'] respectively. 
    
    Hdf container filepath values are a 1D array of filepath parts from which
    filepaths are rebuilt as pathlib objects.
    
    :param adata: an anndata object with filepaths to high_plex_hdfs stored in 
        adata.uns['high_plex_hdfs'].
    :param hph_keys: list of high_plex_hdf in adata.uns['high_plex_hdfs'] for 
        which to get filepaths.
    """    
    hdf_fpaths = {}
    
    if isinstance(hph_keys, list):
        for key in hph_keys:
            hdf_path = pathlib.Path(*adata.uns['high_plex_hdf_paths'][key])
            hdf_fpaths[key] = hdf_path
    else:                
        for image in adata.uns['high_plex_hdf_paths']:
            hdf_path = pathlib.Path(*adata.uns['high_plex_hdf_paths'][image])
            hdf_fpaths[image] = hdf_path
    return hdf_fpaths
    
def adata_hph_mask_cluster_plot(
    adata, 
    cluster_name,
    mask_key=None, 
    hdf_ids=None,
    plot_names=None,
    scale_bar=True,
    scale_length=None,
    save_name=None,
    crop_dims=None,
):
    """
    A function to plot cluster annotated adata events back to the masks from 
    which the adata events were originally derived. 
    
    :param adata: an anndata object that must include links to high plex hdf 
        containers in adata.uns['high_plex_hdf_paths'] 
    :param cluster_name: the cluser annotation to plot back to masks, as named
        in an adata.obs[cluster_name] 
    :param mask_key: channels to plot from adata[var][plot_channel_var]
    :param hdf_ids: the subset of adata high_plex_hdf masks to plot
    :param save_name: if provided the plot will be saved to the given save name
    """
        
    # Make dictionary of masks used for segmentation based on the hdf image 
    # mask container filepaths held in `adata.uns['high_plex_hdf_paths']`
    hdf_fpaths = adata_uns_hdf_fpaths(adata)
        
    if mask_key == None:
        mask_key = adata.obs['hph_mask_key'].unique()[0]

    masks = {}
    for image in hdf_fpaths:
        with h5py.File(hdf_fpaths[image], 'r') as hdf:
            masks[image] = hdf[f'masks/{mask_key}'][()]
   
    # Erode masks so black pixel space in between
    eroded_masks = {}
    for sample in masks:
        eroded_masks[sample] = mask_erode(masks[sample])
        
    # Collect total cluster colors from anndata (with black for 0)
    cluster_colors = ['#000000'] + adata.uns[f'{cluster_name}_colors']
    
    # Set mask label_ids to annotated chosen cluster
    cluster_masks = eroded_masks.copy()
    for sample in cluster_masks:
        mask_df = adata.obs[[cluster_name, 'label_id']][adata.obs['sample_name_ref'] == sample]
        
        if len(np.unique(cluster_masks[sample][cluster_masks[sample] != 0])) != len(mask_df):
            print(len(np.unique(cluster_masks[sample][cluster_masks[sample] != 0])))
            print(len(mask_df))

        for i in np.unique(masks[sample][masks[sample] != 0]):
            cluster_masks[sample][cluster_masks[sample] == i] = int(
                mask_df.loc[mask_df['label_id'] == str(i), 'leiden'].iloc[0]
            ) + 1 # +1 needed to seperate cluster 0 from background 0    
            
    # Set up figure and image grid
    fig = plt.figure(figsize=(3*len(masks), 3), dpi=300)
       
    grid = ImageGrid(fig, 111,          # as in plt.subplot(111)
                     nrows_ncols=(1, len(masks)),
                     axes_pad=0.1,
                     share_all=True,
                     )

    grid_counter = 0
    for i, sample in enumerate(cluster_masks):
        
        # Change plot call order if hdf_ids present
        if isinstance(hdf_ids, list):
            sample = hdf_ids[i]           
        
        # Generate a cmap per mask that is 1 larger than max label present
        cmap = matplotlib.colors.ListedColormap(
            cluster_colors[: np.max(cluster_masks[sample]) + 1])
        boundaries = [i+1 for i in np.unique(cluster_masks[sample])] + [
            np.max(cluster_masks[sample]) + 1] 
        norm = matplotlib.colors.BoundaryNorm(boundaries, cmap.N)

        # Crop dimensions if asked
        if isinstance(crop_dims, list):
            img = cluster_masks[sample][crop_dims[0]:crop_dims[1],crop_dims[2]:crop_dims[3]]
        else:
            img = cluster_masks[sample]
           
        grid[grid_counter].imshow(img, cmap=cmap)        
        
        if isinstance(plot_names, list):
            sample = plot_names[i]
        grid[grid_counter].text(
                0.03, 0.97, 
                s = sample,
                horizontalalignment='left', 
                verticalalignment='top', 
                transform=grid[grid_counter].transAxes,
                fontsize=20,
                color='white',
                bbox=dict(facecolor='black', alpha=0.5)
                )

        grid[grid_counter].axis('off')
        
        if scale_bar == True:
            pixel = adata.obs['step_um'].unique()[0]

            if scale_length != None:
                scale_length = scale_length

            scalebar = ScaleBar(
                pixel, 
                "um", 
                fixed_value=scale_length, 
                location='lower right', 
                width_fraction=0.04,
                color='white', 
                scale_loc='top',
                box_alpha=0.5,
                box_color='black',
                font_properties={'size':14},

            )
                
        grid[grid_counter].add_artist(scalebar)
        
        grid_counter = grid_counter+1
        
    if isinstance(save_name, str):
        plt.savefig(save_name, bbox_inches='tight')
    plt.show()

def _forward(x, power=2):
    return np.power(x, (1/power), where=(x!=0))

def _inverse(x):
    return np.power(x, (power), where=(x!=0))


def adata_hph_channels_plot(
    adata, 
    plot_channel_var='channel_names',
    plot_channels=None, 
    plot_name=None,
    plot_name_split=None,
    image_names=False,
    hdf_ids=None,
    hdf_image_key='raw',
    lower_p_abs_0 = False,
    color_bar=True,
    color_bar_label = 'XRF detector counts',
    scale_bar=True,
    scale_length=None,
    save_name=None,
    legend_fontsize=None,
    power=None,
    ticks=None
):
    """
    A function to plot channels from hph_containers with filepaths stored in 
    adata['uns']['high_plex_hdf_paths'].
    
    :param adata: an anndata object that must include links to high plex hdf 
        containers in adata.uns['high_plex_hdf_paths'] 
    :param plot_channels: channels to plot from adata[var][plot_channel_var]
    :param plot_channel_var: an adata variable name that stores those 
        plot_channels to plot (in adata[var][plot_channel_var])
    :param plot_name: an adata variable name to use for output plot names. Can
        be None.
    :param hdf_ids: the subset of adata high_plex_hdfs to plot
    :param hdf_image_key: the image key in which images are held in the hdf, 
        from which images to plot will be pulled.
    """
    
    # Collect high_plex_hdf filepaths for chosen high_plex_hdf keys
    hdf_fpaths = adata_uns_hdf_fpaths(adata, hph_keys=hdf_ids)

    # Collect complete image stacks from chosen high_plex_hdf filepaths
    image_stacks = {}
    for sample in hdf_fpaths:
        hph_object = hph.hdf_to_HighPlex_Hdf(hdf_fpaths[sample])
        image_stacks[sample] = hph_object.images[hdf_image_key]
        
        if isinstance(power, int):
            image_stacks[sample] = _forward(hph_object.images[hdf_image_key], power)
        
        complete_channels_df = pd.read_hdf(
            hdf_fpaths[sample], 
            key='images/df_channel_metadata'
        )
        
    # Make dataframe of select channels from complete image stack
    select_channels = adata.var.copy()
    select_channels['channel_names'] = select_channels.index
    select_channels.index = range(len(select_channels))
    
    select_channels['complete_image_stack_index'] = complete_channels_df.index[
        complete_channels_df['shortname'].isin(
            select_channels['shortname'])].astype(int)
       
    # Set up figure and image grid
    fig = plt.figure(figsize=(3*len(image_stacks), 3*len(plot_channels)), dpi=100)

    if color_bar is True:
        cbar_mode = "edge"
    else:
        cbar_mode = None
        
    grid = ImageGrid(fig, 111,          # as in plt.subplot(111)
                     nrows_ncols=(len(plot_channels), len(image_stacks)),
                     axes_pad=0.1,
                     share_all=True,
                     cbar_mode=cbar_mode,
                     cbar_location="right",
                     cbar_size="5%",
                     cbar_pad=0.1,
                     )

    grid_counter = 0
    for row, channel in enumerate(plot_channels):  
        
        channel_imgs = dict.fromkeys(list(image_stacks.keys()))   
        img_max = dict.fromkeys(list(image_stacks.keys()))   
        img_min = dict.fromkeys(list(image_stacks.keys()))       
        lower_p=0.5
#         lower_p=lower_p
        upper_p=99.5        
        
        if isinstance(power, int):
            lower_p=0.5
            upper_p=99.5
   
        for sample in channel_imgs:
            idx = select_channels['complete_image_stack_index'][
                select_channels['channel_names'] == channel
            ].iloc[0]            
            img = image_stacks[sample][:,:,idx]
            img_min[sample] = np.percentile(img, lower_p)
            img_max[sample] = np.percentile(img, upper_p)        
            channel_imgs[sample] = img
            
        cax = grid.cbar_axes[row]
        
        # If want a 0 min need this
        if lower_p_abs_0 == True:
            img_min = {sample: 0 for sample in img_min}
#             img_min.values = 0  
        
        if isinstance(power, int) != True:
            cbar = plt.colorbar(
                plt.cm.ScalarMappable(
                        norm=plt.Normalize(
                            vmin=min(img_min.values()), 
                            vmax=max(img_max.values())
                        ), 
                        cmap='gray'
                    ), 
                cax=cax,
            )
            
        if isinstance(power, int):
            norm = PowerNorm(gamma=1/power, vmin=np.power(min(img_min.values()), power), vmax=np.power(max(img_max.values()), power))            
            cbar = plt.colorbar(
                plt.cm.ScalarMappable(
                        norm=norm, 
                        cmap='gray'
                    ), 
                cax=cax,
                ticks=ticks,
            )
            
    
        # Pad tick labels for consistent width
        tick_labels = [tick for tick in cbar.ax.get_yticks()]
        if np.max(tick_labels) < 10:
            tick_labels_padded = ['{0: <4}'.format(round(lbl,1)) for lbl in tick_labels]
        else:
            tick_labels_padded = ['{0: <4d}'.format(int(lbl)) for lbl in tick_labels]
        cbar.ax.set_yticklabels(tick_labels_padded)  # vertically oriented colorbar
        
        cbar.set_label(
            color_bar_label, 
            rotation=270,
            labelpad=17,
            fontsize = 14
        )
            
        for col, sample in enumerate(channel_imgs):
            img = channel_imgs[sample] - min(img_min.values())
            img = filters.gaussian(img, sigma=0.8)
            vmin = min(img_min.values()) - min(img_min.values())
            vmax = max(img_max.values()) - min(img_min.values())
            
            if lower_p_abs_0 == True:
                print("ITS TRUE")
                vmin = 0
            
            grid[grid_counter].imshow(
                img,
                vmin=vmin,
                vmax=vmax,
                cmap='gray'
            )
            
            grid[grid_counter].axis('off')
            
            if isinstance(legend_fontsize, int):
                pass
            else:
                legend_fontsize = 20

            if image_names is True:
                grid[grid_counter].text(
                    0.02, 0.02, 
                    s = sample,
                    horizontalalignment='left', 
                    verticalalignment='bottom', 
                    transform=grid[grid_counter].transAxes,
                    fontsize=legend_fontsize,
                    color='white',
                    bbox=dict(facecolor='black', alpha=0.5)
                    )

            if isinstance(plot_name, str):
                top_left_text = select_channels[plot_name][
                    select_channels['channel_names'] == channel
                    ].iloc[0]
                
                if isinstance(plot_name_split, str):
                    if plot_name_split in top_left_text:
                        split_text = top_left_text.split(plot_name_split)
                        top_left_text = split_text[0]
                        if len(split_text) < 2:
                            bottom_left_text = ' '
                        else:
                            bottom_left_text = split_text[1]
                        grid[grid_counter].text(
                            0.03, 0.03, 
                            s = f"({bottom_left_text}",
                            horizontalalignment='left', 
                            verticalalignment='bottom', 
                            transform=grid[grid_counter].transAxes,
                            fontsize=legend_fontsize,
                            color='white',
                            bbox=dict(facecolor='black', alpha=0.5)
                        )

                grid[grid_counter].text(
                    0.03, 0.97, 
                    s = top_left_text,
                    horizontalalignment='left', 
                    verticalalignment='top', 
                    transform=grid[grid_counter].transAxes,
                    fontsize=legend_fontsize,
                    color='white',
                    bbox=dict(facecolor='black', alpha=0.5)
                    )
                
            # Scale Bar
            if scale_bar == True:
                pixel = adata.obs['step_um'].unique()[0]

                if scale_length != None:
                    scale_length = scale_length

                scalebar = ScaleBar(
                    pixel,
                    "um", 
                    fixed_value=scale_length,
                    location='lower right',
                    width_fraction=0.04,
                    color='white',
                    scale_loc='top',
                    box_alpha=0.5,
                    box_color='black',
                    font_properties={'size':14},
                )
                grid[grid_counter].add_artist(scalebar)

            grid_counter = grid_counter+1
    plt.show()


    # Add a scale bar to each row
    if isinstance(save_name, str):
        plt.savefig(save_name, bbox_inches='tight')           
    plt.show()