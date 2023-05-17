# -*- coding: utf-8 -*-
"""
Created on Mon Jun  7 08:42:25 2021

This file contains helper functions for XRF analysis and handling of hdf files

@author: MerrickS
"""
import h5py
import pandas as pd
import numpy as np

def aggregate_images(image_stack, img_idx, agg_function='mean'):
    """
    From [x, y, c] summarise a list of [x, y] images specified by a list of 
    channel indexes. Summary can be either sum or mean. 
    
    :param image_stack: a 3D image stack [x, y, c]
    :param img_idx: a list of image stack indexes to aggregate
    :param img_idx: a list of image stack indexes to aggregate

    """
    for num, idx in enumerate(img_idx):
        if num == 0:
            aggregate_img = image_stack[:, :, idx]
            aggregate_img = np.expand_dims(aggregate_img, -1)
        else:
            nxt_img = image_stack[:, :, idx]
            nxt_img = np.expand_dims(nxt_img, -1)
            aggregate_img = np.concatenate([aggregate_img, nxt_img], -1) 

    if agg_function=='mean':
        aggregate_img = np.mean(aggregate_img, axis=-1)
    if agg_function=='sum':
        aggregate_img = np.sum(aggregate_img, axis=-1)
    
    return aggregate_img

def aggregate_images_norm(image_stack, img_idx, agg_function='mean'):
    """
    From [x, y, c] summarise a list of [x, y] images specified by a list of 
    channel indexes. Summary can be either sum or mean. 
    
    :param image_stack: a 3D image stack [x, y, c]
    :param img_idx: a list of image stack indexes to aggregate
    :param img_idx: a list of image stack indexes to aggregate

    """
    for num, idx in enumerate(img_idx):
        if num == 0:
            aggregate_img = image_stack[:, :, idx]

            p99 = np.percentile(aggregate_img, 99)
            aggregate_img[aggregate_img > p99] = p99
            aggregate_img = (aggregate_img / p99) * 255
            
            aggregate_img = np.expand_dims(aggregate_img, -1)
        else:
            nxt_img = image_stack[:, :, idx]
            
            p99 = np.percentile(nxt_img, 99)
            nxt_img[nxt_img > p99] = p99
            nxt_img = (nxt_img / p99) * 255
            
            nxt_img = np.expand_dims(nxt_img, -1)
            aggregate_img = np.concatenate([aggregate_img, nxt_img], -1) 

    if agg_function=='mean':
        aggregate_img = np.mean(aggregate_img, axis=-1)
    if agg_function=='sum':
        aggregate_img = np.sum(aggregate_img, axis=-1)
    
    return aggregate_img

def get_hdf_metadata(hdf_fpath, df_hdf_metadata, panel_dir):
    """
    Gets sample metadata for the provided hdf_fpath from the df_hdf_metadata. 
    Requires that the hdf name is represented in either a 'scanset' or 
    'hdf_filename' column in df_hdf_metadata. 

    :param hdf_fpath: path to hdf from which the sample name is extracted
    :type hdf_fpath: string
    :param df_hdf_metadata: path to hdf from which the sample name is extracted
    :type df_hdf_metadata: pandas.Dataframe
    
    :return sample_metadata: metadata related to sample properties
    :type sample_metadata: pandas.Series
    :return df_panel: channel specfic information, e.g. antibody, target, isotope 
    :type df_panel: pandas.Dataframe
    """
    if 'stitch' in str(hdf_fpath):
        print(hdf_fpath.name, 'is a stitch file')
        sample_metadata = df_hdf_metadata[df_hdf_metadata['scanset'] == hdf_fpath.stem].iloc[0]
    else:
        print(hdf_fpath.name, 'is a non stitch file')
        sample_metadata = df_hdf_metadata[df_hdf_metadata['hdf_filename'] == hdf_fpath.stem].iloc[0]
    
    sample_metadata = sample_metadata[1:].to_dict() # [1:] removes index

    panel = sample_metadata['AB_panel']   
    
    if pd.isna(panel):
        df_panel = False
        print(f'No panel for {hdf_fpath.name}')

    else:
        print(panel, f'is matched to sample {hdf_fpath.name}')
        df_panel = pd.read_csv(panel_dir / f'{panel}.csv')  
    
    return sample_metadata, df_panel

def get_hdf_full_plot_df(df_plots, df_panel):
    """
    Incorporates panel information as new columns to the df_plots indexed 
    dataframe produced by `hdf_plots_import`. The expanded dataframe generated
    has columns useful for referencing plots from an acquisition. For instance
    element, isotope, antibody clone or shortnames for plot channels. Columns 
    in the expanded dataframe can also store information indicating key 
    channels. For instance, if channels should be used for generating figures, 
    or if channels are to be used for later nuclear/cytoplasmic segmentation.
    
    :param df_plots: an indexed dataframe of channel names
    :param df_panel: a dataframe of information related to channel names
    
    :return df_plots_panel: a complete dataframe of key channel information
        with indexes related to an image stack [x, y, index]. 

    """
    if df_panel is False:
        df_plots_panel = df_plots # just use the df_plots when no panel
        
    else:
        df_panel = df_panel.drop_duplicates(subset='xrf_emission')    
        df_plots_panel = pd.merge(df_plots, df_panel, left_on='plot_channel', right_on='xrf_emission', how='left')

        if df_plots.shape[0] == df_plots_panel.shape[0]:
            pass
        else:
            print('df_plots and df_plots panel different length no longer related to image')
        
    df_plots_panel = force_object_columns_to_string(df_plots_panel)

    return df_plots_panel


def force_object_columns_to_string(df):    
    if isinstance(df, pd.Series):
        df = df.astype(str)
    if isinstance(df, pd.DataFrame):
        # Ensure string type in object columns
        string_columns = df.columns[df.dtypes == object]
        df.loc[:, string_columns] = df[string_columns].applymap(str)
    return df

def hdf_scalar(hdf, node):
    """
    This function returns a scalar dataset from an hdf at a specified node
    """
    with h5py.File(scan, 'r') as hdf:
        dset = hdf[node][()]
    return dset

def hdf_array(hdf, node):
    """
    This function returns an array dataset from an hdf at a specified node
    """
    with h5py.File(scan, 'r') as hdf:
        dset = hdf[node][:]
    return dset

def hdf_dataset(hdf, node):
    """
    This function checks the dataset shape at a specified node, then extract that 
    dataset from an hdf at the specified node with the appropriate function. 
    
    :param hdf: hdf filepath
    :type hdf: str
    :param node : hdf path pointing to the dataset in hdf
    :type node: str
    
    :return dset: scalar or array type dataset held in hdf node
    """
    with h5py.File(scan, 'r') as hdf:
        try:
            hdf[node]
            data_shape = hdf[node].shape
            
            if len(data_shape) == 0:
                dset = hdf_scalar(hdf, node)
            else:
                dset = hdf_array(hdf, node)
                
        except KeyError:
            dset = []
            print(f'Could not find hdf dataset for {hdf}')
            
    return dset

def rand_cmap(nlabels, type='bright', first_color_black=True, last_color_black=False, verbose=True, rseed=False):
    """
    The random color map function (rand_cmap) is useful for making a colourmap 
    to distinguish the labels of an image segmentation. Code provided by delestro
    https://github.com/delestro/rand_cmap   
    
    Creates a random colormap to be used together with matplotlib. Useful for segmentation tasks
    :param nlabels: Number of labels (size of colormap)
    :param type: 'bright' for strong colors, 'soft' for pastel colors
    :param first_color_black: Option to use first color as black, True or False
    :param last_color_black: Option to use last color as black, True or False
    :param verbose: Prints the number of labels and shows the colormap. True or False
    :param rseed: Sets a seed value to maintain colours the number of labels and shows the colormap. True or False

    :return: colormap for matplotlib
    """
    from matplotlib.colors import LinearSegmentedColormap
    import colorsys
    import numpy as np

    if isinstance(rseed, int):
        np.random.seed(seed=rseed)

    if type not in ('bright', 'soft'):
        print ('Please choose "bright" or "soft" for type')
        return

    if verbose:
        print('Number of labels: ' + str(nlabels))

    # Generate color map for bright colors, based on hsv
    if type == 'bright':
        randHSVcolors = [(np.random.uniform(low=0.0, high=1),
                          np.random.uniform(low=0.2, high=1),
                          np.random.uniform(low=0.9, high=1)) for i in range(nlabels)]

        # Convert HSV list to RGB
        randRGBcolors = []
        for HSVcolor in randHSVcolors:
            randRGBcolors.append(colorsys.hsv_to_rgb(HSVcolor[0], HSVcolor[1], HSVcolor[2]))

        if first_color_black:
            randRGBcolors[0] = [0, 0, 0]

        if last_color_black:
            randRGBcolors[-1] = [0, 0, 0]

        random_colormap = LinearSegmentedColormap.from_list('new_map', randRGBcolors, N=nlabels)

    # Generate soft pastel colors, by limiting the RGB spectrum
    if type == 'soft':
        low = 0.6
        high = 0.95
        randRGBcolors = [(np.random.uniform(low=low, high=high),
                          np.random.uniform(low=low, high=high),
                          np.random.uniform(low=low, high=high)) for i in range(nlabels)]

        if first_color_black:
            randRGBcolors[0] = [0, 0, 0]

        if last_color_black:
            randRGBcolors[-1] = [0, 0, 0]
        random_colormap = LinearSegmentedColormap.from_list('new_map', randRGBcolors, N=nlabels)

    # Display colorbar
    if verbose:
        from matplotlib import colors, colorbar
        from matplotlib import pyplot as plt
        fig, ax = plt.subplots(1, 1, figsize=(15, 0.5))

        bounds = np.linspace(0, nlabels, nlabels + 1)
        norm = colors.BoundaryNorm(bounds, nlabels)

        cb = colorbar.ColorbarBase(ax, cmap=random_colormap, norm=norm, spacing='proportional', ticks=None,
                                   boundaries=bounds, format='%1i', orientation=u'horizontal')

    return random_colormap