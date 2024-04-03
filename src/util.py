# util.py
import os, io, math, scipy, numpy as np, pandas as pd, matplotlib.pyplot as plt, seaborn as sns
from PIL import Image

def prompt_user_input():

    ''' prompt user inputs for analysis and directly store inputs as global variables '''
    global input_dir, keyword, output_dir, experiment_date, mouse_ID, slice_number, drug_frame, frame_rate, output_path

    # input input and output directories
    input_dir = input('Enter the path to the folder containing the input files: ')
    keyword = input('Enter common keyword for the input files, hit enter if not applicable: ')
    output_dir = input('Enter the path to the folder containing the output files: ')

    # input metadata related to the experiment
    experiment_date = input('Enter the experiment date (e.g. 2023-11-06): ')
    mouse_ID = input('Enter the mouse ID: ')
    slice_number = int(input('Enter the slice number: '))
    drug_frame = float(input('Enter the frame number of drug application (enter NA if no drug was applied): '))
    frame_rate = int(input("Enter the frame rate (in Hz) of the recording (e.g. 20): "))
    output_filename = input('Enter the output file name: ') + '.xlsx'
    output_path = os.path.join(output_dir, output_filename)
 
def find_files(input_dir, keyword):
    '''takes an input directory and find the trace csv file, ROA mask tif, cell mask tif, and metadata path'''

    input_files = os.listdir(input_dir)
    input_files = [f.lower() for f in input_files]

    exp_index = [keyword.lower() in f for f in input_files]
    exp_files = np.take(input_files, np.where(exp_index)[0])
    
    csv_index = ['.csv' in f for f in exp_files]
    ROA_index = ['roa' in f for f in exp_files]
    cell_index = ['cell' in f for f in exp_files]
    
    csv_path = os.path.join(input_dir, exp_files[csv_index.index(True)])
    ROA_mask_path = os.path.join(input_dir, exp_files[ROA_index.index(True)])
    cell_mask_path = os.path.join(input_dir, exp_files[cell_index.index(True)])
    
    print("Found the following files: \n")
    print("CSV file: ", csv_path)
    print("ROA mask: ", ROA_mask_path)
    print("Cell mask: ", cell_mask_path)

    return csv_path, ROA_mask_path, cell_mask_path

def read_masks(ROA_mask_path, cell_mask_path):
    ROA_map_labeled, ROA_map_count = read_tif(ROA_mask_path)
    cell_map_labeled, cell_count = read_tif(cell_mask_path)

    print("ROA mask contains", str(ROA_map_count), "ROAs")
    print("Cell mask contains", str(cell_count), "cells")

    return ROA_map_labeled, ROA_map_count, cell_map_labeled, cell_count

def raw_to_filtered(csv_path, order = 4, cutoff = 0.4):
    '''
    read in and convert raw traces to filtered traces
    '''
    
    global ROA_count, frame_count
    # read in the csv data file in as a dataframe and transpose
    # the csv file has no header and first column is the number of ROA
    # in the original data frame, each row represents a frame and each column represents a ROA
    print("Reading in file: ", csv_path, "\n\n")
    raw_data = pd.read_csv(csv_path, header = None, index_col = 0)
    # transpose data so that each row represents a ROA and each column represents a frame
    # and transform to a numpy 2D array containing # of ROA lists,
    # each list containing signal in each frame for this ROA
    raw_traces = raw_data.transpose().to_numpy() 
    ROA_count, frame_count = check_traces(raw_traces)
    print()

    # apply a lowpass Butterworth filter with a 4th order filter at the cutoff of 0.4 Hz
    print("Applying a lowpass Butterworth filter with a", str(order), "th order filter at the cutoff of", str(cutoff), "Hz")
    b, a = scipy.signal.butter(order, cutoff, 'low', analog=False)
    filtered_traces = scipy.signal.filtfilt(b, a, raw_traces)

    return filtered_traces

def check_traces(traces):
    ''' 
    check the raw traces and print out the number of ROA and the number of frames 
    '''
    
    (ROA_count, frame_count) = traces.shape

    print("The current file contains: ")
    print("Number of ROA: ", ROA_count)
    print("Number of frames: ", frame_count)
    return ROA_count, frame_count

def correct_shift(filtered_traces):
    
    '''
    'correct' the shift of baseline level in traces using linear regression for subsequent baseline determination
    specifically, corrected traces are generated by subtracting 0.5 * slope at the current x location from the input traces
    
    note that this correction is not necessary if the baseline is stable
    and that this correction should not be applied to recordings when there is obvious z shift

    prints out the distribution of slopes
    return the corrected traces and dataframe containing slope and intercept of regression
    '''

    (ROA_count, frame_count) = filtered_traces.shape

    slope = []
    intercept = []
    x = np.arange(frame_count)

    for i_ROA in range(ROA_count):
        lm = scipy.stats.linregress(x, filtered_traces[i_ROA])
        slope.append(lm.slope)
        intercept.append(lm.intercept)
    reg = pd.DataFrame({'slope': slope, 'intercept': intercept})
    slope = np.array(slope)
    corrected_traces = filtered_traces - 0.5 * x * slope[:,None]

    # visualize distribution of lm slope
    sns.histplot(reg.slope).set(title = 'Slope distribution histogram')
    return corrected_traces, reg

def check_correction(filtered_traces, corrected_traces, reg):

    '''visual checkpoint for correction'''

    vis_threshold = float(input("Enter the slope cutoff (absolute value) for visualizing correction:"))

    (ROA_count,frame_count) = filtered_traces.shape
    x = np.arange(frame_count)

    for i_ROA in range(ROA_count):
        if abs(reg.slope[i_ROA]) > vis_threshold:
            plt.figure()
            plt.plot(x,filtered_traces[i_ROA],'b-', label = 'original trace')
            plt.plot(x, x * reg.slope[i_ROA] + reg.intercept[i_ROA], 'r-', label = 'regression line')
            plt.plot(x, corrected_traces[i_ROA], 'g-', label = 'corrected trace')
            plt.legend(loc = "upper left")
            plt.title('ROA ' + str(i_ROA + 1)) # ROA ID starts from 1
            plt.show()

def mean_abs_dist(array_2d):
    ''' take a two dimensional numpy array as input
    calculate row-wise average absolute distance from the row mean
    return a one dimensional numpy array '''
    return np.mean(np.absolute(array_2d - np.mean(array_2d, axis=1)[:,None]), axis=1)

def find_roots(trace, threshold):
    # this function does not account for when one point of signal is AT the threshold
    # which I guess rarely happens anyway

    '''input trace (1D array) and threshold (one number)
    look for where the signal crosses threshold
    return the x axis intercepts where the signal crosses threshold '''
    trace_adjusted = trace - threshold
    cross_bool = np.abs(np.diff(np.sign(trace_adjusted))).astype(bool)
    
    y1 = trace_adjusted[:-1][cross_bool]
    y2 = trace_adjusted[1:][cross_bool]
    x1 = np.where(cross_bool)[0]
    x0 = x1 - y1/(y2-y1) # solve for the x axis intercept (AKA when y=0)
    return x0

def find2points(number, array):
    '''find the two most close points'''
    for index in range(len(array)):
        if number - array[index] < 0:
            lpoint = array[index-1]
            rpoint = array[index]
            break
    return (lpoint, rpoint)

def find_signal_boundary(trace, signal_threshold, baseline_threshold):
    ''' 
    takes a 1D trace, signal_threshold, and baseline_threshold (to call beginning and end of an event)
    return tuples of event start and end points (frames) 
    '''

    signal_intercept = find_roots(trace, signal_threshold)
    baseline_intercept = find_roots(trace, baseline_threshold)

    signal_list = []
    for i in signal_intercept:
        if i > baseline_intercept[0] and i < baseline_intercept[-1]:
            m = find2points(i,baseline_intercept)
            if m not in signal_list:
                signal_list.append(m)
    signal_boundary = [(math.floor(lpoint),math.ceil(rpoint)) for (lpoint,rpoint) in signal_list]
    return signal_boundary, signal_intercept, baseline_intercept

def calc_dff(traces, signal_frames = None, baseline_start = 0, baseline_end = -1):

    ''' 
    generates df/f traces based on signal frames 
    designed to work with filtered/smoothed traces but can also work on raw traces
    returns two two-dimensional arrays: dff_traces and dff_traces_nosignal

    Args:
    traces: a two-dimensional array of traces, each row represents a ROA and each column represents a frame
    signal_frames: a two-dimensional boolean array corresponding to each ROA and each frame, true if frame is considered as signal
    baseline_start: optional, the starting frame for baseline calculation, default 0
    baseline_end: optional, the ending frame for baseline calculation, default -1 (end of the trace)

    Returns:
    dff_traces: a two-dimensional array of delta F/F trace based on the provided thresholds
    dff_traces_nosignal: a two-dimensional array of delta F/F trace based on the provided thresholds but free of signal frames (signal frames are set as NaN)
    '''

    # generate a two-dimensional array from traces but free of signal frames
    if signal_frames is None: 
        signal_frames = np.zeros(traces.shape) # if signal_frames are not provided, assume no frame is signal
    traces_nosignal = np.multiply(traces, signal_frames == False) # set signal frames to 0
    traces_nosignal[traces_nosignal == 0] = np.nan # change 0 to nan to remove signal frames in calculation

    baselines = np.nanmean(traces_nosignal[:,baseline_start:baseline_end], axis = 1) # calculate the baseline averages from signal-removed traces (ignoring nan)
    dff_traces = (traces - baselines[:,None])/baselines[:,None] 

    dff_traces_nosignal = np.multiply(dff_traces, signal_frames == False) # set signal frames to 0
    dff_traces_nosignal[dff_traces_nosignal == 0] = np.nan # change 0 to nan to remove signal frames in calculation

    return dff_traces, dff_traces_nosignal

def find_signal_frames(filtered_traces, signal_frames = None, signal_threshold = 3, baseline_start = 0, baseline_end = -1):
   
    '''
    finds signals in the filtered traces using signal_threshold and onset_threshold
    signal is defined as any event that is larger than signal_threshold * SD of the dF/F trace (default 2SD)

    Args:
        filtered_traces: a two-dimensional array of filtered traces
        signal_frames: a two-dimensional boolean array corresponding to each ROA and each frame, true if frame is considered as signal
        signal_threshold (int or float): optional, signal_threshold (3 as default) * dF/F baseline SD as signal threshold
        baseline_start: optional, the starting frame for baseline calculation, default 0 (beginning of the recording)
        baseline_end: optional, the ending frame for baseline calculation, default -1 (end of the recording)


    Returns:
        dff_traces: a two-dimensional array of delta F/F trace based on the provided thresholds
        signal_frames: a two-dimensional boolean array corresponding to each ROA and each frame, true if frame is considered as signal 
        signal_boundaries: each ROA has a list of tuples of event start and end points
    '''
    
    dff_traces, dff_traces_nosignal = calc_dff(filtered_traces, signal_frames, baseline_start, baseline_end) 
    baselines = np.nanmean(dff_traces_nosignal, axis = 1) # calculate the baseline averages from signal-removed dff traces (ignoring nan)
    thresholds = signal_threshold * np.nanstd(dff_traces_nosignal, axis = 1) # calculate the baseline standard deviations from signal-removed dff traces (ignoring nan)
    
    # initialize a new array to store signal frames
    ROA_count, frame_count = filtered_traces.shape
    new_signal_frames = np.zeros((ROA_count, frame_count))
    signal_boundaries = []

    # iterate through each ROA to find signals
    for i_ROA in range(0, ROA_count):
        signal_boundary, _ , _ = find_signal_boundary(dff_traces[i_ROA,], thresholds[i_ROA], baselines[i_ROA])
        signal_boundaries.append(signal_boundary)

        for j_frame in range(0, frame_count):
            for (l,r) in signal_boundary:
                if j_frame >= l and j_frame <= r:
                    new_signal_frames[i_ROA, j_frame] = 1
                    
    return dff_traces, baselines, thresholds, new_signal_frames.astype(bool), signal_boundaries

def iterative_baseline(filtered_traces, signal_frames = None):

    '''
    determine the bassline iteratively
    '''
    signal_frames_previous = signal_frames
    n_iteration = int(input("Enter the number of iterations for signal detection (e.g. 3): "))
    print()

    for iter in range(n_iteration):
        print(f"Processing round {iter+1} of signal detection...")
        if iter < n_iteration - 1: # only store signal_frames for the initial iterations
            _ , _, _, signal_frames, _ = find_signal_frames(filtered_traces, signal_frames_previous)
        else:
            dff_traces, baselines, thresholds, signal_frames, signal_boundaries = find_signal_frames(filtered_traces, signal_frames_previous)
        
        check_ROA(signal_frames) # check current signal detection results
        signal_frames_previous = signal_frames
        print()
    
    return dff_traces, baselines, thresholds, signal_frames, signal_boundaries

def check_ROA(signal_frames):
    '''check and report how many ROAs have signals'''
    ROAs_with_signal = np.sum(signal_frames, axis = 1) > 0

    print("ROAs with signal: ", np.sum(ROAs_with_signal))
    print("ROAs without signal: ", np.sum(ROAs_with_signal == False))

def analyze_signal(dff_traces, signal_frames, signal_boundaries, frame_rate, drug_frame = 'NA'):

    '''analyze the signals and return a dataframe with signal stats (columns) for each individual signal/event (rows)'''
    

    ROA_count = dff_traces.shape[0] # number of ROAs
    noise = np.multiply(dff_traces, signal_frames == False).max(axis = 1) # extract baseline noise (max amplitude within baseline) for each ROA

    # initialize lists to store signal stats
    ROA_ID = []
    start_index = []
    end_index = []
    AUC = []
    amplitude = []
    signal_to_noise = []
    peak_index = []
    peak_time = []
    rise_time = []
    decay_time = []
    half_width = []
    duration = []
    inter_event_interval = []

    # iterate through each ROA to extract signal stats
    for i_ROA in range(0, ROA_count):

        for j_signal in range(0,len(signal_boundaries[i_ROA])):

            (lpoint, rpoint) = signal_boundaries[i_ROA][j_signal]

            ROA_ID.append(i_ROA + 1)
            start_index.append(lpoint)
            end_index.append(rpoint)

            event_trace = dff_traces[i_ROA,lpoint:rpoint+1] # subset out only the signal/event

            AUC.append(scipy.integrate.simpson(event_trace, dx = 1/frame_rate)) # area under the curve using Simpson's rule
            amplitude.append(max(event_trace))  # signal amplitude
            signal_to_noise.append(max(event_trace)/noise[i_ROA]) # signal to noise ratio

            max_index = np.array([np.argmax(event_trace)]) # max dff index/frame number within the signal/event range
            peak_index.append(lpoint + max_index[0]) # max dff index (frame number) within the whole trace 
            peak_time.append(peak_index[-1]/frame_rate) # peak time in seconds

            half = scipy.signal.peak_widths(event_trace, max_index, rel_height=0.5)
            prct_10 = scipy.signal.peak_widths(event_trace, max_index, rel_height=0.1)
            prct_90 = scipy.signal.peak_widths(event_trace, max_index, rel_height=0.9)

            rise_time.append((prct_10[2][0] - prct_90[2][0])/frame_rate) # rise time in seconds
            decay_time.append((prct_90[3][0] - prct_10[3][0])/frame_rate) # decay time in seconds
            half_width.append(half[1][0]/frame_rate) # full width at half maximum in seconds
            duration.append((rpoint - lpoint)/frame_rate) # duration of the signal/event in seconds

            if j_signal == 0:
                inter_event_interval.append(None) # initialize inter-event interval for thr first signal as None
            else:
                inter_event_interval.append((lpoint - end_index[-2])/frame_rate) # calculate inter-event interval from the last event
        
    
    signal_stats = pd.DataFrame({'ROA_ID': ROA_ID, 
                                 'signal_start_index': start_index, 
                                 'signal_end_index': end_index, 
                                 'peak_index': peak_index, 
                                 'peak_time': peak_time, 
                                 'AUC': AUC, 
                                 'amplitude': amplitude, 
                                 'signal_to_noise': signal_to_noise,
                                 'rise_time': rise_time, 
                                 'decay_time': decay_time, 
                                 'half_width': half_width, 
                                 'duration': duration,
                                 'inter_event_interval': inter_event_interval})
    
    # add column to indicate if the signal starts before drug application (True) or after (False)
    if drug_frame != 'NA':
        signal_stats['Drug'] = np.where(signal_stats['signal_start_index'] < drug_frame, 'Before', 'After')
    else:
        signal_stats['Drug'] = 'NA'
                                 
    return signal_stats

def read_tif(tif_path):
    '''
    read a binary tif file and return the labeled matrix and the number of labels

    ROAs are marked as 1 in the binary matrix
    one ROA is defined as a group of connected pixels (including diagonally connected pixels)
    '''

    image = open(tif_path, 'rb').read()
    map_tif = Image.open(io.BytesIO(image))
    map_array = np.asarray(map_tif)
    map_labeled, map_num = scipy.ndimage.label(map_array, structure = [[1,1,1],[1,1,1],[1,1,1]]) # determine ROA and number of ROAs
    return map_labeled, map_num

def align_ROA_cell(ROA_map_labeled, cell_map_labeled, ROA_map_count):
    ''' aline ROA_ID and cell_ID based on the labeled map'''

    ROA = range(1, ROA_map_count+1)
    ROA_cell = []

    for i_ROA in ROA:
        cell_assigned = cell_map_labeled[ROA_map_labeled == i_ROA] # find corresponding cell ID for each pixel inside a ROA
        ROA_cell.append(scipy.stats.mode(cell_assigned, keepdims = False).mode) # assign the most common cell ID for the ROA as its cell registration

    df_ROA_cell = pd.DataFrame({'ROA_ID': ROA, 'cell_ID': ROA_cell})
    return df_ROA_cell

def ROA_analysis(signal_stats, df_ROA_cell):
    '''analyze the signals based on the ROA and return a dataframe with ROA stats (columns) for each individual ROA (rows)'''

    # calculate the signal stats based on ROA
    cols = ['AUC','amplitude','signal_to_noise','rise_time','decay_time','half_width','duration','inter_event_interval']
    ROA_based_count = signal_stats.groupby(['ROA_ID', 'Drug'], as_index = False).count()
    ROA_based = signal_stats.groupby(['ROA_ID', 'Drug'], as_index = False).mean()
    ROA_based['signal_count'] = ROA_based_count['AUC']

    # identify ROA type based on activity before and after drug application
    ROA = range(1, ROA_count + 1)
    ROA_type = []
    for i_ROA in ROA:
        df = ROA_based[ROA_based.ROA_ID == i_ROA]
        
        if len(df) == 2:
            ROA_type.append('stable')
        elif len(df) == 1:
            if int(df.Drug == 'Before'): # force into a int to bypass series ambiguity
                ROA_type.append('off')
            else:
                ROA_type.append('on')
        else:
            ROA_type.append('inactive')
    df_ROA_cell['ROA_type'] = ROA_type
    ROA_based = pd.merge(df_ROA_cell, ROA_based, on = ['ROA_ID', 'cell_ID'], how = 'left')

    return ROA_based, df_ROA_cell

def inspect_trace(ROA_ID, dff_traces, baselines, thresholds):
    
    ''' 
    inspect trace visually with baseline, signal threshold and drug application time indicated
    
    Args:
    ROA_ID: input ROA_ID (int) to visually inspect
    '''
    ROA_count, frame_count = dff_traces.shape
    x = np.arange(frame_count)

    # indexing with ROA_ID - 1 because python...
    plt.figure(figsize=(10,5))
    plt.plot(x,dff_traces[ROA_ID-1],color = 'grey')
    plt.axhline(y = baselines[ROA_ID-1], color = 'r', linestyle = '-', label = "baseline")
    plt.axhline(y = thresholds[ROA_ID-1], color = 'g', linestyle = '-', label = "threshold")
    plt.axvline(x = drug_frame, color = 'b', alpha = 0.5, label = "drug application")
    plt.legend(loc = 'upper left')
    plt.xlabel("Frame")
    plt.ylabel("dF/F")
    plt.title('ROA ID: ' + str(ROA_ID))
    plt.show(block = False)