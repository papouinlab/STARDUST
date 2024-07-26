# util.py
import scipy.signal
import os, io, warnings, math, scipy, numpy as np, pandas as pd, matplotlib, PIL, seaborn as sns
import matplotlib.pyplot as plt
from PIL import Image
 
def check_modules():
    '''
    Check module version for STARDUST. 
    '''
    tracking = 1
    if scipy.__version__ < '1.11.0':
        print('Please update SciPy to 1.11.0 or higher.')
        tracking = 0
    if np.__version__ < '1.26.4':
        print('Please update NumPy to 1.11.0 or higher.')
        tracking = 0
    if pd.__version__ < '2.1.4':
        print('Please update Pandas to 2.1.4 or higher.')
        tracking = 0
    if matplotlib.__version__ < '3.8.0':
        print('Please update Matplotlib to 3.8.0 or higher.')
        tracking = 0
    if sns.__version__ < '0.12.2':
        print('Please update Seaborn to 0.12.2 or higher.')
        tracking = 0
    if PIL.__version__ < '10.2.0':
        print('Please update pillow to 10.2.0 or higher.')
        tracking = 0
    
    if tracking == 1:
        print('All module version requirements are met for STARDUST. You\'re good to go!')
         

def check_path(dir):
    if "/" in dir: # mac
        if dir[-1] != "/":
            dir = dir + "/"
    elif '\\' in dir: # win
        if dir[-1] != "\\":
            dir = dir + "\\"
    return dir
 
def prompt_input(analysis_type = "full"):
    '''
    Prompt user input for file paths and output file name.
    '''
    # input
    if analysis_type == "full":
        input_dir = input("Enter the input folder path: ")
        time_series_filename = input("Enter the time series file name: ")
        ROA_mask_filename = input("Enter the ROA mask file name: ")
        cell_mask_filename = input("Enter the cell mask file name: ")
        
        if '.tif' not in ROA_mask_filename:
            ROA_mask_filename = ROA_mask_filename + '.tif'
        if '.tif' not in cell_mask_filename:
            cell_mask_filename = cell_mask_filename + '.tif'   
    else:
        input_dir = input("Enter the input folder path: ")
        time_series_filename = input("Enter the time series file name: ")
        if analysis_type == "ROA-based":
            ROA_mask_filename = input("Enter the ROA mask file name: ")
            if '.tif' not in ROA_mask_filename:
                ROA_mask_filename = ROA_mask_filename + '.tif'
        elif analysis_type == "cell-based":
            cell_mask_filename = input("Enter the cell mask file name: ")
            if '.tif' not in cell_mask_filename:
                cell_mask_filename = cell_mask_filename + '.tif'
    
    input_dir = check_path(input_dir) # make sure the path ends with '/'
    if '.csv' not in time_series_filename:
        time_series_filename = time_series_filename + '.csv'
    
    
    # check if the files are in the input directory
    input_files = os.listdir(input_dir)
    if analysis_type != "full":
        if time_series_filename not in input_files:
            warnings.warn("Time series file not found in the input directory. Please check the file name and input path.")
        else:
            time_series_path = os.path.join(input_dir, time_series_filename)
        
        if analysis_type == "ROA-based":
            if ROA_mask_filename not in input_files:
                warnings.warn("ROA mask file not found in the input directory. Please check the file name and input path.")
            else:
                ROA_mask_path = os.path.join(input_dir, ROA_mask_filename)
            cell_mask_path = None # placeholder
            
        elif analysis_type == "cell-based":
            if cell_mask_filename not in input_files:
                warnings.warn("Cell mask file not found in the input directory. Please check the file name and input path.")
            else:
                cell_mask_path = os.path.join(input_dir, cell_mask_filename)
            ROA_mask_path = None # placeholder
        
    elif analysis_type == "full":
        if time_series_filename not in input_files:
            warnings.warn("Time series file not found in the input directory. Please check the file name and input path.")
        else:
            time_series_path = os.path.join(input_dir, time_series_filename)
        if ROA_mask_filename not in input_files:
            warnings.warn("ROA mask file not found in the input directory. Please check the file name and input path.")
        else:
            ROA_mask_path = os.path.join(input_dir, ROA_mask_filename)
        if cell_mask_filename not in input_files:
            warnings.warn("Cell mask file not found in the input directory. Please check the file name and input path.")
        else:
            cell_mask_path = os.path.join(input_dir, cell_mask_filename)
    
    # output
    output_dir = input('Enter the path to the folder for output files (Hit Enter if wish to save to the input directory): ')
    output_filename = input('Enter the output file name: ')
    if output_dir == "":
        output_dir = input_dir
    output_path = os.path.join(output_dir, output_filename) 
    
    return time_series_path, ROA_mask_path, cell_mask_path, output_path

def get_metadata():
    # input metadata related to the experiment
    drug_frame = float(input('Enter the frame number of drug application (enter 0 if no drug was applied): '))
    frame_rate = float(input("Enter the frame rate (in Hz) of the recording (e.g. 1): "))
    spatial_resolution = float(input('Enter the spatial resolution of the recording (in µm/pixel, if one pixel size is 0.884 x 0.884 µm^2, enter 0.884): '))
    
    return drug_frame, frame_rate, spatial_resolution

def find_files(input_dir, keyword):
    '''
    Takes an input directory and find the trace csv file, ROA mask tif, cell mask tif.
    '''

    input_files = os.listdir(input_dir)
    input_files = [f.lower() for f in input_files]

    exp_index = [keyword.lower() in f for f in input_files]
    exp_files = np.take(input_files, np.where(exp_index)[0])
    
    # find the time series raw csv file
    csv_index = ['.csv' in f for f in exp_files]
    csv_files = np.take(exp_files, np.where(csv_index)[0])
    time_series_index = ['raw' in f for f in csv_files]

    # find the ROA and cell mask tif files
    tif_index = ['.tif' in f for f in exp_files]
    tif_files = np.take(exp_files, np.where(tif_index)[0])
    ROA_index = ['roa' in f for f in tif_files]
    cell_index = ['cell' in f for f in tif_files]
    
    time_series_path = os.path.join(input_dir, csv_files[time_series_index.index(True)])
    ROA_mask_path = os.path.join(input_dir, tif_files[ROA_index.index(True)])
    cell_mask_path = os.path.join(input_dir, tif_files[cell_index.index(True)])
    
    print("Found the following files: \n")
    print("CSV file: ", time_series_path)
    print("ROA mask: ", ROA_mask_path)
    print("Cell mask: ", cell_mask_path)

    return time_series_path, ROA_mask_path, cell_mask_path

def read_tif(tif_path, type):

    '''
    Read a binary tif file and return the labeled matrix and the number of labels.

    Args:
    tif_path: path to the tif file.
    type: 'ROA' or 'cell' to determine the structure of the labeled matrix.

    Returns:
    map_array: a binary numpy array of the tif file.
    map_labeled: a labeled matrix of the tif file.
    map_count: the number of labels in the labeled matrix.

    ROAs and cells are marked as 1 in the binary matrix. 
    One ROA is defined as a group of connected pixels (including diagonally connected pixels). 
    One cell is defined as a group of connected pixels (NOT including diagonally connected pixels).
    '''
    print("Reading in file for ", type, "mask: ", tif_path, "\n")
    image = open(tif_path, 'rb').read()
    map_tif = Image.open(io.BytesIO(image))
    map_array = np.asarray(map_tif)

    # determine ROA and number of ROAs
    if type == "ROA":
        map_labeled, map_count = scipy.ndimage.label(map_array, structure = [[1,1,1],[1,1,1],[1,1,1]])
        print("This ROA mask contains", str(map_count), "ROAs.")
    elif type == "cell":
        map_labeled, map_count = scipy.ndimage.label(map_array) 
        print("This cell mask contains", str(map_count), "cells.")
    print("\n")

    return map_array, map_labeled, map_count

def visualize_map(ROA_map_array = np.empty(0), cell_map_array = np.empty(0)):

    '''
    Visualizes the ROA and cell masks from binary numpy array. 
    '''
    if ROA_map_array.any() and cell_map_array.any():
        fig, axs = plt.subplots(1,2, figsize = (10, 10))
        axs[0].imshow(ROA_map_array)
        axs[0].set_title('ROA mask')
        axs[1].imshow(cell_map_array)
        axs[1].set_title('Cell mask')
        plt.show()
    elif ROA_map_array.any():
        plt.imshow(ROA_map_array)
        plt.title('ROA mask')
        plt.show()
    elif cell_map_array.any():
        plt.imshow(cell_map_array)
        plt.title('Cell mask')
        plt.show()

def raw_to_filtered(csv_path, order = 4, cutoff = 0.4):
    '''
    Read in and convert raw traces to filtered traces using Butterworth filter. 

    Args:
    csv_path: path to the csv file containing raw traces.
    order: optional, order of the Butterworth filter (default 4).
    cutoff: optional, cutoff frequency of the Butterworth filter in Hz (default 0.4).

    Returns:
    raw_traces: a 2D numpy array of raw traces.
    filtered_traces: a 2D numpy array of filtered traces.
    '''
    
    # read in the csv data file in as a dataframe and transpose
    # the csv file has one row of header (Mean1, Mean2...) and first column is the number of ROA
    # in the original data frame, each row represents a frame and each column represents a ROA
    print("Reading in file: ", csv_path, "\n\n")
    raw_data = pd.read_csv(csv_path, header = 0, index_col = 0)
    if raw_data.columns[0] != "Mean1": # just in case when there's no header... 
        raw_data = pd.read_csv(csv_path, header = None, index_col = 0)
    # transpose data so that each row represents a ROA and each column represents a frame
    # and transform to a numpy 2D array containing # of ROA lists,
    # each list containing signal in each frame for this ROA
    raw_traces = raw_data.transpose().to_numpy() 

    # apply a lowpass Butterworth filter with a 4th order filter at the cutoff of 0.4 Hz
    print("Applying a lowpass Butterworth filter with a", str(order), "th order filter at the cutoff of", str(cutoff), "Hz")
    b, a = scipy.signal.butter(order, cutoff, 'low', analog = False)
    filtered_traces = scipy.signal.filtfilt(b, a, raw_traces)

    return raw_traces, filtered_traces

def check_traces(traces):
    ''' 
    Check the dimension of traces and print out the number of ROAs and the number of frames.
    '''
    
    (ROA_count, frame_count) = traces.shape

    print("The current file contains: ")
    print("Number of ROA: ", ROA_count)
    print("Number of frames: ", frame_count)
    return ROA_count, frame_count

def correct_shift(filtered_traces, correction_factor = 0.5):
    
    '''
    'Correct' the shift of baseline level in traces using linear regression for subsequent F0 determination. 
    Specifically, corrected traces are generated by subtracting correction_factor (default 0.5) * slope at the current x location from the input traces. 
    Prints out the distribution of slopes. 
    Note that this correction is not necessary if the baseline is stable and that this correction should not be applied to recordings when there is obvious z shift.

    Args:
    filtered_traces: a 2D numpy array of filtered traces.
    correction_factor: optional, the factor to multiply the slope to correct the traces (default 0.5).

    Returns:
    corrected_traces: a 2D numpy array of corrected traces.
    reg: linear regression information containing the slope and intercept of the linear regression for each ROA.

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
    corrected_traces = filtered_traces - correction_factor * x * slope[:,None]

    # visualize distribution of lm slope
    sns.histplot(reg.slope).set(title = 'Slope distribution histogram')
    return corrected_traces, reg

def check_correction(filtered_traces, corrected_traces, reg):

    '''
    Visual checkpoint for correction.
    '''

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

def pull_largeslope_traces(ROA_count, reg):
    '''
    Pull out ROA IDs that have user-defined large slope for further inspection. 
    
    Prompted input: reg_threshold (absolute value of the slope cutoff)
    '''

    reg_threshold = float(input("Enter the slope cutoff (absolute value) to pull out for check later:"))
    check_ROAs = []

    for i_ROA in range(ROA_count):
        if abs(reg.slope[i_ROA]) > reg_threshold:
            check_ROAs.append(i_ROA + 1) # ROA ID starts from 1
    
    print(f"The following ROAs have a slope larger than {reg_threshold} or smaller than {-reg_threshold}:", check_ROAs)
    return check_ROAs

def find_roots(trace, threshold):
    '''
    Finds where the trace crosses threshold and returns the x axis intercepts here the signal crosses threswhold.

    Args:
    trace: a 1D numpy array of trace.
    threshold: the threshold to find the intercepts.

    Returns:
    x0: a 1D numpy array of intercepts.
    '''

    trace_adjusted = trace - threshold
    cross_bool = np.abs(np.diff(np.sign(trace_adjusted))).astype(bool)
    
    y1 = trace_adjusted[:-1][cross_bool]
    y2 = trace_adjusted[1:][cross_bool]
    x1 = np.where(cross_bool)[0]
    x0 = x1 - y1/(y2-y1) # solve for the x axis intercept (AKA when y=0)
    return x0

def find2points(x, array):
    '''
    Given a specific number x and an array containing numbers (arranged from smallest to largest), find the closest two flanking numbers of x inside the array.
    '''
    for index in range(len(array)):
        if x - array[index] < 0:
            lpoint = array[index-1]
            rpoint = array[index]
            break
    return (lpoint, rpoint)

def find_signal_boundary(trace, signal_threshold, baseline_threshold, include_end_incomplete = False):
    ''' 
    Takes a 1D trace, signal_threshold, and baseline_threshold (to call beginning and end of an event) and finds tuples of event start and end points (frames). 

    Args:
        trace: a 1D numpy array of trace.
        signal_threshold: the threshold to find the signal.
        baseline_threshold: the threshold to find the baseline.
        include_end_incomplete: optional, include incomplete signals at the end of the trace, default False.

    Returns:
        signal_boundary: a list of tuples of event start and end points.
        signal_intercept: a 1D numpy array of signal intercepts.
        baseline_intercept: a 1D numpy array of baseline intercepts.
    '''
    # find intercepts for signal and baseline
    signal_intercept = find_roots(trace, signal_threshold)
    baseline_intercept = find_roots(trace, baseline_threshold)

    signal_list = [] # stores the start and end points of each signal event
    for i in signal_intercept:
        if i > baseline_intercept[0] and i < baseline_intercept[-1]: 
            (lpoint, rpoint) = find2points(i,baseline_intercept)
            if (lpoint, rpoint) not in signal_list:
                signal_list.append((lpoint, rpoint))
        elif include_end_incomplete == True and i > baseline_intercept[-1]: # incomplete signal at the end of the trace
            (lpoint, rpoint) = (baseline_intercept[-1], len(trace)) 
            if (lpoint, rpoint) not in signal_list:
                signal_list.append((lpoint, rpoint))

    # find frame numbers for start and end of signals
    #signal_boundary = [(math.floor(lpoint),math.ceil(rpoint)) for (lpoint,rpoint) in signal_list] 
    signal_boundary = signal_list
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
    
    # calculate the baseline averages from signal-removed traces (ignoring nan) and generate dF/F traces
    baselines = np.nanmean(traces_nosignal[:,baseline_start:baseline_end], axis = 1) 
    dff_traces = (traces - baselines[:,None])/baselines[:,None] 

    dff_traces_nosignal = np.multiply(dff_traces, signal_frames == False) # set signal frames to 0
    dff_traces_nosignal[dff_traces_nosignal == 0] = np.nan # change 0 to nan to remove signal frames in calculation

    return dff_traces, dff_traces_nosignal

def find_signal_frames(filtered_traces, signal_frames = None, signal_threshold = 3, baseline_start = 0, baseline_end = -1, include_incomplete = False):
   
    '''
    Finds signals in the filtered traces using signal_threshold and onset_threshold.
    Signal is defined as any event that is larger than signal_threshold * SD of the dF/F0 trace (default 3SD).

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
    baselines = np.nanmean(dff_traces_nosignal[:,baseline_start:baseline_end], axis = 1) # calculate the baseline averages from signal-removed dff traces (ignoring nan)
    thresholds = signal_threshold * np.nanstd(dff_traces_nosignal[:,baseline_start:baseline_end], axis = 1) # calculate the baseline standard deviations from signal-removed dff traces (ignoring nan)
    
    # initialize a new array to store signal frames
    ROA_count, frame_count = filtered_traces.shape
    new_signal_frames = np.zeros((ROA_count, frame_count))
    signal_boundaries = []

    # iterate through each ROA to find signals
    for i_ROA in range(0, ROA_count):
        signal_boundary, _ , _ = find_signal_boundary(dff_traces[i_ROA,], thresholds[i_ROA], baselines[i_ROA], include_incomplete)
        signal_boundaries.append(signal_boundary)

        for j_frame in range(0, frame_count):
            for (l,r) in signal_boundary:
                if j_frame >= l and j_frame <= r:
                    new_signal_frames[i_ROA, j_frame] = 1
                    
    return dff_traces, baselines, thresholds, new_signal_frames.astype(bool), signal_boundaries

def iterative_baseline(traces, signal_frames = None, baseline_start = 0, baseline_end = -1, include_incomplete = False):

    '''
    Determine the bassline iteratively.

    Args:
    traces: a two-dimensional array of traces
    signal_frames: optional, a two-dimensional boolean array corresponding to each ROA and each frame, true if frame is considered as signal
    baseline_start: optional, the starting frame for baseline calculation, default 0 (beginning of the recording)
    baseline_end: optional, the ending frame for baseline calculation, default -1 (end of the recording)
    include_incomplete: optional, include incomplete signals at the end of the trace, default False

    Prompted inputs:
    n_iteration: number of iterations for signal detection.
    signal_threshold: signal_threshold * dF/F baseline SD as signal threshold.

    Returns:
    dff_traces: a two-dimensional array of delta F/F trace based on the provided thresholds.
    baselines: a one-dimensional array of baseline values.
    thresholds: a one-dimensional array of signal thresholds.
    signal_frames: a two-dimensional boolean array corresponding to each ROA and each frame, true if frame is considered as signal
    signal_boundaries: a list of tuples of event start and end points of each detected activity in each ROA
    '''

    signal_frames_previous = signal_frames
    n_iteration = int(input("Enter the number of iterations for signal detection. We found that most baseline determination stabilize after 5 iterations: "))
    signal_threshold = float(input("Enter the signal threshold for signal detection. Trace fractions where the fluorescence is larger than signal_threshold * SD (suggested value 2-3) are determined as active signal: "))
    print(f"Using signal threshold of {signal_threshold}* SD and detecting baseline from frame {baseline_start} to {baseline_end}.\n")

    for iter in range(n_iteration):
        print(f"Processing round {iter+1} of signal detection...")
        if iter < n_iteration - 1: # only store signal_frames for the initial iterations
            _ , _, _, signal_frames, _ = find_signal_frames(traces, signal_frames_previous, signal_threshold, baseline_start, baseline_end, include_incomplete)
        else:
            dff_traces, baselines, thresholds, signal_frames, signal_boundaries = find_signal_frames(traces, signal_frames_previous, signal_threshold, baseline_start, baseline_end, include_incomplete)
        
        check_ROA(signal_frames) # check current signal detection results
        signal_frames_previous = signal_frames
        print()
    
    return dff_traces, baselines, thresholds, signal_frames, signal_boundaries, signal_threshold

def check_ROA(signal_frames):
    '''check and report how many ROAs have signals'''
    ROAs_with_signal = np.sum(signal_frames, axis = 1) > 0

    print("ROAs with signal: ", np.sum(ROAs_with_signal))
    print("ROAs without signal: ", np.sum(ROAs_with_signal == False))

def analyze_signal(dff_traces, signal_frames, signal_boundaries, frame_rate, drug_frame):

    '''
    Analyze the signals and return a dataframe with signal stats (columns) for each individual signal/event (rows).
    
    Args:
    dff_traces: a two-dimensional array of delta F/F trace.
    signal_frames: a two-dimensional boolean array corresponding to each ROA and each frame, true if frame is considered as signal.
    signal_boundaries: a list of tuples of event start and end points of each detected activity in each ROA. 
    frame_rate: the frame rate of the recording.
    drug_frame: the frame number of drug application.
    
    Returns:
    signal_features: a dataframe with signal stats/features for each individual signal/event.
    '''

    (ROA_count, frame_count) = dff_traces.shape 
    noise = np.multiply(dff_traces, signal_frames == False).max(axis = 1) # extract baseline noise (max amplitude within baseline) for each ROA

    # initialize lists to store signal stats
    ROA_ID = []
    start_frame = []
    start_time = []
    end_frame = []
    end_time = []
    AUC = []
    amplitude = []
    signal_to_noise = []
    peak_frame = []
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
            
            if rpoint + 1 > frame_count: # if the current signal is an incomplete signal at the end of the trace
                
                start_frame.append(lpoint)
                start_time.append(lpoint/frame_rate)
                end_frame.append(np.nan)
                end_time.append(np.nan)
                event_trace = dff_traces[i_ROA,math.ceil(lpoint):]

                AUC.append(np.nan) # AUC not calculatable
                amplitude.append(max(event_trace))  # signal amplitude
                signal_to_noise.append(max(event_trace)/noise[i_ROA]) # signal to noise ratio

                max_index = np.array([np.argmax(event_trace)]) # max dff index/frame number within the signal/event range
                peak_frame.append(lpoint + max_index[0]) # max dff index (frame number) within the whole trace 
                peak_time.append(peak_frame[-1]/frame_rate) # peak time in seconds

                prct_10 = scipy.signal.peak_widths(event_trace, max_index, rel_height=0.1)
                prct_90 = scipy.signal.peak_widths(event_trace, max_index, rel_height=0.9)

                rise_time.append((prct_10[2][0] - prct_90[2][0])/frame_rate) # rise time in seconds
                decay_time.append(np.nan) # decay time in seconds
                half_width.append(np.nan) # full width at half maximum in seconds
                duration.append(np.nan) # duration of the signal/event in seconds
            
            else:
                
                start_frame.append(lpoint)
                start_time.append(lpoint/frame_rate)
                end_frame.append(rpoint)
                end_time.append(rpoint/frame_rate)

                event_trace = dff_traces[i_ROA,math.ceil(lpoint):math.floor(rpoint)+1] # subset out only the signal/event

                AUC.append(scipy.integrate.simpson(event_trace, dx = 1/frame_rate)) # area under the curve using Simpson's rule
                amplitude.append(max(event_trace))  # signal amplitude
                signal_to_noise.append(max(event_trace)/noise[i_ROA]) # signal to noise ratio

                max_index = np.array([np.argmax(event_trace)]) # max dff index/frame number within the signal/event range
                peak_frame.append(lpoint + max_index[0]) # max dff index (frame number) within the whole trace 
                peak_time.append(peak_frame[-1]/frame_rate) # peak time in seconds

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
                inter_event_interval.append((lpoint + 1 - end_time[-2])/frame_rate) # calculate inter-event interval from the last event
        
    
    signal_features = pd.DataFrame({'ROA_ID': ROA_ID, 
                                    'signal_start_frame': start_frame, 
                                    'signal_start_time': start_time,
                                    'signal_end_frame': end_frame, 
                                    'signal_end_time': end_time,
                                    'peak_frame': peak_frame, 
                                    'peak_time': peak_time, 
                                    'AUC': AUC, 
                                    'amplitude': amplitude, 
                                    'signal_to_noise': signal_to_noise,
                                    'rise_time': rise_time, 
                                    'decay_time': decay_time, 
                                    'half_width': half_width, 
                                    'duration': duration,
                                    'inter_event_interval': inter_event_interval})
    
    # add column to indicate if the signal peaks before drug application (baseline) or after (drug)
    if drug_frame == 0:
        signal_features['epoch'] = 'NA'
    else:
        signal_features['epoch'] = np.where(signal_features['peak_frame'] < drug_frame, 'baseline', 'drug')
                                 
    return signal_features

def align_ROA_cell(ROA_map_labeled, cell_map_labeled, ROA_map_count, spatial_resolution):
    ''' 
    Aline ROA_ID and cell_ID based on the labeled map and return a dataframe ROA_info with ROA and cell alignment.
    '''

    ROA = range(1, ROA_map_count+1)
    ROA_cell = []
    ROA_not_assigned = []
    size_pixel = []

    for i_ROA in ROA:
        # find the size of the ROA in number of pixels
        i_size = np.sum(ROA_map_labeled == i_ROA) 
        size_pixel.append(i_size)
        
        # find the most common cell ID assigned to the ROA
        cell_assigned = cell_map_labeled[ROA_map_labeled == i_ROA] # find corresponding cell ID for each pixel inside a ROA
        most_frequent = scipy.stats.mode(cell_assigned, keepdims = False).mode # find the most common cell ID inside the ROA
        if most_frequent == 0: # if the most common cell ID is 0, then find the second most common cell ID
            if len(np.unique(cell_assigned[cell_assigned != 0])) != 0:
                most_frequent = scipy.stats.mode(cell_assigned[cell_assigned != 0], keepdims = False).mode
            else:
                ROA_not_assigned.append(i_ROA)
        ROA_cell.append(most_frequent) # assign the most common cell ID for the ROA as its cell registration

    ROA_info = pd.DataFrame({'ROA_ID': ROA, 'cell_ID': ROA_cell, 'size_pixel' : size_pixel})
    ROA_info['size_um2'] = ROA_info['size_pixel'] * spatial_resolution**2

    if len(ROA_not_assigned) != 0:
        print(f"There are {len(ROA_not_assigned)} ({round(len(ROA_not_assigned)/ROA_map_count * 100,2)}%) ROAs not assigned to any cell.")
        print("ROA IDs:", ROA_not_assigned)
        print("Please double check the cell mask registration. \n")
    print("ROA and cell alignment completed.")
    return ROA_info

def size_ROA(ROA_map_labeled, ROA_map_count, spatial_resolution):
    '''
    Extract size information in pixel and convert to um^2 for each ROA. 
    Use this function instead of align_ROA_cell when doing ROA-based analysis (without cell segmentation). 
    '''
    
    ROA = range(1, ROA_map_count+1)
    size_pixel = []
    
    for i_ROA in ROA:
        i_size = np.sum(ROA_map_labeled == i_ROA)
        size_pixel.append(i_size)
    ROA_info = pd.DataFrame({'ROA_ID': ROA, 'size_pixel' : size_pixel})
    ROA_info['size_um2'] = ROA_info['size_pixel'] * spatial_resolution**2
    
    return ROA_info

def ROA_analysis(signal_stats, ROA_info, frame_count, frame_rate, drug_frame, cell_segmentation = False):
    '''
    Analyze signals based on ROA and return a dataframe with ROA stats (columns) for each individual ROA (rows).
    
    Args:
    signal_stats: 
    '''

    # calculate the signal features based on ROA
    ROA_based_count = signal_stats.groupby(['ROA_ID', 'epoch'], as_index = False).count()
    ROA_based = signal_stats.groupby(['ROA_ID', 'epoch'], as_index = False).mean()
    ROA_based['signal_count'] = ROA_based_count['AUC']

    # identify ROA type (inactive, stable, on, off, NA) based on activity during baseline and after drug application
    ROA = ROA_info.ROA_ID
    ROA_type = []

    for i_ROA in ROA:
        df = ROA_based[ROA_based.ROA_ID == i_ROA]
        
        if 'NA' in df.epoch.unique():
            ROA_type.append('NA')
        elif 'baseline' in df.epoch.unique() and 'drug' in df.epoch.unique():
            ROA_type.append('stable')
        elif 'baseline' in df.epoch.unique():
            ROA_type.append('off')
        elif 'drug' in df.epoch.unique():
            ROA_type.append('on')
        else:
            ROA_type.append('inactive')

    ROA_info['ROA_type'] = ROA_type
    if cell_segmentation == True:
        ROA_based = pd.merge(ROA_based, ROA_info[['ROA_ID', 'cell_ID', 'ROA_type', 'size_um2']], 
                            on = ['ROA_ID', 'cell_ID'], how = 'right')
    else:
        ROA_based = pd.merge(ROA_based, ROA_info[['ROA_ID', 'ROA_type', 'size_um2']],
                            on = ['ROA_ID'], how = 'right')

    ROA_based['signal_count'] = np.where(ROA_based['ROA_type'] == 'inactive', 0, ROA_based['signal_count'])

    # calculate recording total, baseline and drug length (in minutes) for frequency calculation
    total_length_min = frame_count/(frame_rate*60)
    baseline_length_min = drug_frame/(frame_rate*60)
    drug_length_min = total_length_min - baseline_length_min

    ROA_based['rec_length'] = np.nan
    ROA_based['rec_length'] = np.where(ROA_based['epoch'] == 'NA', total_length_min, ROA_based['rec_length'])
    ROA_based['rec_length'] = np.where(ROA_based['ROA_type'] == 'inactive', total_length_min, ROA_based['rec_length'])
    ROA_based['rec_length'] = np.where(ROA_based['epoch'] == 'baseline', baseline_length_min, ROA_based['rec_length'])
    ROA_based['rec_length'] = np.where(ROA_based['epoch'] == 'drug', drug_length_min, ROA_based['rec_length'])
    ROA_based['frequency_permin'] = ROA_based['signal_count']/ROA_based['rec_length']   
    ROA_based = pd.merge(ROA_based, ROA_info[['ROA_ID', 'size_um2']], on = 'ROA_ID', how = 'right')
    
    if cell_segmentation:
        cols = ['ROA_ID', 'cell_ID', 'ROA_type', 'size_um2', 'epoch',
                'AUC','amplitude','signal_to_noise','rise_time','decay_time','half_width','duration',
                'inter_event_interval', 'signal_count', 'frequency_permin']
    else:
        cols = ['ROA_ID', 'ROA_type', 'size_um2', 'epoch',
                'AUC','amplitude', 'signal_to_noise', 'rise_time', 'decay_time', 'half_width', 'duration',
                'inter_event_interval', 'signal_count', 'frequency_permin']
    return ROA_based[cols], ROA_info

def ROA_type_summary(ROA_info):
    '''
    Count the number of ROAs that belongs to each category (inactive, stable, on, off, and NA). 
    '''
    
    ROA_sub = ROA_info[['ROA_ID', 'ROA_type']]
    ROA_summary = ROA_sub.groupby('ROA_type').count()
    ROA_summary.rename(columns = {'ROA_ID':'count'}, inplace = True)
    ROA_summary['percentage'] = ROA_summary['count']/ROA_summary['count'].sum() * 100
    
    return ROA_summary

def cell_analysis(signal_stats, ROA_info):
    '''
    Analyze the signal stats at the cell level. Cell-assigned ROA analysis. 
    '''
    cell_based_ROA_count = ROA_info.groupby(['cell_ID'], as_index = False).count()
    cell_based_ROA_count = cell_based_ROA_count[['cell_ID', 'ROA_ID']]
    cell_based_ROA_count.columns = ['cell_ID', 'ROA_count']
    cell_based_signal_count = signal_stats.groupby(['cell_ID', 'epoch'], as_index = False).count()
    cell_based = signal_stats.groupby(['cell_ID', 'epoch'], as_index = False).mean()
    cell_based['signal_count'] = cell_based_signal_count['ROA_ID']

    cols = ['cell_ID', 'epoch', 'signal_count', 'AUC', 'amplitude', 'signal_to_noise', 'rise_time', 'decay_time', 'half_width', 'duration']
    cell_based = pd.merge(cell_based_ROA_count, cell_based[cols], on = 'cell_ID', how = 'right')
    return cell_based

def cell_based(signal_features, frame_count, frame_rate, drug_frame):
    
    '''
    Calculate the signal stats at the cell level for cell-based analysis. 
    '''
    
    cell_based_count = signal_features.groupby(['ROA_ID', 'epoch'], as_index = False).count()
    cell_based = signal_features.groupby(['ROA_ID', 'epoch'], as_index = False).mean()
    cell_based['signal_count'] = cell_based_count['AUC']
    cell_based['cell_ID'] = cell_based['ROA_ID']
    
    # calculate recording total, baseline and drug length (in minutes) for frequency calculation
    total_length_min = frame_count/(frame_rate*60)
    baseline_length_min = drug_frame/(frame_rate*60)
    drug_length_min = total_length_min - baseline_length_min
    
    cell_based['rec_length'] = np.nan
    cell_based['rec_length'] = np.where(cell_based['epoch'] == 'NA', total_length_min, cell_based['rec_length'])
    cell_based['rec_length'] = np.where(cell_based['epoch'] == 'baseline', baseline_length_min, cell_based['rec_length'])
    cell_based['rec_length'] = np.where(cell_based['epoch'] == 'drug', drug_length_min, cell_based['rec_length'])
    cell_based['frequency_permin'] = cell_based['signal_count']/cell_based['rec_length']
    
    cols = ['cell_ID', 'epoch', 'signal_count', 'AUC', 'amplitude', 'signal_to_noise', 'rise_time', 'decay_time', 'half_width', 'duration',
            'inter_event_interval', 'frequency_permin']
    return cell_based[cols]

def inspect_trace(ROA_IDs, dff_traces, baselines, thresholds, drug_frame):
    
    ''' 
    Inspect trace visually with baseline, signal threshold and drug application time indicated.
    
    Args:
        ROA_IDs: input ROA_IDs (iterable) to visually inspect.
        dff_traces: 2D array of dF/F traces.
        baselines: 1D array of baseline values.
        thresholds: 1D array of signal thresholds.
        drug_frame: frame number of drug application (int).
    '''
    frame_count = dff_traces.shape[1]
    x = np.arange(frame_count)

    for i_ROA in ROA_IDs:
        
        plt.figure(figsize=(10,5))
        plt.plot(x,dff_traces[i_ROA-1],color = 'grey') # indexing with i_ROA - 1 because python...
        plt.axhline(y = baselines[i_ROA-1], color = 'r', linestyle = '-', label = "baseline")
        plt.axhline(y = thresholds[i_ROA-1], color = 'g', linestyle = '-', label = "threshold")
        if drug_frame != 0:
            plt.axvline(x = drug_frame, color = 'b', alpha = 0.5, label = "drug application")
        plt.legend(loc = 'upper left')
        plt.xlabel("Frame")
        plt.ylabel("dF/F")
        plt.title('ROA ID: ' + str(i_ROA))
        plt.show(block = False)

def output_data(output_path, metadata, dff_traces, signal_features, save_as = 'csv', 
                ROA_based = None, ROA_info = None, ROA_summary = None, cell_based = None, cell_assigned = None):

    '''
    Save the output dataframes to csv or excel.
    '''
    
    if save_as.lower() == 'csv':
        print("Saving outputs as csv files...")
        metadata.to_csv(output_path + '_metadata.csv', index = False)
        np.savetxt(output_path +'_dff_traces.csv', dff_traces, delimiter=",")
        signal_features.to_csv(output_path + '_signal_features.csv', index = False)
        if ROA_based is not None:
            ROA_based.to_csv(output_path + '_ROA_based.csv', index = False)
        if ROA_info is not None:
            ROA_info.to_csv(output_path + '_ROA_cell_key.csv', index = False)
        if ROA_summary is not None:
            ROA_summary.to_csv(output_path + '_ROA_type_summary.csv', index = True)
        if cell_based is not None:
            cell_based.to_csv(output_path + '_cell_based.csv', index = False)
        if cell_assigned is not None:
            cell_assigned.to_csv(output_path + '_cell_assigned_ROA_analysis.csv', index = False)
    
    elif save_as.lower() == 'excel':
        print("Saving outputs as an excel file...")
        with pd.ExcelWriter(output_path + '.xlsx') as writer:
            metadata.to_excel(writer, sheet_name = 'metadata')
            signal_features.to_excel(writer, sheet_name ='signal features')
            if ROA_based is not None:
                ROA_based.to_excel(writer, sheet_name = 'ROA based')
            if ROA_info is not None:
                ROA_info.to_excel(writer, sheet_name = 'ROA cell key')
            if ROA_summary is not None:
                ROA_summary.to_excel(writer, sheet_name = 'ROA type summary')
            if cell_based is not None:
                cell_based.to_excel(writer, sheet_name = 'cell based')
            if cell_assigned is not None:
                cell_assigned.to_excel(writer, sheet_name = 'cell-assigned ROA analysis')
    else:
        Warning('Invalid file format. Please choose csv or excel.')
    
    print("Outputs saved to" + output_path)
