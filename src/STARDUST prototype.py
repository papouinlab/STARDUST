import io
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal, interpolate
from numpy import mean, absolute, var
import math
from scipy.signal import find_peaks
import pandas as pd
from scipy.integrate import simpson
from scipy import stats
from scipy.ndimage import label
from PIL import Image
import os
import seaborn as sns

#input three data files, spefify the input file name 
data_path = os.path.join(os.path.dirname(__file__),'..','Experimental files','test.csv')
ROAmask_path = os.path.join(os.path.dirname(__file__),'..','Experimental files','test_ROA_Mask.tif')
Cellmask_path = os.path.join(os.path.dirname(__file__),'..','Experimental files','test_Cell_Mask.tif')
input_dir, input_filename = os.path.split(data_path) 


frame_rate = 1 # sampling rate (Para_0)
drug_frame = 0 # The frame when drug application (Para_1)
drug_Time = drug_frame*frame_rate


#create output file name based on the input fileyyyy
opFile = input_filename.split('.')[0] + '_Result.xlsx'
output_path = os.path.join(input_dir,opFile) #output file is saved to the same folder as input files
sheet1 = 'Mastersheet data_Signal'
sheet2 = 'ROA_summary'
sheet3 = "ROA_based details"
sheet4 = 'Cell_based details'
sheet5 = "ROA trace"
sheet6 = 'Cell trace'
sheet7 = 'ROA align to drug trace'
sheet8 = 'Cell align to drug trace'
sheet9 = 'Inter-signal intervals'

#find the total ROI number in this recording 
#find the frame number in this recording
TimeSeriesRawData = np.loadtxt(data_path,delimiter = ',')
num_rows, num_columns = np.shape(TimeSeriesRawData)
ROI_TotalNumber = num_columns -1
framenumber = num_rows
frameIndex =framenumber - 1 


#Store the dataset into a 2D array

Total_Traces = []
for i in range(ROI_TotalNumber):
    Total_Traces.append(TimeSeriesRawData[:,i+1])


#Select the frame to be analyzed (Para_2)
parTotal_Traces = []
for i in range(len(Total_Traces)):
    parTotal_Traces.append(Total_Traces[i][:]) #e.g. to analyze the first 300 frames (Total_Traces[i][0:300])

framenumber = len(parTotal_Traces[0])
print("The frame number is",framenumber)
frameIndex = framenumber -1

# filtered_traces = parTotal_Traces
#Define a low pass filter and apply it to the raw trace
times = framenumber
b, a = signal.butter(4,0.4,'lowpass') #(Para_3)
filtered_traces = signal.filtfilt(b,a, parTotal_Traces)


#if there is no obvious shift in the recording, don't need to use this correction. try to be conservative when using this correction.
#correct the baseline for shifted signals
Total_correctline = []
for i in range(ROI_TotalNumber):
    x = np.arange(framenumber)
    regress_line = stats.linregress(x,filtered_traces[i])
    corrected_line = filtered_traces[i]-0.5*x*regress_line.slope #(Para_4, default is 0.5, 0<= factor <=1)
    Total_correctline.append(corrected_line)
    # Checkpoint for correction
#     # plt.figure()
#     # plt.plot(x,filtered_traces[i],'b-', label = 'orignal trace')
#     # plt.plot(x,x*regress_line.slope+regress_line.intercept,'r-',label ='regression line')
#     # plt.plot(x,corrected_line,'g-',label = 'corrected trace')
#     # plt.xlabel("Time(s)")
#     # plt.ylabel("Fluorescence intensity")
#     # plt.legend()
#     # plt.show()

filtered_traces = Total_correctline


#Calculate the absolute standard deviation of the traces 
def AbsoluteMeanD(ROA_OI):
    AMD = mean(absolute(ROA_OI-mean(ROA_OI)))
    return AMD

def stdev(ROA_OII):
    var_of_ROA = var(ROA_OII)
    stdev = math.sqrt(var_of_ROA)
    return stdev

#find the cross points of the trace with a spefified horizontal line

def find_roots(x,y):
    s = np.abs(np.diff(np.sign(y))).astype(bool)
    return x[:-1][s] + np.diff(x)[s]/(np.abs(y[1:][s]/y[:-1][s])+1)

#find the two most close points
def find2points(number, array):
    global lpoint,rpoint
    for index in range(len(array)):
        if number - array[index] < 0:
            lpoint = array[index-1]
            rpoint = array[index]
            break
    return lpoint, rpoint            



FinalOutput = []
NumSigROA = []
TraceHasSignal = []
Noise = []
Total_dff_line = []
Total_AlignTrace = []
TraceCount = 0
IEI = []

# Trace analysis section
for T in range(1, ROI_TotalNumber+1): # Para_5: select the traces for analysis
    Noise_trace = []
    Trace_number = T
    Trace_idx = Trace_number - 1
    baseline = mean(filtered_traces[Trace_idx][:]) # Para_6: Set the segment for baseline estimation, e.g. [0:200], the goal is to get an approximation of the baseline
    dff_line = (filtered_traces[Trace_idx]-baseline)/baseline
    x = np.arange(framenumber)
    Timesof_AMD = 3 * AbsoluteMeanD(dff_line) # Either absolute mean of deviation or standard deviation can be used to estimate threhold. Para_7: 2 or 3 times of SD/AMD are usually used.
    Timesof_STDEV = 3 * stdev(dff_line)


    #use the find roots functions to identify the intersection of several lines and the traces 

    dff_baseline =0.25* mean(dff_line[:])
    #dff_baseline = 0.25*mean(dff_line[:]) + 0.1 # for manual adjustment to find the peak if needed

    z1 = find_roots(x,dff_line-Timesof_STDEV)
    z2 = find_roots(x,dff_line-dff_baseline)

    if len(z1) == 0:
        print("there is no signal in this ROA in the first round")
        Total_dff_line.append(dff_line)
        AlignTrace = np.roll(dff_line,100-drug_frame) #Para_8: Roll the trace for certain farmes to align the drug timepoints across experiments, e.g. 100 is the drug timepoint after rolling
        Total_AlignTrace.append(AlignTrace)
        continue
    


    closest_points = []
    for i in z1:
        if i > z2[0] and i < z2[-1]:
            m = find2points(i,z2)
            closest_points.append(m)


    Signal_Segment = []
    for x in closest_points:
        if x not in Signal_Segment:
            Signal_Segment.append(x) 


    if len(Signal_Segment) == 0:
        print("there is no signal in this ROA")
        Total_dff_line.append(dff_line)
        AlignTrace = np.roll(dff_line,100-drug_frame) # set the drug_frame at 100s in the final output trace
        Total_AlignTrace.append(AlignTrace)
        continue


    Total_Signal = Signal_Segment
    print("The Total_Signal is:",Total_Signal)


    Total_Signal_Round = []
    for i in Total_Signal:
        Signal_Round = []
        signalL = int(i[0])-1
        signalR = int(i[1])+1
        if signalL < 0:
                signalL == 0
        if signalR > framenumber:
            signalR = framenumber
        Signal_Round = ([signalL,signalR])
        Total_Signal_Round.append(Signal_Round)

    print("The total signal round is", Total_Signal_Round)
    amplitude = []
    halfwidth = []
    count = 0

    xseg_total = []

    for i in Total_Signal_Round:
        a = i[0]
        b = i[1]
        if b == framenumber:
            b = frameIndex
        x = dff_line[a+1:b+1]
        xseg = np.arange(a+1,b+1)
        xseg_total.append(xseg)


    xseg_total = np.concatenate(xseg_total)
    full_array = np.arange(1,framenumber+1)
    truncated_array = full_array.copy()

    for segment in xseg_total:
        segment_array = np.array(segment)
        mask = ~np.all(truncated_array.reshape(-1, 1) == segment_array, axis=1)
        truncated_array = truncated_array[mask]


    split_indices = np.where(np.diff(truncated_array) != 1)[0] + 1
    output_array = np.split(truncated_array, split_indices)


    '''Loop iteration'''
    iteration_times = 5 #Para_9
    for Time in range(iteration_times-1):
        real_baseline = output_array
        baselineSum = []
        for i in real_baseline: 
            a = i[0]
            b = i[-1]
            for p in filtered_traces[Trace_idx][a:b]:
                baselineSum.append(p)
        baseline = mean(baselineSum)
        print("the baseline is",baseline)
        dff_line = (filtered_traces[Trace_idx]-baseline)/baseline
        x = np.arange(framenumber)
        baselineSumFiltered = []
        for i in real_baseline:
            a = i[0]
            b = i[-1]
            for r in dff_line[a:b]:
                baselineSumFiltered.append(r)
        
        dff_baseline = mean(baselineSumFiltered)
        Timesof_AMD = 3 * AbsoluteMeanD(baselineSumFiltered)
        Timesof_STDEV = 3 * stdev(baselineSumFiltered)
        z1 = find_roots(x,dff_line-Timesof_STDEV)
        z2 = find_roots(x,dff_line-dff_baseline)

        if len(z1) == 0:
            print("there is no signal in this ROA in this iteration")
            continue

        z2 = np.sort(z2)

        closest_points = []
        for i in z1:
            if z2[0] < i < z2[-1]:
                m = find2points(i,z2)
                closest_points.append(m)
        
        #exclude the replicate
        Signal_Segment = []
        for x in closest_points:
            if x not in Signal_Segment:
                Signal_Segment.append(x) 

        if len(Signal_Segment) == 0:
            print("there is no signal in this ROA after excluding replicates")
            continue

        Total_Signal = Signal_Segment


        Total_Signal_Round = []
        for i in Total_Signal:
            Signal_Round = []
            signalL = int(i[0])-1
            signalR = int(i[1])+1
            if signalL < 0:
                    signalL == 0
            if signalR > framenumber:
                signalR = framenumber
            Signal_Round = ([signalL,signalR])
            Total_Signal_Round.append(Signal_Round)

        amplitude = []
        halfwidth = []
        count = 0
        xseg_total = []

        for i in Total_Signal_Round:
            a = i[0]
            b = i[1]
            if b == framenumber:
                b = frameIndex
            x = dff_line[a+1:b+1]
            xseg = np.arange(a+1,b+1)
            xseg_total.append(xseg)


        xseg_total = np.concatenate(xseg_total)
        full_array = np.arange(1,framenumber+1)
        truncated_array = full_array.copy()

        for segment in xseg_total:
            segment_array = np.array(segment)
            mask = ~np.all(truncated_array.reshape(-1, 1) == segment_array, axis=1)
            truncated_array = truncated_array[mask]


        split_indices = np.where(np.diff(truncated_array) != 1)[0] + 1
        output_array = np.split(truncated_array, split_indices)



    '''Final iteration'''

    
    real_baseline = output_array
    baselineSum = []
    for i in real_baseline:
        a = i[0]
        b = i[-1]
        for p in filtered_traces[Trace_idx][a:b]:
            baselineSum.append(p)

    baseline = mean(baselineSum)
    print("the baseline is",baseline)
    dff_line = (filtered_traces[Trace_idx]-baseline)/baseline
    x = np.arange(framenumber)
    #determine the dff_line baseline
    baselineSumFiltered = []
    for i in real_baseline:
        a = i[0]
        b = i[-1]
        for r in dff_line[a:b]:
            baselineSumFiltered.append(r)

    dff_baseline = mean(baselineSumFiltered)
    Noise_trace.append(max(baselineSumFiltered))
    Timesof_AMD = 3 * AbsoluteMeanD(baselineSumFiltered)
    Timesof_STDEV = 3 * stdev(baselineSumFiltered)

    z1 = find_roots(x,dff_line-Timesof_STDEV)
    z2 = find_roots(x,dff_line-dff_baseline)

    if len(z1) == 0:
        print("there is no signal in this ROA after all iterations")
        continue
    

    z2 = np.sort(z2)

    closest_points = []
    for i in z1:
        if z2[0] < i < z2[-1]:
            m = find2points(i,z2)
            closest_points.append(m)

    #exclude the replicate
    Signal_Segment = []
    for x in closest_points:
        if x not in Signal_Segment:
            Signal_Segment.append(x) 
    #print('The signal segment is: ',Signal_Segment)

    if len(Signal_Segment) == 0:
        print("there is no signal in this ROA after excluding replicates")
        continue

    Total_Signal = Signal_Segment
    print("The Total Signal after second iteration is:",Total_Signal)
    IEI_trace = []
    if len(Total_Signal)>1:
        for i in range (len(Total_Signal)-1):
            IEI_trace.append(Total_Signal[i+1][0]-Total_Signal[i][1])
    else:
        IEI_trace = 'NaN'
    print("The interval for this trace is: ", IEI_trace)
    IEI.append(IEI_trace)



    Total_Signal_Round = []
    for i in Total_Signal: 
        signalL = math.floor(i[0])
        signalR = math.ceil(i[1])
        if signalL < 0:
            signalL == 0
        if signalR > framenumber:
            signalR == framenumber
        Signal_Round = ([signalL,signalR])
        Total_Signal_Round.append(Signal_Round)

    print('The total signal bound after second iteration is: ', Total_Signal_Round)
    print("The ROA # is ",T)

    HiPoint = np.where(dff_line == dff_line.max()) #find the peak response, which could be helpful for experiments with large signal wave

    Total_dff_line.append(dff_line)

    amplitude = []
    halfwidth = []
    AUC = []
    peaktime =[]
    signalonset = []
    decay = []
    rise = []
    signal_to_noise = []
    count = 0

    xseg_total = []
    def findpeaks(segment):
        max = segment[0]
        for i in segment:
            if i > max:
                max = i
        return max

    onset = 0

    for i in Total_Signal_Round:
        #plt.figure()
        a = i[0]
        b = i[1]
        if HiPoint[0] > a  and  HiPoint[0] < b:
            onset = a #onset of the largest signal in the trace
        x = dff_line[a:b+1]
        if b == framenumber:
            b = frameIndex   
        xseg = np.arange(a,b+1)
        xseg_total.append(xseg)
        peaks = findpeaks(x)
        peaksx = np.where(x == peaks)
        Amplitude_Signal = peaks
        half_width_s = 0.5 * Amplitude_Signal
        zHW = find_roots(xseg,x - half_width_s)
        AUCsignal = simpson(x,dx = 1)
        amplitude.append(Amplitude_Signal)
        signal_to_noise.append(Amplitude_Signal/max(baselineSumFiltered))
        AUC.append(AUCsignal)
        halfwidth.append((zHW[-1]-zHW[0])*frame_rate)
        tenpercent = 0.1 * Amplitude_Signal
        nintypercent = 0.9 * Amplitude_Signal
        ztenoercenpeak = find_roots(xseg,x-tenpercent)
        zninetypercent = find_roots(xseg, x-nintypercent)
        risetime = zninetypercent[0] -ztenoercenpeak[0]
        signaldecay = ztenoercenpeak[-1]-zninetypercent[-1]
        decay.append(signaldecay*frame_rate)
        rise.append(risetime*frame_rate)
        peaktimepoint = peaksx[0] + a
        signalonset.append(a*frame_rate)
        peaktime.append(peaktimepoint*frame_rate)



    AlignTrace = np.roll(dff_line,100-drug_frame)
    Total_AlignTrace.append(AlignTrace)

    Noise.append(Noise_trace)
    SignalNumber = len(Total_Signal)
    NumSigROA.append(SignalNumber)
    TraceHasSignal.append(Trace_number)
    TraceCount +=1

    duration = []
    for i in Total_Signal:
        duration.append(np.diff(i)[0]*frame_rate)

    
    TraceNumber = np.array(np.zeros(SignalNumber)+Trace_number).reshape(-1,SignalNumber)
    amplitude = np.array(amplitude).reshape(-1,SignalNumber)
    duration = np.array(duration).reshape(-1,SignalNumber)
    halfwidth = np.array(halfwidth).reshape(-1,SignalNumber)
    AUC = np.array(AUC).reshape(-1,SignalNumber)
    peaktime = np.array(peaktime).reshape(-1,SignalNumber)
    signalonset = np.array(signalonset).reshape(-1,SignalNumber)
    decay = np.array(decay).reshape(-1,SignalNumber)
    rise = np.array(rise).reshape(-1,SignalNumber)
    signal_to_noise = np.array(signal_to_noise).reshape(-1,SignalNumber)
    Data_signal = np.concatenate((TraceNumber.T, amplitude.T,duration.T,halfwidth.T,AUC.T,rise.T,decay.T,signalonset.T,peaktime.T,signal_to_noise.T),axis =1)
    #print (Data_signal)

    
    if len(FinalOutput) == 0:
        FinalOutput = Data_signal
    else:
        FinalOutput = np.concatenate((FinalOutput,Data_signal),axis = 0)

    '''This is a section that allows you to visualize the traces, comment the section if it's not needed 
        Red dots mark the intersection of baseline, green dots mark the intersection of threshold'''
    
    #print the label of all of the intersection
    # plt.figure(figsize=(10,5))
    # x = np.arange(framenumber)
    # plt.plot(x,dff_line,color = 'grey')
    # plt.plot(z1, np.zeros(len(z1))+Timesof_STDEV, marker="o", ls="", ms=2,color = 'green', label = 'threshold')
    # plt.plot(z2, np.zeros(len(z2))+ dff_baseline, marker="o", ls="", ms=2, color = 'red',label = 'baseline')
    # plt.legend(loc = 'upper left')
    # plt.xlabel("Time(s)")
    # plt.ylabel("Amplitude(dF/F)")
    # plt.show(block = False)


To_write = input("Do you want to write the results? please enter 'y' or 'n': ")
if To_write == 'y':
    NumSigROA = np.array(NumSigROA).reshape(-1, TraceCount)
    TraceHasSignal = np.array(TraceHasSignal).reshape(-1,TraceCount)
    Summary = np.concatenate((TraceHasSignal.T,NumSigROA.T), axis = 1)
    Total_dff_line = np.array(Total_dff_line)
    Total_AlignTrace = np.array(Total_AlignTrace)


    ##from this line below can be commented out depends on what kind of experiments analyzed 
    image_string_ROAs = open(ROAmask_path, 'rb').read()
    image_string_Cells = open(Cellmask_path, 'rb').read()
    map_ROAs = Image.open(io.BytesIO(image_string_ROAs))
    map_Cells = Image.open(io.BytesIO(image_string_Cells))
    ROAs = np.asarray(map_ROAs)
    Cells = np.asarray(map_Cells)
    print(ROAs)
    print(Cells)

    # label connected components in map1 and map2
    labeled_ROAs, num_labels_ROAs = label(ROAs,structure= [[1,1,1],[1,1,1],[1,1,1]])
    labeled_Cells, num_labels_Cells = label(Cells,structure = [[1,1,1],[1,1,1],[1,1,1]])
    print(num_labels_ROAs)
    print(num_labels_Cells)


    #initialize dictionary to store mapping between small and large regions 
    region_map = {}

    # #Initialize list to store assigned regions for each small region in ROAs
    assigned_regions = [[]for i in range(labeled_ROAs.max())]

    # #Iterate over each small region in ROAs
    for i in range(1,labeled_ROAs.max()+1):
        corresponding_regions = list(labeled_Cells[labeled_ROAs == i])
        if len(corresponding_regions) == 0:
            region_map[i] = np.nan
            assigned_regions[i-1] = []
        elif len(corresponding_regions) ==1:
            region_map[i] = corresponding_regions.pop()
            assigned_regions[i-1] = [region_map[i]]
        else:
            region_map[i] = max(corresponding_regions)

    #generate dataframe from the region map to get the list of numbers of ROA/cell
    df_ROA = pd.DataFrame.from_dict(region_map,orient = 'index').reset_index()
    df_ROA = df_ROA.rename(columns = {'index':'ROA number', 0: 'Cell number'})
    df_ROAcount = df_ROA.groupby(['Cell number']).agg(count = ('ROA number','count')).reset_index()

    #generate dataframe for the entire dataset
    df_final = pd.DataFrame(FinalOutput,columns = ['Trace Number','Amplitude','Duration','Halfwidth','Area Under Curve','Rise Time','Decay Time','Signal Onset','Peak Time','SNR'])
    df_final['Cell number'] = df_final['Trace Number'].map(region_map)
    
    # Para_10: Optional section for selecting the framees of interest for analysis.
    # df_final = df_final[(df_final['Signal Onset'].between (drug_frame-100, drug_frame )) | (df_final['Signal Onset']).between(drug_frame+100, drug_frame+200)] 

    df_final['Drug'] = np.where(df_final['Signal Onset']< drug_frame,'Before Drug',"After Drug")

    def Timelapse(Drugtype): #Para_11
        if Drugtype == 'Before Drug':
            return drug_Time # or according to the duration that is selected in Para_10
        elif Drugtype == 'After Drug':
            return (framenumber-drug_frame)*frame_rate # or according to the duration that is selected in Para_10

    # Count the number of signals in before and after drugs respectively
    df_ROAbased_count = df_final.groupby(['Trace Number','Drug']).size().reset_index(name = 'Signals in drug') 
    df_ROAbased = df_final.groupby(['Trace Number','Drug']).mean().reset_index()
    df_ROAbased['ROA signal type'] = df_ROAbased.groupby('Trace Number')['Drug'].transform(lambda x: 'stable' if ('Before Drug' in x.values) and ('After Drug' in x.values) else 'Off' if 'Before Drug' in x.values else 'On')
    df_ROAbased['Number of signals'] = df_ROAbased_count['Signals in drug']
    df_ROAbased['Time lapse'] = df_ROAbased['Drug'].apply(Timelapse)
    df_ROAbased['Frequency_min'] = df_ROAbased['Number of signals']/df_ROAbased['Time lapse']*60

    #generate dataframe for the grouped dataset according to the group and merge with the number of signals/cell and ROAs/cell
    df_SignalList = df_final.groupby(['Drug','Cell number']).agg(count = ('Cell number','count')).reset_index()
    df_final_mean = df_final.groupby(['Drug','Cell number']).mean()
    df_merged = pd.merge(df_SignalList,df_final_mean,on = ['Drug','Cell number'])
    df_merged['Trace Number'] = df_merged['Cell number'].map(df_ROAcount.set_index('Cell number')['count'])
    df_Cellbased = df_merged.rename(columns={'count':'# of signals','Trace Number':'# of ROAs'})
    df_Summary = pd.DataFrame(Summary,columns = ['Trace Number','# of signal per ROA'])

    # #choose either total dff_line or realigned trace, depends on what you need
    df_ROA_trace = pd.DataFrame(Total_dff_line)
    df_cell_trace = df_ROA_trace.reset_index(drop=True)
    df_cell_trace.index +=1
    df_cell_trace = df_cell_trace.reset_index()
    df_cell_trace = df_cell_trace.rename(columns = {'index':'ROA number'})
    df_cell_trace['Cell number'] = df_cell_trace['ROA number'].map(region_map)
    df_cell_trace = df_cell_trace.groupby(['Cell number']).mean()
    df_cell_trace = df_cell_trace.drop('ROA number',axis=1)

    # #choose either total dff_line or realigned trace, depends on what you need
    df_ROA_aligntrace = pd.DataFrame(Total_AlignTrace)
    df_cell_aligntrace = df_ROA_aligntrace.reset_index(drop=True)
    df_cell_aligntrace.index +=1
    df_cell_aligntrace = df_cell_aligntrace.reset_index()
    df_cell_aligntrace = df_cell_aligntrace.rename(columns = {'index':'ROA number'})
    df_cell_aligntrace['Cell number'] = df_cell_aligntrace['ROA number'].map(region_map)
    df_cell_aligntrace = df_cell_aligntrace.groupby(['Cell number']).mean()
    df_cell_aligntrace = df_cell_aligntrace.drop('ROA number',axis=1)

    #Get the intervals between signals for each trace 
    #Create a dataframe for the inhomogeneous list 
    data_list = IEI
    max_length = max(len(item) if isinstance(item, list) else 1 for item in data_list)
    data_list = [item if isinstance(item, list) else ['NaN'] * max_length for item in data_list]
    df_IEI = pd.DataFrame(data_list, columns=[f'col_{i}' for i in range(max_length)])
    df_Summary = pd.concat([df_Summary,df_IEI], axis = 1)
    flattened_IEI = []
    for sublist in data_list:
        flattened_IEI.extend(sublist)
    flattened_IEI = [x for x in flattened_IEI if x != 'NaN']
    df_flatten_IEI = pd.DataFrame(flattened_IEI,columns = ['Intervals'])

    with pd.ExcelWriter(output_path,engine = 'xlsxwriter') as writer:
        df_final.to_excel(writer, sheet_name=sheet1)
        df_Summary.to_excel(writer,sheet_name=sheet2)
        df_ROAbased.to_excel(writer,sheet_name=sheet3)
        df_Cellbased.to_excel(writer,sheet_name=sheet4)
        df_ROA_trace.to_excel(writer,sheet_name=sheet5) # All traces are included even those that don't have signals by default
        df_cell_trace.to_excel(writer,sheet_name =sheet6)# All traces are included even those that don't have signals by default
        df_ROA_aligntrace.to_excel(writer,sheet_name =sheet7) # All traces are included even those that don't have signals by default
        df_cell_aligntrace.to_excel(writer,sheet_name =sheet8)# All traces are included even those that don't have signals by default
        df_flatten_IEI.to_excel(writer, sheet_name= sheet9) # All the inter-event intervals are included

if To_write == 'n':
    exit()
