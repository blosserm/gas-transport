

# Example Useage
# -----
#python /Users/Matt/Documents/uscCode/gasTransport/guvTrack.py --flatfield  /Users/Matt/Documents/malmstadt\ lab/38msBackground512.tif /Users/Matt/Documents/malmstadt\ lab/GUVFLow_0.2ulmin_18_2.15.19.czi
#python /Users/Matt/Documents/uscCode/gasTransport/guvTrack.py /Users/Matt/Documents/malmstadt\ lab/GUVFLow_0.2ulmin_18_2.15.19.czi

import pims
import argparse
import bottleneck as bn
from skimage import io
import numpy as np
from tqdm import tqdm # adds a taskbar to any loop, e.g. for i in tqdm(range(10000)):
from matplotlib import pyplot as plt # for plotting
import os #for making directories
from scipy import signal
import csv
import skimage
from skimage import morphology, util, filters
from skimage.filters import threshold_otsu
from skimage.morphology import disk

# Get some variables from command line (nb switch raw_input for input in python v3 here it's v2.7)
parser = argparse.ArgumentParser()
parser.add_argument('filename', help='Name of input file, or the directory containing images')
parser.add_argument('--flatfield', help='subtract off the median value' )
parser.add_argument('--darkCurrent', help = 'dark current, make sure median is already dark current corrected')
args = parser.parse_args()
filename = args.filename
#constants
timeStep = 0.01 #in seconds
#channelPos1 = 252 #left channel left edge, in pixels
#channelPos2 = 32 #right channel left edge, in pixels
#channelPos3 = 146
channelPos1 = 0 #first channel left edge, in pixels
channelPos2 = 216 #second channel left edge, in pixels
channelPos3 = 419

#channelPos1 = 422 #left channel left edge, in pixels
#channelPos2 = 209 #right channel left edge, in pixels
#channelPos3 = 0

backgroundPos = 102 #empty comparison area to control for illumination fluctuations
channelW = 90 #channel width in pixels
medFilt = 501 #number of frames for rolling median filter, must be odd
vertPos = 100
vertH = 320
minGUVsize = 80 #in sqpixels
nBrightest = 5


print 'opening file'
# open the file using the magic of PIMS.
input_images = pims.Bioformats(filename)
#if neccessary, unwrap. Differs between tifs and czis
if len(input_images)==1:
    input_images = input_images[0]
stackL = len(input_images)
print stackL
if stackL==1:
    input_images = input_images[0]

print 'file open'
#input_images = np.float32(input_images)

if args.darkCurrent:
    print 'subtracting dark current'
    darkCurrent = np.float32(pims.Bioformats(args.darkCurrent))
    if len(darkCurrent)==1:
        darkCurrent = darkCurrent[0]
    #input_images = input_images-darkCurrent

if args.flatfield:
    print 'flat field correcting'
    backgroundFile = args.flatfield
    background = pims.Bioformats(backgroundFile)
    if len(background)==1:
        background = background[0]
    background = background/np.average(background)



#Find the median background, using 1000 frames over the course of the stack
print 'subtracting the median background'
subLen = 1000
iter = stackL/subLen
subStack = input_images[1:stackL:iter]
median = np.nanmedian(subStack, axis=0)

if args.darkCurrent:
    median = np.float32(median - darkCurrent)
if args.flatfield:
    median = median/background


#initialize variables to measure
int1 = np.float32(np.ones(stackL))
int2 = np.float32(np.ones(stackL))
int3 = np.float32(np.ones(stackL))
backInt = np.float32(np.ones(stackL))
times = np.uint32(range(stackL))*timeStep
radius1 = np.float32(np.zeros(stackL))
radius2 = np.float32(np.zeros(stackL))
radius3 = np.float32(np.zeros(stackL))


print 'measuring intensities'
for i in tqdm(range(stackL)):
#for i in tqdm([1,2,3,4,5,6,7,8,9,10]):
    #apply processing to the frame. doing so inside the loop is more memory efficient
    frame = input_images[i]
    frame = np.float32(frame)
    #print frame
    if args.darkCurrent:
        frame = np.float32(frame - darkCurrent)
    if args.flatfield:
        frame = frame/background
    frame = frame - median


    region1 = np.ndarray.flatten(frame[vertPos:vertPos+vertH, channelPos1:channelPos1+channelW])
    ind1 = np.argpartition(region1, -nBrightest)[-nBrightest:]
    int1[i] = np.mean(region1[ind1])
    region2 = np.ndarray.flatten(frame[vertPos:vertPos+vertH, channelPos2:channelPos2+channelW])
    ind2 = np.argpartition(region2, -nBrightest)[-nBrightest:]
    int2[i] = np.mean(region2[ind2])
    region3 = np.ndarray.flatten(frame[vertPos:vertPos+vertH, channelPos3:channelPos3+channelW])
    ind3 = np.argpartition(region3, -nBrightest)[-nBrightest:]
    int3[i] = np.mean(region3[ind3])
    #check to see if there's anything bright in this frame. If not, record the
    #background and exit this iteration. This saves a LOT of time.
    if np.max(frame) <400:
        continue
    frameToMask = filters.gaussian(frame, 3)
    #calculate the threshold
    thresh = threshold_otsu(frameToMask, nbins=2)
    #apply threshold to get a binary mask
    mask = frameToMask > thresh
    #clean up the mask by removing small regions, removing internal holes,
    #smoothing the edges, and dilating (get better agreement with subj radius)
    mask = morphology.remove_small_objects(mask, min_size=minGUVsize)
    mask = morphology.binary_closing(mask,disk(3))
    mask = morphology.convex_hull_image(mask)
    #mask = morphology.binary_dilation(mask)

    #label the different regions (GUVs) in the image
    label_image, number_of_labels = skimage.measure.label(mask, return_num=True)
    #extract measurements from the regions
    region = skimage.measure.regionprops(label_image, intensity_image=frame)

    #backInt[i] = sum(sum(frame[vertPos:vertPos+vertH, backgroundPos:backgroundPos+channelW]))
    #print sum(sum( backgroundPos:backgroundPos+channelW]))
    backInt = 1
    #find the intensity of the three channels, and normalize it to the control region
    #int1[i] = np.float32(np.divide(sum(sum(frame[vertPos:vertPos+vertH, channelPos1:channelPos1+channelW])), backInt[i]))
    #int2[i] = np.divide(sum(sum(frame[vertPos:vertPos+vertH, channelPos2:channelPos2+channelW])), backInt[i])
    #int3[i] = np.divide(sum(sum(frame[vertPos:vertPos+vertH, channelPos3:channelPos3+channelW])), backInt[i])
    #print frame[vertPos:vertPos+vertH, channelPos1:channelPos1+channelW]
    #int1[i] = np.max(np.float32(frame[vertPos:vertPos+vertH, channelPos1:channelPos1+channelW]))
    #int2[i] = np.max(frame[vertPos:vertPos+vertH, channelPos2:channelPos2+channelW])
    #int3[i] = np.max(frame[vertPos:vertPos+vertH, channelPos3:channelPos3+channelW])

    for j in range(len(region)):
        width = region[j].minor_axis_length
        bottom, left, top, right = region[j].bbox
        if bottom > 5 and top < 512-5 and region[j].major_axis_length < 256:
            if left>channelPos1 and right<(channelPos1+channelW):
                radius1[i] = width/2
                region1 = np.ndarray.flatten(region[j].intensity_image)
                ind1 = np.argpartition(region1, -nBrightest)[-nBrightest:]
                int1[i] = np.mean(region1[ind1])
            if left>channelPos2 and right<(channelPos2+channelW):
                radius2[i] = width/2
                region1 = np.ndarray.flatten(region[j].intensity_image)
                ind1 = np.argpartition(region1, -nBrightest)[-nBrightest:]
                int2[i] = np.mean(region1[ind1])
            if left>channelPos3 and right<(channelPos3+channelW):
                radius3[i] = width/2
                region1 = np.ndarray.flatten(region[j].intensity_image)
                ind1 = np.argpartition(region1, -nBrightest)[-nBrightest:]
                int3[i] = np.mean(region1[ind1])


    #region1 = np.ndarray.flatten(frame[vertPos:vertPos+vertH, channelPos1:channelPos1+channelW])
    #region1 = np.ndarray.flatten(region[0].intensity_image)
    #ind1 = np.argpartition(region1, -nBrightest)[-nBrightest:]
    #int1[i] = np.mean(region1[ind1])
    #region2 = np.ndarray.flatten(frame[vertPos:vertPos+vertH, channelPos2:channelPos2+channelW])
    #ind2 = np.argpartition(region2, -nBrightest)[-nBrightest:]
    #int2[i] = np.mean(region2[ind2])

filt1 = signal.medfilt(int1, medFilt)
filt2 = signal.medfilt(int2, medFilt)
filt3 = signal.medfilt(int3, medFilt)

peaks1 = int1 - filt1
peaks2 = int2 - filt2
peaks3 = int3 - filt3

#finds off set to line up traces of each channel
corr = signal.correlate(peaks2, peaks1, 'same')
maxCorr = np.where(corr == max(corr))
offset = maxCorr[0][0]-stackL/2

corr2 = signal.correlate(peaks3, peaks1, 'same')
maxCorr2 = np.where(corr2 == max(corr2))
offset2 = maxCorr2[0][0]-stackL/2

#identifies peaks


#janky transform to 2 pixel wide 'image' because i like the way skimage deals with regions
#intDub's are the intensities that will be measured. So need to decide if they're
#going to be the background sutracted or not.
isPeak1 = radius1 > 0
isPeak1 = np.c_[isPeak1,isPeak1]
isPeak1 = morphology.binary_closing(isPeak1,disk(4))
intDub1 = np.c_[peaks1, peaks1]
radDub1 = np.c_[radius1, radius1]
labelIm1, dumb = skimage.measure.label(isPeak1, return_num=True)
peakReg1 = skimage.measure.regionprops(labelIm1, intensity_image=intDub1)
radReg1 = skimage.measure.regionprops(labelIm1, intensity_image=radDub1)

isPeak2 = radius2 > 0
isPeak2 = np.c_[isPeak2,isPeak2]
isPeak2 = morphology.binary_closing(isPeak2,disk(4))
intDub2 = np.c_[peaks2, peaks2]
radDub2 = np.c_[radius2, radius2]
labelIm2, dumb = skimage.measure.label(isPeak2, return_num=True)
peakReg2 = skimage.measure.regionprops(labelIm2, intensity_image=intDub2)
radReg2 = skimage.measure.regionprops(labelIm2, intensity_image=radDub2)

isPeak3 = radius3 > 0
isPeak3 = np.c_[isPeak3,isPeak3]
isPeak3 = morphology.binary_closing(isPeak3,disk(4))
intDub3 = np.c_[peaks3, peaks3]
radDub3 = np.c_[radius3, radius3]
labelIm3, dumb = skimage.measure.label(isPeak3, return_num=True)
peakReg3 = skimage.measure.regionprops(labelIm3, intensity_image=intDub3)
radReg3 = skimage.measure.regionprops(labelIm3, intensity_image=radDub3)

filenameprefix = os.path.splitext(filename)[0]
if not os.path.isdir(filenameprefix):
    os.mkdir(filenameprefix)

with open(filenameprefix+'/' +'peaks.csv', 'wb') as peakFile:
    wr = csv.writer(peakFile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    wr.writerow(['position', 'event num', 'frame', 'pred. frame', 'mean intensity', 'intensity stdev', 'radius', 'radius stdev'])
    for i in range(len(peakReg1)):
        wr.writerow(['1', i, peakReg1[i].centroid[0],peakReg1[i].centroid[0], peakReg1[i].mean_intensity, np.std(peakReg1[i].intensity_image),
            radReg1[i].mean_intensity, np.std(radReg1[i].intensity_image)])
    for i in range(len(peakReg2)):
        wr.writerow(['2', i, peakReg2[i].centroid[0], peakReg1[i].centroid[0]-offset, peakReg2[i].mean_intensity, np.std(peakReg2[i].intensity_image),
            radReg2[i].mean_intensity, np.std(radReg2[i].intensity_image)])
    for i in range(len(peakReg3)):
        wr.writerow(['3', i, peakReg3[i].centroid[0],peakReg3[i].centroid[0]-offset2, peakReg3[i].mean_intensity, np.std(peakReg3[i].intensity_image),
            radReg3[i].mean_intensity, np.std(radReg3[i].intensity_image)])

peakLoc1 = np.ones(len(peakReg1))
peakLoc2 = np.ones(len(peakReg2))
peakLoc3 = np.ones(len(peakReg3))
for i in range(len(peakReg1)):
    peakLoc1[i] = int(peakReg1[i].centroid[0])
for i in range(len(peakReg2)):
    peakLoc2[i] = int(peakReg2[i].centroid[0])
for i in range(len(peakReg3)):
    peakLoc3[i] = int(peakReg3[i].centroid[0])

#separate and process movies of just guvs
peakLoc1 = peakLoc1.astype(int)
peakLoc2 = peakLoc2.astype(int)
peakLoc3 = peakLoc3.astype(int)
if len(peakLoc1)>0:
    peakIm1 = input_images[peakLoc1]
    if args.darkCurrent:
        peakIm1 = np.float32(peakIm1 - darkCurrent)
    if args.flatfield:
        peakIm1 = peakIm1/background
    peakIm1 = peakIm1 - median
    for i in range(len(peakIm1)):
        skimage.io.imsave(filenameprefix+'/' +'guvPos1.tif', peakIm1[i], append=True)

if len(peakLoc2)>0:
    peakIm2 = input_images[peakLoc2]
    if args.darkCurrent:
        peakIm2 = np.float32(peakIm2 - darkCurrent)
    if args.flatfield:
        peakIm2 = peakIm2/background
    peakIm2 = peakIm2 - median
    for i in range(len(peakIm2)):
        skimage.io.imsave(filenameprefix+'/' +'guvPos2.tif', peakIm2[i], append=True)

if len(peakLoc3)>0:
    peakIm3 = input_images[peakLoc3]
    if args.darkCurrent:
        peakIm3 = np.float32(peakIm3 - darkCurrent)
    if args.flatfield:
        peakIm3 = peakIm3/background
    peakIm3 = peakIm3 - median
    for i in range(len(peakIm3)):
        skimage.io.imsave(filenameprefix+'/' +'guvPos3.tif', peakIm3[i], append=True)

#save images





with open(filenameprefix+'/' +'intensity.csv', 'wb') as csvfile:
    wr = csv.writer(csvfile, quoting=csv.QUOTE_ALL)
    wr.writerows(zip(*[['frame number'], ['time (s)'], ['raw Int 1'], ['raw Int2'],
    ['raw Int 3'], ['median filter1'], ['median filter 2'], ['median filter 3'],
    ['filtered Int 1'], ['filtered Int 2'], ['filtered Int 3'],
    ['radius 1'], ['radius 2'], ['radius 3']]))
    wr.writerows(zip(*[times/timeStep, times, int1, int2, int3, filt1, filt2, filt3,
    peaks1, peaks2, peaks3, radius1, radius2, radius3]))





####plots####
#convert times back to frame number


fig1 = plt.figure(1)
plt.plot(times, int1)
plt.plot(times, int2)
plt.plot(times, int3)
plt.xlabel('t (s)')
plt.ylabel('Intensity (arb.)')
plt.savefig(filenameprefix+'/' +'rawI.pdf', bbox_inches='tight')

fig2 = plt.figure(2)
plt.plot(times, filt1)
plt.plot(times, filt2)
plt.plot(times, filt3)
plt.xlabel('t (s)')
plt.ylabel('Intensity (arb.)')
plt.savefig(filenameprefix+'/' +'filter.pdf', bbox_inches='tight')

fig3 = plt.figure(3)
ax = fig3.add_subplot(111)
plt.plot(times, peaks1)
plt.plot(times, peaks2)
plt.plot(times, peaks3)
plt.xlabel('t (s)')
plt.ylabel('Intensity (arb.)')
peakLoc1 = peakLoc1.astype(int)
for i in peakLoc1:
    ax.annotate(i, (i, peaks1[i]))
plt.savefig(filenameprefix+'/' +'filteredI.pdf', bbox_inches='tight')

fig4 = plt.figure(4)
plt.plot(times, peaks1+max(peaks2))
plt.plot(times - offset*timeStep, peaks2)
plt.plot(times - offset2*timeStep, peaks3-max(peaks3)*1.1)
plt.xlabel('t (s)')
plt.ylabel('Intensity (arb.)')
plt.savefig(filenameprefix+'/' +'offsetIinTandI.pdf', bbox_inches='tight')

fig5 = plt.figure(5)
plt.plot(times, peaks1)
plt.plot(times - offset*timeStep, peaks2)
plt.plot(times - offset2*timeStep, peaks3)
plt.xlabel('t (s)')
plt.ylabel('Intensity (arb.)')
plt.savefig(filenameprefix+'/' +'offsetIinT.pdf', bbox_inches='tight')

fig6 = plt.figure(8)
plt.plot(times, radius1)
plt.plot(times, radius2)
plt.plot(times, radius3)
plt.xlabel('frame')
plt.ylabel('Radius (pixels)')
plt.savefig(filenameprefix+'/' +'radius.pdf', bbox_inches='tight')





with open(filenameprefix+'/' +'params.txt', 'wb') as txt:
    txt.write("file name = " + filename)
    if args.flatfield:
        txt.write("\nflatfield = " + backgroundFile)
    else:
        txt.write("\nflatfield = False")
    if args.darkCurrent:
        txt.write("\ndark current = "+ args.darkCurrent)
    else:
        txt.write("\ndark current = False")
    txt.write("\nchannel 1 left edge = " + str(channelPos1))
    txt.write("\nchannel 2 left edge = " + str(channelPos2))
    txt.write("\nchannel 3 left edge = " + str(channelPos3))
    txt.write("\nbackground left edge = " + str(backgroundPos))
    txt.write("\nchannel width = " + str(channelW))
    txt.write("\nroi top edge = " + str(vertPos))
    txt.write("\nroi height = " + str(vertH))
    txt.write("\ntime step = " + str(timeStep))
    txt.write("\nmedian filter size = " + str(medFilt))
    txt.write("\nestimated lag time (frames) = " + str(offset))
    txt.write("\nestimated lag time (timses) = " + str(offset*timeStep))
    txt.write("\nestimated lag time 2 (frames) = " + str(offset2))
    txt.write("\nestimated lag time (timses) = " + str(offset2*timeStep))



#plt.show()



        #io.imsave(filename + '.output.tif',im_out)

#    if args.normalise :
#        overall_median = bn.nanmedian(im_out, axis=0)
#        im_out = np.float32(im_out / overall_median)
        #print "overall_median", overall_median
