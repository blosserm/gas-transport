# gas-transport
code for analyzing co2 transport experiments on guvs, used in the paper https://www.biorxiv.org/content/10.1101/2020.11.16.384958v1.abstract

guvTrack extracts intensity values from a tiff stack (czi is the default), expecting a moving guv flowing through a serpentine channel
A intensity trace is generated, and discrete events corresponding to individual guv traversals are created.

gasTransportFitFunc reads the output of guvTrack, combines them, and fits an exponential decay


