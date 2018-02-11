# entropicOfCancer
k-nearest neighbor method of finding early signs of cervical cancer from anonymous bioimaging data

Written in C, this algorithm breaks down .jpegs into manageable chunks, derives some data from the pixel color data,
and then uses given guidelines to delineate between stage 1, stage 2, and stage 3 cancer. 

Essentially this is a classifier, and while I recognized the potential for using deep learning, after the effectiveness I saw
in my work at the CCRC in clustering algorithms I decided that clustering the pixels into classes based on certain RGB thresholds would
be the most effective way to locate cancerous tissue, as there is a clear threshold after which the tissue is recognizably cancerous. 

The location of the 'reddish' clusters is also important when determining the stage of the cancer, so I had to build a 2-D array from the .jpeg's 
pixel data so I could understand where in the image each cluster of similarly hued pixels are. I ended up using various terminology that was consistent
with the biological nomenclature (lobe, node) to help translate the data into an organically driven analog.

Also of use was the Hamming distance between nodes of pixels that I used to recognize larger clusters that were beyond the scope of the smallish nodes.
By this method I was able to create lobes of data that were either recognizable cancerous or not. This was the way I translated the k-nearest neighbor technique
into my program, because if there is a nearly zero Hamming distance from one node to another I know that those nodes are within a larger node. I'm sure there might be more robust ways of doing this
but I wanted to incorporate k-NN into my program and so I just shoehorned it via the Hamming distance. Once the program is finished I'll be able to better understand
how I migth have done it differently, and whether or not it was useful to use the Hamming distane at all. 

Note that this program does not yet compile and is a work in progress. I tinker with it every once in a while to hepl keep my C chops up to snuff
and to stay fresh in the bioinformatic domain. 
