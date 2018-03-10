# Skeleton to Graph transformation using segments labeling technique.

The second way to create a graph from a 3D skeleton.

The method is based on the initial segmentation of a binary array into the regions: nodes and branches. This can be done using labeling function from the scikit-image library and the matrix of neighborhood, showing the number n of neighbors for each pixel. For the nodes n>3 and for the branches n=1,2. End nodes (n=1) are created during the connection process, when the edges are added to the graph. The connections can be established using binary dilation of the nodes. The output graph is a NetworkX object.
