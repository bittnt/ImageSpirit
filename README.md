# ImageSpirit Source code
This is the software bundle "ImageSpirit", created and maintained by:
  Shuia Zheng <szheng@robots.ox.ac.uk>
  Ming Ming Cheng <cmm.thu@gmail.com>


1. Dependencies
----------------------
This code is tested under Windows 8.1 system with Visual Studio 2010 and Visual Studio 2012. You are required to have OpenCV 2.4 and Qt 4.8 installed. For Qt 4.8, you may find BlueGo will be useful to help installing Qt4.8. The code may not work well with Qt 5.0+ at the moment.
>https://bitbucket.org/Vertexwahn/bluego



2. Build Process
----------------------
You open the sln file with visual studio, then compile it, if you set up all the paths.

3. Functions
----------------------
This software contains the implementation of using factorial DenseCRF in the application of semantic image segmentation. Further image editing application can be done with the CmmCode library.

4. Reference
----------------------
If you found the software is useful, please refer to:
```tex
@article{cheng2014imagespirit,
	author = {Cheng, Ming-Ming and Zheng, Shuai and Lin, Wen-Yan and Vineet, Vibhav and Sturgess, Paul and Crook, Nigel and Mitra, Niloy J. and Torr, Philip},
	title = {ImageSpirit: Verbal Guided Image Parsing},
	journal = {ACM Trans. Graph.},
	issue_date = {November 2014},
	volume = {34},
	number = {1},
	month = dec,
	year = {2014},
	issn = {0730-0301},
	pages = {3:1--3:11},
	articleno = {3},
	numpages = {11},
	doi = {10.1145/2682628},
	acmid = {2682628},
}
```
```tex
@inproceedings{DenseObjAtt_CVPR2014,
author = {Shuai Zheng and Ming-Ming Cheng and Jonathan Warrell and Paul Sturgess and Vibhav Vineet and Carsten Rother and Philip H. S. Torr},
title = {Dense Semantic Image Segmentation with Objects and Attributes},
booktitle = { IEEE International Conference on Computer Vision and Pattern Recognition (CVPR)},
address = {Columbus, Ohio, United States},
year = {2014}
}
```

5. License
----------------------
This software is free for non-commercial use. Any commercial use is strictly prohibited. 

This software has used DenseCRF software from Philipp Krähenbühl et al. and Permutohedral lattice software from Andrew Adams et al. Please also refer to their publications if this software is in use.
http://www.philkr.net/
http://people.csail.mit.edu/abadams/

Some part of the higher order potential implementation has also related to the implementation from Vibhav Vineet et al. Please refer to http://www.robots.ox.ac.uk/~tvg/publications.php

You can use the unary potential from ALE. 
Stable version: http://www.robots.ox.ac.uk/~phst/ale.htm
Latest version: http://www.inf.ethz.ch/personal/ladickyl/ALE.zip






