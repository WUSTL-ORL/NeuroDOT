NEURODOT README

Table of Contents:

1. Installation
2. Dependencies
3. Getting Started
4. Reference Material
5. Feedback

1. Installation

Note: NeuroDOT is optimized for use in MATLAB 2020b on Windows and Unix. If you do not have MATLAB, contact your IT department or go to www.mathworks.com for more information on MATLAB. Your results may vary with other versions or operating systems.

To install NeuroDOT, copy the unzipped folder into the directory of your choosing. Then open MATLAB and either click on the "Set Path" button on the "Home" tab of the main MATLAB console, or copy and paste the following command into your command line, substituting the bracketed words for the path into which you unzipped NeuroDOT:

addpath(genpath('[INSTALL DRIVE AND PATH]\NeuroDOT\'))

2. Dependencies

To use NeuroDOT, the following dependencies are necessary to be installed:
	1. Matlab 2020b (https://www.mathworks.com/products/new_products/release2020b.html)
		1. Matlab 2020b and add-on toolboxes:
		2. Signal Processing Toolbox
		3. Deep Learning Toolbox
		4. Image Processing Toolbox
		5. Statistics and Machine Learning Toolbox
		6. Parallel Computing Toolbox
   	2. NIRFASTer (https://github.com/nirfaster/NIRFASTer) 
   	3. SNIRF (https://github.com/fNIRS/snirf) 
  	4. easyh5 (https://github.com/NeuroJSON/easyh5) 
  	5. jsnirfy (https://github.com/NeuroJSON/jsnirfy) 
   	6. GIFTI (https://github.com/gllmflndn/gifti) 
   	7. FreeSurfer 7.2 (https://surfer.nmr.mgh.harvard.edu/fswiki/rel7downloads) 
   	8. Connectome Workbench (https://humanconnectome.org/software/get-connectome-workbench) 

3. Getting Started

The toolbox contains 4 folders: Data, Documentation, Functions, and Support_Files. In the Documentation folder you will find the User Manual, the Overview, and the  Tutorials. These will give you all the basic information needed to work with fNIRS data in NeuroDOT. Additionally, there are scripts in the Documentation/Scripts folder that cover all of the material in the various Tutorials.


4. Reference Material

Also within the "Documentation" folder are several other tutorial appendices for beginner and advanced users, function reference files, and extended versions of each tutorial that can be viewed in the MATLAB Help Viewer. Additionally, for any function in NeuroDOT, you can type "help [function name]" into the MATLAB command line and view the function reference that way.


5. Feedback

If you have any suggestions or questions regarding NeuroDOT, please submit your feedback through this Google Form: https://forms.gle/iEYfEZhfj99FVEs29

NeuroDOT is constantly undergoing new developments to adapt to the evolving needs of the fNIRS community. Please do not hesitate to reach out if you have questions: neurodot-support@wustl.edu
