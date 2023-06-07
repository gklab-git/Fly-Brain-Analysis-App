# Fly Brain Analysis GUI

This GUI is a user-friendly platform for analyzing whole brain activity patterns alongside their corresponding time series and saving the respective information.

> This project started off as a fun idea to be able to look at whole brain imaging data in a way that was engaging and interactive. Thus, emerged the idea of a GUI.


## Getting Started

The app allows instant saving of images of all/selected independent components of activity mapped onto the brain alongside their  time activity curves. In addition, it allows specification of time intervals where components show peaks,
calculation of those peaks, entry of additional information (such as brain region and user comments) and instant saving of the collated information of component peaks and their corresponding time series 
respectively.

### Prerequisites
The GUI is based entirely in Python 3 and includes the following dependencies
> PIL >= 8.3.1
> 
> matplotlib >= 3.4.3
> 
> numpy >= 1.19.5
> 
> nibabel >= 3.2.1
> 
> scipy >= 1.7.1
> 
> pandas >= 1.3.2
> 
> xlwt
> 
> pyqt5
> 

## Usage
1. Select the main IC image containing components to be visualized. The GUI would automatically load corresponding time series files and a cross sectional view of the brain underlying the independent components (ICs)

![python_ZVNiaW8Wku](https://github.com/gklab-git/fly-brain-analysis-app/assets/118174779/93c3b589-8300-4d39-b011-2e19776064a5)

2. Navigate to a component of choice using the slider (or alternatively by entering the component number)

![python_nE8zSte1OH](https://github.com/gklab-git/fly-brain-analysis-app/assets/118174779/26d2020d-55be-48ed-80eb-67e1261287bd)

3. Visualization of the corresponding brain activation and time series

![python_R0hLzGEucR](https://github.com/gklab-git/fly-brain-analysis-app/assets/118174779/9926efcb-9af9-48ee-88c5-0594cd31b22f)

4. Entry of peak information into table and saving respective plots and table data

![python_TcgiPyYItk](https://github.com/gklab-git/fly-brain-analysis-app/assets/118174779/da1e2e6b-cd57-409e-9303-fcb73246c554)


## Contributing

The basic framework of the app would allow adaptation to other domains of neuroscience and neuroimaging and could be a simple and intuitive platform to navigate through these activity patterns and save information.

## Acknowledgments

The GUI was developed by S. Parhi as part of Grunwald Kadow Lab under the mentorship of K. P. Siju and in collaboration with P. Bandow.

