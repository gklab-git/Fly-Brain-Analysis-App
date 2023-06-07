## Brain Activity Analysis
"""This Qt-based application allows neuroscientists to semi-automatically load components of brain activation,
their respective time activity curves and source time-volume images. Simultaneous visualization of these set
of inputs along with easy navigation and application-specific algorithms for optimistically prioritizing a 
majority/all of the components that show relevant activity, allows them to find components of interest
effectively and efficiently and to save necessary information to corresponding output files.

The app allows instant saving of plots showing all/selected components mapped onto the brain alongside their 
time activity curves. In addition, it allows specification of time intervals where components show peaks,
calculation of those peaks, entry of additional information (such as brain region and user comments) and
instant saving of the collated information of component peaks and their corresponding time series 
respectively.  
* School of Life Sciences, Technical University of Munich
"""

## Import necessary packages
# Utils and general computation
from __future__ import unicode_literals
import sys
import os
import random
import io
import copy
import openpyxl
from PIL import Image
import numpy as np
import nibabel as nib
import scipy.io
from scipy import ndimage
import math
import pandas as pd

# Plotting
import matplotlib
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import matplotlib._color_data as mcd

# GUI
# Make sure that we are using QT5
matplotlib.use('Qt5Agg')
from PyQt5 import QtCore, QtWidgets
from scipy.signal import find_peaks
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtCore import *
from PyQt5.QtWidgets import QApplication, QTableView
from PyQt5.QtCore import QAbstractTableModel, Qt

## Initialize table that would store component entries and corresponding peak occurences 
df = pd.DataFrame({'Brain Region': ['Whole'],
                   'Component Number':[0],
                   'Start Time': [0],
                   'End Time': [500],
                   'Peak Value': [1],
                   'Comments': [None]})

dfr=df.copy()

## Adjunct GUI parameters
progname = os.path.basename(sys.argv[0])
progversion = "1.0"

#Segmentation algorithm
class RegionGrow2D:   
    def __init__(self, image):
        self.image = image
        self.x=image.shape[1]
        self.y=image.shape[0]
        self.outputMask = np.zeros_like(self.image)
        self.lowerThreshold = image.max()/3
        self.neighborMode = "8n"
        self.queue = []
    
    def main(self, seed):
        newItem = seed
        
        self.outputMask[newItem[0], newItem[1]] = 1
        self.queue.append((seed[0], seed[1]))
        
        while len(self.queue) != 0:
            newItem = self.queue.pop()
            if self.neighborMode == "8n":
                neighbors = [[newItem[0]-1, newItem[1]-1],
                             [newItem[0]-1, newItem[1]],
                             [newItem[0]-1, newItem[1]+1],
                             [newItem[0], newItem[1]-1],
                             [newItem[0], newItem[1]],
                             [newItem[0], newItem[1]+1],
                             [newItem[0]+1, newItem[1]-1],
                             [newItem[0]+1, newItem[1]],
                             [newItem[0]+1, newItem[1]+1]]
                for neighbor in neighbors:
                    self.checkNeighbour(neighbor[0], neighbor[1])
            elif self.neighborMode == "4n":
                self.checkNeighbour(newItem[0], newItem[1]-1)
                self.checkNeighbour(newItem[0], newItem[1]+1)
                self.checkNeighbour(newItem[0]-1, newItem[1])
                self.checkNeighbour(newItem[0]+1, newItem[1])
        return self.outputMask
        
    def checkNeighbour(self, y, x):
        if (x < self.x and y < self.y  
            and x > -1 and y > -1):
            intensity = self.image[y, x]
            if self.isIntensityAcceptable(intensity) and self.outputMask[y,x] == 0:
                self.outputMask[y,x] = 1
                self.queue.append((y, x))
    
    def isIntensityAcceptable(self, intensity):
        if intensity > self.lowerThreshold:
            return True
        return False 

def hex_to_rgb(value):
    value=value.lstrip('#')
    lv=len(value)
    return list(int(value[i:i+lv//3],16)/256 for i in range(0,lv,lv//3))

code_list=[]
i=0

for count,colour in enumerate(mcd.XKCD_COLORS):
    code=mcd.XKCD_COLORS[colour]
    rgb=hex_to_rgb(code)
    #if rgb[0]>0.95 or rgb[1]>0.95 or rgb[2]>0.95:
    #    continue
    code_list.append(rgb)
    if count>199:
        break

## Table model along with suitable functions/variables for adapting tabular data to a form suitable for usage in the Qt-based GUI
class pandasModel(QAbstractTableModel):

    def __init__(self, data):
        QAbstractTableModel.__init__(self)
        self._data = data

    def rowCount(self, parent=None):
        return self._data.shape[0]

    def columnCount(self, parent=None):
        return self._data.shape[1]

    def data(self, index, role=Qt.DisplayRole):
        if index.isValid():
            if role == Qt.DisplayRole:
                return str(self._data.iloc[index.row(), index.column()])
        return None

    def headerData(self, col, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return self._data.columns[col]
        return None

## Define base template class initializing the main plotting parameters and variables for subsequently visualizing data in the GUI
class MyMplCanvas(FigureCanvas):
    """Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.)."""

    def __init__(self, parent=None, width=7, height=7, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.canvas=FigureCanvas(self.fig)
        self.axes = self.fig.add_subplot(121)
        self.axes.axis('off')
        self.axes2 = self.fig.add_subplot(122)
        self.fig.tight_layout()

        #Parameter to control appearance of components overlaid on brain
        self.thr =0.2

        self.image_path=[]
        self.ica_path=[]
        self.image_data=[]
        self.ica_data=[]
        self.mat_path=[]
        self.D=[]
        self.fft=[]
        self.selective_order=[]

        self.compute_figure()

        FigureCanvas.__init__(self, self.fig)
        # self.setParent(parent)

        # # FigureCanvas.setSizePolicy(self,
        # #                            QtWidgets.QSizePolicy.Expanding,
        # #                            QtWidgets.QSizePolicy.Expanding)

        # FigureCanvas.updateGeometry(self)


## Define main class extending base template class for real-time visualization of imaging data and the corresponding time activity curves in the GUI
class MyStaticMplCanvas(MyMplCanvas):
    """Simple canvas with a sine plot."""
    def compute_figure(self, ica_path=[], value=0,n_proj=[]):        
        if len(ica_path)>0:
            self.ica_path=ica_path
            ica=nib.load(ica_path)
            self.ica_data=ica.get_fdata()
            for i in range(self.ica_data.shape[3]):
                for s in range(self.ica_data.shape[2]):           
                    self.ica_data[:,:,s,i]=ndimage.gaussian_filter(self.ica_data[:,:,s,i],2)
                    self.ica_data[:,:,s,i]/=self.ica_data[:,:,:,i].max()
            self.mat_path=ica_path.rstrip('IC.nii') + 'TS.mat'
            self.image_path=ica_path.rpartition('dFF')[0]+'.nii'
            self.axes.cla()
            self.axes2.cla()
            self.slider_value=0
            if os.path.isfile(self.image_path):
                image=nib.load(self.image_path)
                self.image_data=image.get_fdata()

        if n_proj!=[] and n_proj==0 and len(self.ica_path)>0:
            ica=nib.load(self.ica_path)
            self.ica_data=ica.get_fdata() 

            for i in range(self.ica_data.shape[3]):
                for s in range(self.ica_data.shape[2]):           
                    self.ica_data[:,:,s,i]=ndimage.gaussian_filter(self.ica_data[:,:,s,i],2)
                    self.ica_data[:,:,s,i]/=self.ica_data[:,:,:,i].max()

        if n_proj!=[] and n_proj!=0 and self.ica_data.shape[2]>=n_proj:
            print(self.ica_data.shape)
            ica=nib.load(self.ica_path)
            self.ica_data=ica.get_fdata()
            ica_eff=np.zeros((self.ica_data.shape[0],self.ica_data.shape[1],n_proj,self.ica_data.shape[3]))
            channels=int(self.ica_data.shape[2]/n_proj)
            k=self.ica_data.shape[2]%channels
            k=int(k/2)  
            for i in range(self.ica_data.shape[3]):
                for s in range(self.ica_data.shape[2]):           
                    self.ica_data[:,:,s,i]=ndimage.gaussian_filter(self.ica_data[:,:,s,i],2)
                    self.ica_data[:,:,s,i]/=self.ica_data[:,:,:,i].max()
                for s in range(n_proj):
                    ica_eff[:,:,s,i] = np.average(self.ica_data[:,:,k+channels*s:k+channels*(s+1),i],axis=2)
            self.ica_data=ica_eff       
         

        if (len(self.ica_path)>0):           
            if len(self.ica_data)>0:
                if len(self.image_data)==0:
                    self.image_data=np.zeros(np.shape(self.ica_data))
                    self.slice=np.ones((np.shape(self.ica_data)[0],np.shape(self.ica_data)[1],3))
                else:
                    self.slice=self.image_data[:,:,int(self.image_data.shape[2]/3),0].T
                    
                self.axes.imshow(self.slice)

                mid_point=int(np.shape(self.ica_data)[2]/2)
                ICA_mid=np.amax(self.ica_data[:,:,:,:],axis=2)
                ICA3D=np.squeeze(ICA_mid)
                se=np.ones((3,3))

                subslice=np.sqrt(self.slice)
                rg=RegionGrow2D(subslice)
                seed_pts=np.where(subslice==np.max(subslice))
                seed=[min(seed_pts[0]),min(seed_pts[1])]
                self.submask=rg.main(seed)

                if not os.path.isfile(self.mat_path):
                    mat={'TSo':np.zeros((1,ICA3D.shape[2]))}
                else:
                    mat=scipy.io.loadmat(self.mat_path)
                D=mat['TSo']
                self.D = D
                self.relevant_comps=[]
                self.remove_comps=[]

                for i in range(self.D.shape[1]):
                    Df=np.convolve(self.D[:,i],np.ones(40),'same')/40
                    peak_coords,_=find_peaks(Df, height=0)
                    peak_coords_min,_=find_peaks(-Df, height=0)
                    peak_coords_new=[p for p in peak_coords if Df[p]>0.5]+[p for p in peak_coords_min if -Df[p]>0.5]
                    if len(peak_coords_new)>0:
                        peak_coords_new=[p for p in peak_coords if abs(Df[p]-min(Df[max(0,p-50):p]))>0.5]
                    peak_coords_new=np.array(peak_coords_new)
                    if len(peak_coords_new)<7 and len(peak_coords_new)>2 and len(peak_coords)<20:
                        self.relevant_comps.append(i)
                n_comps=ICA3D.shape[2]

                C=np.zeros((3,n_comps))
                for i in range(n_comps):
                    if i<200:
                        C[:,i]=code_list[i]
                    elif i<n_comps:
                        C[:,i]=C[:,199]+(1-C[:,199])*(i-199)/(n_comps-199)

                self.C=C
                img_size=self.image_data.shape
                #value = uint8(app.BrainRegionSlider.Value)
                bcg=np.zeros((n_comps,img_size[1],img_size[0],4))

 
                for i in range(n_comps):
                    colour=i%12
                    ICAvali=ICA3D[:img_size[0],:img_size[1],i].T
                    ICAval=ICAvali/np.max(ICAvali)
                    ICAval[ICAval<0.4]=0
                    ICAval[ICAval>0.4]=1
                    ICAopenval=ndimage.binary_opening(ICAval,se)
                    #ICAopenval=ndimage.gaussian_filter(ICAopenval,1)
                    #ICAopenval=ICAopenval/ICAopenval.max()
                    corr=int(int(i/6)>0)
                    bcg[i,:,:,0]=C[0,i]*ICAopenval
                    bcg[i,:,:,1]=C[1,i]*ICAopenval
                    bcg[i,:,:,2]=C[2,i]*ICAopenval
                    bcg[i,:,:,3]=ICAopenval

                    # if (ICAopenval*self.submask).sum()<10*(ICAopenval*(1-self.submask)).sum():
                    #     self.remove_comps.append(i)
                    
                    #ICAval(ICAval(:,:,1)==1)=colour;
                    
                self.relevant_comps=[x for x in self.relevant_comps if x not in self.remove_comps]
                self.selective_order=np.zeros((np.shape(D)[1]))
                count=0
                for i in range(np.shape(D)[1]+len(self.relevant_comps)):
                    if (i<len(self.relevant_comps)):
                        self.selective_order[i]=self.relevant_comps[i]
                    else:
                        if (i-len(self.relevant_comps)) in self.relevant_comps:
                            count+=1
                        else:
                            self.selective_order[i-count]=i-len(self.relevant_comps)
                self.bcg = bcg                
                
                for i in range(n_comps):
                    self.axes.imshow(bcg[i,:,:,:]) 
                    self.axes.axis('off')                
                    self.axes2.plot(D[:,i]/np.max(abs(D[:,i]))+2*i,c=tuple(C[:,i].T))                          
                self.axes2.set_xlabel('Time(sec)')
                self.axes2.set_ylabel(r'$\Delta f/f$')
                self.fig.tight_layout()
                self.canvas.draw()
           

    def slider_moved(self, v):
        if len(self.ica_data)==0:
            return
        self.axes2.cla()
        self.slider_value=v
        if v==0:
            pass
        else:
            self.axes.imshow(self.slice)
            self.axes.imshow(self.bcg[v-1,:,:,:]) 
            self.axes.axis('off')                      
            self.axes2.plot(self.D[:,v-1],c=tuple(self.C[:,v-1].T))          
            self.axes2.set_xlabel('Time(sec)')
            self.axes2.set_ylabel(r'$\Delta f/f$')
            self.fig.tight_layout()
            self.canvas.draw()

    def comp_entered(self, v):
        if len(self.ica_data)==0:
            return
        self.axes2.cla()
        self.slider_value=v
        if v==0:
            pass
        else:         
            self.axes.imshow(self.slice)
            self.axes.imshow(self.bcg[v-1,:,:,:]) 
            self.axes.axis('off')                      
            self.axes2.plot(self.D[:,v-1],c=tuple(self.C[:,v-1].T))
            self.axes2.set_xlabel('Time(sec)')
            self.axes2.set_ylabel(r'$\Delta f/f$')
            self.fig.tight_layout()
            self.canvas.draw()
            
    def slider2_moved(self, v2):
        if len(self.ica_data)==0:
            return
        self.axes.cla()
        ICA_mid=self.ica_data[:,:,v2,:]
        ICA3D=np.squeeze(ICA_mid)
        se=np.ones((3,3))
        v=self.slider_value

        #colour=v%12
        ICAvali=ICA3D[:,:,v-1].T
        ICAval=ICAvali/np.max(ICAvali)
        ICAval[ICAval<self.thr]=0
        ICAval[ICAval>self.thr]=1
        ICAopenval=ndimage.binary_opening(ICAval,se)

        self.bcg[v-1,:,:,0]=self.C[0,v-1]*ICAopenval
        self.bcg[v-1,:,:,1]=self.C[1,v-1]*ICAopenval
        self.bcg[v-1,:,:,2]=self.C[2,v-1]*ICAopenval
        self.bcg[v-1,:,:,3]=ICAopenval


        if v==0:
            pass
        else:                               
            self.axes.imshow(self.slice)
            self.axes.imshow(self.bcg[v-1,:,:,:]) 
            self.axes.axis('off')
            self.fig.tight_layout()                      
            self.canvas.draw()

    def IC_entered(self, v2):
        if len(self.ica_data)==0:
            return
        self.axes.cla()
        if v2==-1:
            ICA_mid=np.max(self.ica_data,axis=2)
        else :
            ICA_mid=self.ica_data[:,:,v2,:]
        ICA3D=np.squeeze(ICA_mid)
        se=np.ones((3,3))
        v=self.slider_value

        #colour=v%12
        ICAvali=ICA3D[:,:,v-1].T
        ICAval=ICAvali/np.max(ICAvali)
        ICAval[ICAval<self.thr]=0
        ICAval[ICAval>self.thr]=1
        ICAopenval=ndimage.binary_opening(ICAval,se)

        self.bcg[v-1,:,:,0]=self.C[0,v-1]*ICAopenval
        self.bcg[v-1,:,:,1]=self.C[1,v-1]*ICAopenval
        self.bcg[v-1,:,:,2]=self.C[2,v-1]*ICAopenval
        self.bcg[v-1,:,:,3]=ICAopenval


        if v==0:
            pass
        else:                                  
            self.axes.imshow(self.slice)
            self.axes.imshow(self.bcg[v-1,:,:,:])   
            self.axes.axis('off') 
            self.fig.tight_layout()                   
            self.canvas.draw()

## Main application class based on PyQt5
class ApplicationWindow(QtWidgets.QMainWindow):
    def __init__(self):
        QtWidgets.QMainWindow.__init__(self)
        
        # Base attributes
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setWindowTitle("Application main window")
        
        # Menu bar
        self.menuBar = self.menuBar()
        self.menuBar.setNativeMenuBar(False)
        fileMenu=self.menuBar.addMenu('File')
        editMenu=self.menuBar.addMenu('Edit')
        
        # App window outline
        self.main_widget = QtWidgets.QWidget(self)
        l = QtWidgets.QVBoxLayout(self.main_widget)
        # Add plotting functionality
        self.sc = MyStaticMplCanvas(self.main_widget, width=5, height=4, dpi=100)
        l.addWidget(self.sc.canvas)
        self.main_widget.setFocus()
        self.setCentralWidget(self.main_widget)

        number_projection_menu = QMenu('Number of Projections',self)
        number_projection_menu.triggered.connect(self.select_projection)

        self.projection_action=[]

        #dict={'All (Default)':'0','One':'1','Three':'3','Five':'5','Seven':'7'}

        for item in ['All (Default)','One','Three','Five','Seven']:
            action=number_projection_menu.addAction(item)
            action.setCheckable(True)
            #action.setShortcut('Ctrl+'+dict[item])
            self.projection_action.append(action)
        
        # Add items to 'File' menu for performing load, save, selection and quit operations
        load_ica_action = QAction('Load ICA',self)
        load_ica_action.setShortcut('Ctrl+I')
        load_ica_action.triggered.connect(self.openICA)
        saveplots_action = QAction ('Save Plots',self)
        saveplots_action.setShortcut('Ctrl+P')
        saveplots_action.triggered.connect(self.savePlots)
        savetable_action = QAction ('Save Table',self)
        savetable_action.setShortcut('Ctrl+T')
        savetable_action.triggered.connect(self.saveTable)
        exit_action = QAction ('Exit',self)
        exit_action.setShortcut('Ctrl+Q')
        exit_action.triggered.connect(lambda:QApplication.quit())
        fileMenu.addMenu(number_projection_menu)
        fileMenu.addAction(load_ica_action)
        fileMenu.addAction(saveplots_action)
        fileMenu.addAction(savetable_action)
        fileMenu.addAction(exit_action)
        
        # Edit menu items
        undo_action = QAction('Undo',self)
        undo_action.setShortcut('Ctrl+Z')
        redo_action = QAction('Redo',self)
        redo_action.setShortcut('Ctrl+Y')
        cut_action = QAction ('Cut',self)
        cut_action.setShortcut('Ctrl+X') 
        
        # Flags to store table entries that are being deleted/restored for reverting operation, if needed
        self.undoFlag=0
        self.row=[]
        self.redoRow=[]
        self.storeRow=[]
        self.storeCol=[]
        
        # Add items to edit menu to perform cut, undo and redo operations
        editMenu.addAction(cut_action)
        cut_action.triggered.connect(self.cutEntry)
        editMenu.addAction(undo_action)
        editMenu.addAction(redo_action)
        undo_action.triggered.connect(self.undoEntry)
        redo_action.triggered.connect(self.redoEntry)
        

        #self.lcd=QLCDNumber()
        #self.lcd.setStyleSheet("QLCDNumber{background-color:gray;color:black;}")
        #self.lcd2=QLCDNumber()
        #self.lcd2.setStyleSheet("QLCDNumber{background-color:gray;color:black;}")
        # self.dial = QDial()
        # self.dial.setMinimum(0)
        # self.dial.setMaximum(43)
        # self.dial.setValue(0)
        self.ICGroupBox=QGroupBox(" IC and/or Cross Section")
        self.layout_12=QHBoxLayout()
        self.layout1=QVBoxLayout()
        self.layout2=QVBoxLayout()
        # Add slider and dialog box for allowing the user to navigate seamlessly through all components and instantly to desired components
        self.slider=QSlider(Qt.Horizontal)
        self.slider.valueChanged.connect(self.sliderMoved)
        
        # Form items to allow user to specify entries corresponding to peak occurence (namely - brain region, start and end time) and save them to table
        self.comp=QLineEdit()

        self.compBox=QDialogButtonBox(QDialogButtonBox.Ok)
        self.compBox.accepted.connect(self.compEntered)
        self.currentComp=0

        l.addWidget(self.ICGroupBox)
        l.addWidget(self.slider)
        l.addWidget(self.comp)
        l.addWidget(self.compBox)
        
        self.slider2=QSlider(Qt.Horizontal)
        self.slider2.setRange(1,36)
        self.slider2.valueChanged.connect(self.slider2Moved)
        self.ICslice=QLineEdit()
        self.ICsliceBox=QDialogButtonBox(QDialogButtonBox.Ok)
        self.ICsliceBox.accepted.connect(self.ICsliceEntered)
        self.currentICslice=0
        l.addWidget(self.slider2)
        l.addWidget(self.ICslice)
        l.addWidget(self.ICsliceBox)

        self.brainRegion=QLineEdit()
        self.startTime=QLineEdit()
        self.endTime=QLineEdit()
        self.comments=QLineEdit()
        self.buttonBox=QDialogButtonBox(QDialogButtonBox.Ok)
        self.buttonBox.accepted.connect(self.addEntry)
        
        # Create form and add respective items
        self.formGroupBox=QGroupBox(" Label and Process")
        layout=QFormLayout()
        layout.addRow(QLabel("Brain Region : "), self.brainRegion)
        layout.addRow(QLabel("Start Time : "), self.startTime)
        layout.addRow(QLabel("End Time : "), self.endTime)
        layout.addRow(QLabel("Comments : "), self.comments)
        layout.addRow(self.buttonBox)
        self.formGroupBox.setLayout(layout)
        l.addWidget(self.formGroupBox)
        
        # Add table to application window
        self.df=df
        model = pandasModel(self.df)
        self.view = QTableView()
        self.view.setModel(model)
        self.view.setColumnWidth(1,150)
        self.view.setColumnWidth(4,200)
        #self.view.selectionModel().currentChanged.connect(self.on_currentChanged)
        l.addWidget(self.view)
        
        # Variable to store time series
        self.ts=[]

        # Add information of last set of saved files to status bar
        self.statusLabel = QLabel()
        self.statusLabel.setText("")
        self.fileStatus = QLabel()
        self.fileStatus.setText("")
        self.statusBar=QStatusBar()
        self.statusBar.addWidget(self.statusLabel)
        self.statusBar.addPermanentWidget(self.fileStatus)
        l.addWidget(self.statusBar)
        # Variable to store brain region corresponding to selected IC
        self.regid=[]

        self.show()

    # Load IC file and the corresponding 4-D imaging data and time series, assign them to respective variables and display summary visualization and initial dummy table entry
    def openICA(self):
        icaPath, _ = QFileDialog.getOpenFileName()
        if icaPath =="":
            return
        self.statusLabel.setText(icaPath[icaPath.rindex('/')+1:])
        self.fileStatus.setText("No Files Saved")
        self.slider.setValue(0)
        self.slider2.setValue(0)
        self.ts=[]
        self.df=dfr.copy()
        model = pandasModel(self.df)
        self.view.setModel(model)
        # Calculate parameters and display all components mapped onto brain
        self.sc.compute_figure(ica_path=icaPath)
        self.slider2_lim=self.sc.ica_data.shape[2]
        self.slider2.setRange(1,self.slider2_lim)
        # Dummy table
        if len(self.sc.D)>0:
            self.slider.setRange(0,np.shape(self.sc.D)[1])
            self.df.loc[0,'Peak Value']=np.max(self.sc.D)
            #self.df.loc[-1,'Full Time Series']=self.sc.D
            model = pandasModel(self.df)
            self.view.setModel(model)
            self.regid=str(0)

    # Select number of projections
    def select_projection(self,action):
        ind=self.projection_action.index(action)
        for item in (self.projection_action[:ind]+self.projection_action[ind+1:]):
            item.setChecked(False)

        self.slider.setValue(0)
        self.slider2.setValue(0)
        self.ts=[]
        self.df=dfr.copy()
        model = pandasModel(self.df)
        self.view.setModel(model)

        key_value=[0,1,3,5,7]
        components=[self.slider2_lim-1,1,3,5,7]
        self.slider2.setRange(1,components[ind])
        self.sc.compute_figure(n_proj=key_value[ind])
        if len(self.sc.D)>0:
            self.slider.setRange(0,np.shape(self.sc.D)[1])
            self.df.loc[0,'Peak Value']=np.max(self.sc.D)
            #self.df.loc[-1,'Full Time Series']=self.sc.D
            model = pandasModel(self.df)
            self.view.setModel(model)
            self.regid=str(0)
            
    # Save plots containing either all the activated components highlighted over the brain or selected components along with the corresponding time series
    def savePlots(self):
        if len(self.sc.mat_path) == 0:
            return
        reg_id=self.slider.value()
        if reg_id in self.df['Component Number'].values:
            idx=self.df[self.df['Component Number']==reg_id].first_valid_index()
            reg_id=self.df.at[idx,'Brain Region']
        else:
            reg_id=str(reg_id)
        #fig=self.sc.fig.copy()
        self.sc.fig.savefig(os.getcwd()+'/'+self.sc.mat_path[self.sc.mat_path.rindex('/')+1:self.sc.mat_path.index('dFF')]+'_plots_'+reg_id+'.png')
        self.regid=reg_id
        if self.fileStatus.text()=="No Files Saved":
            self.fileStatus.setText("Last Plot Saved : "+self.regid)
        else:
            self.fileStatus.setText("Last Plot Saved : "+self.regid+", Table Saved")
    
    # Actively read slider value and display visualization of selected component and update text box
    def sliderMoved(self,event):
        if len(self.sc.selective_order)==0:
            self.slider.setValue(0)
            self.comp.setText(str(""))
            return

        v=int(self.sc.selective_order[self.slider.value()-1])+1
        self.comp.setText(str(v))
        self.sc.slider_moved(v)
            
    # Update visualization and slider position (and value) when input is given to text box and the 'Ok' button next to it is pressed
    def compEntered(self):
        if len(self.sc.selective_order)==0:
            self.comp.setText(str(""))
            return
        try:
            self.currentComp=int(self.comp.text())
            self.ind=0
            if self.currentComp!=0:
                self.ind=np.where(self.sc.selective_order==(self.currentComp-1))
                self.currentComp=int(self.sc.selective_order[self.ind])
            self.sc.comp_entered(v=self.currentComp)
            self.slider.setValue(int(self.ind[0])+1)
        except:
            self.currentComp=np.shape(self.sc.D)[1]
            self.sc.comp_entered(v=self.currentComp)
            self.slider.setValue(self.currentComp)

    # Actively read 2nd slider value and display visualization of selected volume slice of component and update text box (for 3D data)
    def slider2Moved(self,event):
            if len(self.sc.ica_data)==0:
                self.slider2.setValue(0)
                self.ICslice.setText(str(""))
                return
            self.sc.slider2_moved(v2=self.slider2.value()-1)
            self.ICslice.setText(str(self.slider2.value()))

    # Update IC slice visualization and 2nd slider position (and value) when input is given to text box and the 'Ok' button next to it is pressed (for 3D data)
    def ICsliceEntered(self):
        if len(self.sc.ica_data)==0:
            self.slider2.setValue(0)
            self.ICslice.setText(str(""))
            return
        try:
            self.currentICslice=int(self.ICslice.text())-1
            self.sc.IC_entered(v2=self.currentICslice)
            self.slider2.setValue(self.currentICslice+1)
        except:
            self.currentICslice=np.shape(self.sc.ica_data)[2]-1
            self.sc.IC_entered(v2=self.currentICslice)
            self.slider2.setValue(self.currentICslice+1)
            
    # Calculate peak values, adjust them to account for noise, and collect time series. As entries get added to the table, keep saving 
    # the output at regular intervals to ensure that the entire information is not accidentally lost
    def addEntry(self):
        self.storeRow=[]
        self.undoFlag=0
        t1=int(self.startTime.text())
        t2=int(self.endTime.text())
        if self.comments:
            comments=self.comments.text()
        else:
            comments=""
        i=int(self.sc.selective_order[self.slider.value()-1])+1
        ts=self.sc.D[t1:t2,i-1]
        peak_val=np.max(ts)
        min_val=np.min(ts)
        if abs(min_val)>abs(peak_val):
            peak_val=min_val
        if self.df.at[0,'Brain Region']=='Whole':
            self.df.loc[0]=[self.brainRegion.text(),i,t1,t2,peak_val,comments]
        else:
            self.df.loc[len(self.df)]=[self.brainRegion.text(),i,t1,t2,peak_val,comments]
        self.df.sort_values(['Brain Region'],ascending=True,inplace=True)
        self.df=self.df.reset_index(drop=True)
        idx=self.df[self.df['Component Number']==i].first_valid_index()
        #print(idx)
        reg_id=self.df.at[idx,'Brain Region']
        col=pd.DataFrame({reg_id:pd.Series(self.sc.D[:,i-1])})
        if len(self.ts)==0:
            self.ts=col
        else: 
            if reg_id not in self.ts.columns:
                self.ts=pd.concat([self.ts,col],axis=1)
            else:
                pass
        self.ts=self.ts.reindex(sorted(self.ts.columns),axis=1)
        model = pandasModel(self.df)
        self.view.setModel(model)
        
    # Save output tables with auto-assigned names that allow easy identification and organization
    def saveTable(self):
        if len(self.sc.mat_path) == 0:
            return
        self.df.to_excel(os.getcwd()+'/'+self.sc.mat_path[self.sc.mat_path.rindex('/')+1:self.sc.mat_path.index('.')]+'_region_peaks.xlsx',index=False)
        try:
            self.ts.to_excel(os.getcwd()+'/'+self.sc.mat_path[self.sc.mat_path.rindex('/')+1:self.sc.mat_path.index('.')]+'_complete_series.xlsx',index=False)
        except:
            pass
        if self.fileStatus.text()=="No Files Saved":
            self.fileStatus.setText("Table Saved")
        else:
            self.fileStatus.setText("Last Plot Saved : "+self.regid+", Table Saved")
            
    # Delete entries from table, store entries and set flags for undo/redo operations
    def cutEntry(self):
        self.row=self.view.currentIndex().row() 
        if len(self.df)==1:
            return
        self.undoFlag=1
        self.storeRow=self.df.loc[self.row].copy()
        self.undodf=self.df.copy() 
        
        reg=self.df.at[self.row,'Brain Region']
        self.df.drop([self.row],axis=0,inplace=True)
        if reg in self.df['Brain Region']:
            self.storeCol=[]
            return
        self.storeCol=self.ts[[reg]]
        self.undots=self.ts.copy() 

        self.ts.drop([reg],axis=1,inplace=True)   

        model = pandasModel(self.df)
        self.view.setModel(model)

    # Revert delete/restore operations on table (upto one entry)
    def undoEntry(self):
        if self.undoFlag!=1:
            return
        else:
            self.undoFlag=2   
            self.redodf=self.df.copy()         
            self.df=self.undodf.copy()
            if len(self.storeCol)==0:
                return
            self.redots=self.ts.copy()
            self.ts=self.undots.copy()
            model = pandasModel(self.df)
            self.view.setModel(model)

    # Re-execute delete/restore table options reverted using undo (upto one entry)
    def redoEntry(self):
        if self.undoFlag!=2:
            return
        else:
            self.undoFlag=1           
            self.df=self.redodf.copy()
            if len(self.storeCol)==0:
                return
            self.ts=self.redots.copy()
            model = pandasModel(self.df)
            self.view.setModel(model)
            
    # Exit operations and 'About' information
    def fileQuit(self):
        self.close()

    def closeEvent(self, ce):
        self.fileQuit()

    def about(self):
        QtWidgets.QMessageBox.about(self, "About",
                                    """Application for loading Fly brain data, plotting and saving"""
                                )

## Start the app using necessary application window parameters
qApp = QtWidgets.QApplication(sys.argv)

aw = ApplicationWindow()
aw.resize(700,900)
aw.setWindowTitle("Fly Brain Analysis")
aw.show()
sys.exit(qApp.exec_())
#qApp.exec_()