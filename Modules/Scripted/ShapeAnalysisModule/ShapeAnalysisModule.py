import os, sys
import unittest
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import logging
import csv
from slicer.util import VTKObservationMixin
import platform
import time

#
# ShapeAnalysisModule
#

class ShapeAnalysisModule(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "Shape Analysis Module"
    self.parent.categories = ["Examples"]
    self.parent.dependencies = []
    self.parent.contributors = ["Laura Pascal (Kitware Inc.)"]
    self.parent.helpText = """
    """
    self.parent.acknowledgementText = """
    This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc.
    and Steve Pieper, Isomics, Inc. and was partially funded by NIH grant 3P41RR013218-12S1.
    """ # replace with organization, grant and thanks.

#
# ShapeAnalysisModuleWidget
#

class ShapeAnalysisModuleWidget(ScriptedLoadableModuleWidget):
  """Uses ScriptedLoadableModuleWidget base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setup(self):
    ScriptedLoadableModuleWidget.setup(self)

    #
    #   Global variables
    #
    self.Logic = ShapeAnalysisModuleLogic(self)

    #
    #  Interface
    #
    loader = qt.QUiLoader()
    self.moduleName = 'ShapeAnalysisModule'
    scriptedModulesPath = eval('slicer.modules.%s.path' % self.moduleName.lower())
    scriptedModulesPath = os.path.dirname(scriptedModulesPath)
    path = os.path.join(scriptedModulesPath, 'Resources', 'UI', '%s.ui' % self.moduleName)
    qfile = qt.QFile(path)
    qfile.open(qt.QFile.ReadOnly)
    widget = loader.load(qfile, self.parent)
    self.layout = self.parent.layout()
    self.widget = widget
    self.layout.addWidget(widget)

    # Global variables of the Interface
    #   Group Project IO
    self.GroupProjectInputDirectory = self.getWidget('DirectoryButton_GroupProjectInputDirectory')
    self.GroupProjectOutputDirectory = self.getWidget('DirectoryButton_GroupProjectOutputDirectory')
    self.Debug = self.getWidget('checkBox_Debug')
    #   Post Processed Segmentation
    self.OverwriteSegPostProcess = self.getWidget('checkBox_OverwriteSegPostProcess')
    self.RescaleSegPostProcess = self.getWidget('checkBox_RescaleSegPostProcess')
    self.sx = self.getWidget('SliderWidget_sx')
    self.sy = self.getWidget('SliderWidget_sy')
    self.sz = self.getWidget('SliderWidget_sz')
    self.label_sx = self.getWidget('label_sx')
    self.label_sy = self.getWidget('label_sy')
    self.label_sz = self.getWidget('label_sz')
    self.LabelState = self.getWidget('checkBox_LabelState')
    self.label_ValueLabelNumber = self.getWidget('label_ValueLabelNumber')
    self.ValueLabelNumber = self.getWidget('SliderWidget_ValueLabelNumber')
    #   Generate Mesh Parameters
    self.OverwriteGenParaMesh = self.getWidget('checkBox_OverwriteGenParaMesh')
    self.NumberofIterations = self.getWidget('SliderWidget_NumberofIterations')
    #   Parameters to SPHARM Mesh
    self.OverwriteParaToSPHARMMesh = self.getWidget('checkBox_OverwriteParaToSPHARMMesh')
    self.SubdivLevelValue = self.getWidget('SliderWidget_SubdivLevelValue')
    self.SPHARMDegreeValue = self.getWidget('SliderWidget_SPHARMDegreeValue')
    self.thetaIterationValue = self.getWidget('spinBox_thetaIterationValue')
    self.phiIterationValue = self.getWidget('spinBox_phiIterationValue')
    self.medialMesh = self.getWidget('checkBox_medialMesh')
    #   Advanced Post Processed Segmentation
    self.GaussianFiltering = self.getWidget('checkBox_GaussianFiltering')
    self.label_VarianceX = self.getWidget('label_VarianceX')
    self.VarianceX = self.getWidget('SliderWidget_VarianceX')
    self.label_VarianceY = self.getWidget('label_VarianceY')
    self.VarianceY = self.getWidget('SliderWidget_VarianceY')
    self.label_VarianceZ = self.getWidget('label_VarianceZ')
    self.VarianceZ = self.getWidget('SliderWidget_VarianceZ')
    #   Advanced Parameters to SPHARM Mesh
    self.useRegTemplate = self.getWidget('checkBox_useRegTemplate')
    self.regTemplate = self.getWidget('PathLineEdit_regTemplate')
    self.useFlipTemplate = self.getWidget('checkBox_useFlipTemplate')
    self.flipTemplate = self.getWidget('PathLineEdit_flipTemplate')
    self.MTemplate = self.getWidget('checkBox_MTemplate')
    self.ParaOut1Template = self.getWidget('checkBox_ParaOut1Template')
    self.choiceOfFlip = self.getWidget('comboBox_choiceOfFlip')
    #   Flip Options
    self.tableWidget_ChoiceOfFlip = self.getWidget('tableWidget_ChoiceOfFlip')
    self.changeFlips = self.getWidget('pushButton_changeFlips')
    self.visualizationOfFlipInSPV = self.getWidget('pushButton_visualizationOfFlipInSPV')
    #   Visualization
    #   Apply CLIs
    self.ApplyButton = self.getWidget('applyButton')
    self.progress_layout = self.getWidget('progress_layout')

    # Connections
    #   Group Project IO
    #   Post Processed Segmentation
    self.RescaleSegPostProcess.connect('clicked(bool)', self.onSelectSpacing)
    self.LabelState.connect('clicked(bool)', self.onSelectValueLabelNumber)
    #   Generate Mesh Parameters
    #   Parameters to SPHARM Mesh
    #   Advanced Post Processed Segmentation
    self.GaussianFiltering.connect('clicked(bool)', self.onSelectGaussianVariance)
    #   Advanced Parameters to SPHARM Mesh
    #   Flip Options
    self.changeFlips.connect('clicked(bool)', self.onChangeFlips)
    self.visualizationOfFlipInSPV.connect('clicked(bool)', self.onPreviewFlips)
    #   Visualization
    #   Apply CLIs
    self.ApplyButton.connect('clicked(bool)', self.onApplyButton)

    # Widget Configuration
    #     Table for the Flip Options
    self.tableWidget_ChoiceOfFlip.setColumnCount(2)
    self.tableWidget_ChoiceOfFlip.setHorizontalHeaderLabels([' Files ', ' Choice of Flip '])
    self.tableWidget_ChoiceOfFlip.setColumnWidth(0, 400)
    horizontalHeader = self.tableWidget_ChoiceOfFlip.horizontalHeader()
    horizontalHeader.setStretchLastSection(False)
    horizontalHeader.setResizeMode(0, qt.QHeaderView.Stretch)
    horizontalHeader.setResizeMode(1, qt.QHeaderView.ResizeToContents)
    self.tableWidget_ChoiceOfFlip.verticalHeader().setVisible(False)

    #     Progress Bar
    self.progress_layout.addWidget(self.Logic.ProgressBar)
    self.CLIProgressBars = qt.QWidget()
    self.progress_layout.addWidget(self.CLIProgressBars)

  def cleanup(self):
    pass

  # Functions to recovery the widget in the .ui file
  def getWidget(self, objectName):
    return self.findWidget(self.widget, objectName)

  def findWidget(self, widget, objectName):
    if widget.objectName == objectName:
      return widget
    else:
      for w in widget.children():
        resulting_widget = self.findWidget(w, objectName)
        if resulting_widget:
          return resulting_widget
    return None

  #
  #   Post Processed Segmentation
  #
  def onSelectSpacing(self):
    self.label_sx.enabled = self.RescaleSegPostProcess.checkState()
    self.label_sy.enabled = self.RescaleSegPostProcess.checkState()
    self.label_sz.enabled = self.RescaleSegPostProcess.checkState()
    self.sx.enabled = self.RescaleSegPostProcess.checkState()
    self.sy.enabled = self.RescaleSegPostProcess.checkState()
    self.sz.enabled = self.RescaleSegPostProcess.checkState()

  def onSelectValueLabelNumber(self):
    self.label_ValueLabelNumber.enabled = self.LabelState.checkState()
    self.ValueLabelNumber.enabled = self.LabelState.checkState()

  #
  #   Advanced Post Processed Segmentation
  #
  def onSelectGaussianVariance(self):
    self.label_VarianceX.enabled = self.GaussianFiltering.checkState()
    self.VarianceX.enabled = self.GaussianFiltering.checkState()
    self.label_VarianceY.enabled= self.GaussianFiltering.checkState()
    self.VarianceY.enabled = self.GaussianFiltering.checkState()
    self.label_VarianceZ.enabled = self.GaussianFiltering.checkState()
    self.VarianceZ.enabled = self.GaussianFiltering.checkState()

  #
  #   Apply CLIs
  #
  def onApplyButton(self):
    # Run workflow
    if not self.Logic.Node.IsBusy():
      logging.info('Widget: Running ShapeAnalysisModule')
      self.ApplyButton.setText("Cancel")

      self.Logic.addObserver(self.Logic.Node, slicer.vtkMRMLCommandLineModuleNode().StatusModifiedEvent,
                             self.onLogicModified)

      self.Logic.Node.SetStatus(self.Logic.Node.Scheduled)
      self.Logic.allCaseStartTime = time.time()

      self.Logic.clearFlipOptionsTable()

      self.Logic.ShapeAnalysisCases()

    # Cancel Workflow
    else:
      logging.info("Widget: Cancelling ShapeAnalysisModule")
      self.ApplyButton.setEnabled(False)
      self.Logic.Cancel()

  def onLogicModified(self, logic_node, event):
    status = logic_node.GetStatusString()
    logging.info('-- %s : ShapeAnalysisModule', status)

    # if not busy (completed, error, cancelled)
    if not logic_node.IsBusy():
      self.Logic.removeObserver(logic_node, slicer.vtkMRMLCommandLineModuleNode().StatusModifiedEvent,
                                self.onLogicModified)
      # Create Error Message
      if status == 'Completed with errors':
        logging.error(self.Logic.ErrorMessage)
        qt.QMessageBox.critical(slicer.util.mainWindow(),
                                'ShapeAnalysisModule',
                                self.Logic.ErrorMessage)

      #  Empty lists
      self.Logic.pipeline = {}
      self.Logic.completed = {}

      # Change buttons and ProgressBars
      self.ApplyButton.setEnabled(True)
      self.ApplyButton.setText("Run ShapeAnalysisModule")
      self.CLIProgressBars.setParent(None)
      self.Logic.ProgressBar.show()

    # if running, create toolbars
    elif status == 'Running':
      self.CLIProgressBars = qt.QWidget()
      progressbars_layout = qt.QVBoxLayout(self.CLIProgressBars)
      self.progress_layout.addWidget(self.CLIProgressBars)
      self.Logic.ProgressBar.hide()
      for i in range(len(self.Logic.pipeline)):
        progressbars_layout.addWidget(self.Logic.pipeline[i].ProgressBar)

  #
  #   Flip Options
  #
  def onPreviewFlips(self):
    # Creation of a CSV file to load the vtk files in ShapePopulationViewer
    filePathCSV = slicer.app.temporaryPath + '/' + 'PreviewFlips.csv'
    self.Logic.creationCSVFileForSPV(filePathCSV)

    # Launch the CLI ShapePopulationViewer
    parameters = {}
    parameters["CSVFile"] = filePathCSV
    launcherSPV = slicer.modules.launcher
    slicer.cli.run(launcherSPV, None, parameters, wait_for_completion=True)

  def onChangeFlips(self):
    # Run workflow
    if not self.Logic.Node.IsBusy():
      logging.info('Widget: Running ParaToSPHARMMesh to Change Flips')
      self.changeFlips.setText("Cancel")

      self.Logic.addObserver(self.Logic.Node, slicer.vtkMRMLCommandLineModuleNode().StatusModifiedEvent,
                             self.onLogicModified)

      self.Logic.Node.SetStatus(self.Logic.Node.Scheduled)
      # self.Logic.allCaseStartTime = time.time()

      self.Logic.option_pipeline = 'changeFlip'
      self.Logic.ShapeAnalysisCases()

    # Cancel Workflow
    else:
      logging.info("Widget: Cancelling ParaToSPHARMMesh to Change Flips")
      self.changeFlips.setEnabled(False)
      self.Logic.Cancel()

#
# ShapeAnalysisModuleLogic
#

class ShapeAnalysisModuleLogic(ScriptedLoadableModuleLogic, VTKObservationMixin):
  """This class should implement all the actual
  computation done by your module.  The interface
  should be such that other python code can import
  this class and make use of the functionality without
  requiring an instance of the Widget.
  Uses ScriptedLoadableModuleLogic base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """
  def __init__(self, interface):

    VTKObservationMixin.__init__(self)

    self.interface = interface
    self.InputCases = list()
    self.option_pipeline = None
    self.allCaseStartTime = 0

    # Dictionaries
    self.pipeline = {}
    self.completed = {}

    # Status
    self.Node = slicer.vtkMRMLCommandLineModuleNode()
    self.Node.SetStatus(self.Node.Idle)
    self.Node.SetName("ShapeAnalysisModule")
    self.ProgressBar = slicer.qSlicerCLIProgressBar()
    self.ProgressBar.setCommandLineModuleNode(self.Node)
    self.ErrorMessage = 'Unexpected error'

  def ShapeAnalysisCases(self):

    self.InputCases = list()
    self.pipeline = {}

    # Search cases
    inputDirectory = self.interface.GroupProjectInputDirectory.directory.encode('utf-8')
    for file in os.listdir(inputDirectory):
      if file.endswith(".gipl") or file.endswith(".gipl.gz"):
        self.InputCases.append(file)

    # No cases
    if not len(self.InputCases) > 0:
      slicer.util.errorDisplay("No cases found in " + inputDirectory)
      return -1

    # Create pipelines
    else:
      # Init
      for i in range(len(self.InputCases)):
        self.completed[i] = False
        self.pipeline[i] = ShapeAnalysisModulePipeline(i, self.InputCases[i], self.option_pipeline, self.interface)

        self.addObserver(self.pipeline[i].Node, slicer.vtkMRMLCommandLineModuleNode().StatusModifiedEvent,
                       self.onPipelineModified)
      # Logic ready
      self.Node.SetStatus(self.Node.Running)
      # Launch Workflow
      self.startPipeline(0)
      return 0

  def startPipeline(self, id):
    self.pipeline[id].setup()
    self.pipeline[id].Node.SetStatus(self.Node.Scheduled)
    self.pipeline[id].runFirstCLIModule()

  def onPipelineModified(self, pipeline_node, event):
    pipeline_id = None
    current_pipeline = None
    for key, pipeline in self.pipeline.iteritems():
      if pipeline.Node == pipeline_node:
        pipeline_id = key
        current_pipeline = pipeline
    if pipeline_id is None:
      logging.error('Error: Unidentified pipeline modified')
      return -1

    status = pipeline_node.GetStatusString()
    logging.info('-- %s: Case %d: %s', status, pipeline_id, current_pipeline.CaseInput)

    if not pipeline_node.IsBusy():
      self.removeObserver(pipeline_node, slicer.vtkMRMLCommandLineModuleNode().StatusModifiedEvent,
                          self.onPipelineModified)
      # Report time taken to get to the non-busy state for the pipeline
      logging.info('Case %d (%s) took %d sec to run', pipeline_id, current_pipeline.CaseInput,
                    current_pipeline.getPipelineComputationTime())
      statusForNode = None
      # If canceled, stop everything
      if pipeline_node.GetStatusString() == 'Cancelled':
        statusForNode = pipeline_node.GetStatus()
        # If completed, with errors or not
      else:
        # If completed with errors, inform user
        if pipeline_node.GetStatusString() == 'Completed with errors':
          self.ErrorMessage = current_pipeline.ErrorMessage
          logging.error(current_pipeline.ErrorMessage)
        # Then starts next case if it exists
        self.completed[pipeline_id] = True
        # If there is no anymore case
        if self.areAllPipelineCompleted():
          logging.info('All pipelines took: %d sec to run', time.time() - self.allCaseStartTime)
          statusForNode = pipeline_node.GetStatus()
          self.fillTableForFlipOptions()

      # Remove nodes from scene depending if we are
      #  in multicase, and if we want to keep intermediate files
      # self.pipeline[pipeline_id].cleanup()
      # gc.collect()

      if statusForNode is None:
        # Run next pipeline
        self.startPipeline(pipeline_id + 1)
      else:
        self.Node.SetStatus(statusForNode)

  def areAllPipelineCompleted(self):
    for i in range(len(self.completed)):
      if not self.completed[i]:
        return False
    return True

  def Cancel(self):
    self.Node.SetStatus(self.Node.Cancelling)
    for i in range(len(self.pipeline)):
      self.pipeline[i].Cancel()

  # Function to fill the table of the flip options for all the SPHARM mesh output
  #    - Column 0: filename of the SPHARM mesh output vtk file
  #    - Column 1: combobox with the flip corresponding to the output file
  def fillTableForFlipOptions(self):
    table = self.interface.tableWidget_ChoiceOfFlip
    row = 0

    outputDirectory = self.interface.GroupProjectOutputDirectory.directory.encode('utf-8')
    SPHARMMeshOutputDirectory = outputDirectory + "/SPHARMMesh/"
    for basename in self.InputCases:
      filepath = SPHARMMeshOutputDirectory + os.path.splitext(basename)[0] + "SPHARM.vtk"
      if os.path.exists(filepath):
        table.setRowCount(row + 1)
        # Column 0:
        filename = os.path.splitext(os.path.basename(filepath))[0]
        labelVTKFile = qt.QLabel(filename)
        labelVTKFile.setAlignment(0x84)
        table.setCellWidget(row, 0, labelVTKFile)

        # Column 1:
        widget = qt.QWidget()
        layout = qt.QHBoxLayout(widget)
        comboBox = qt.QComboBox()
        comboBox.addItems(['No Flip',
                           'Flip Along Axis of x and y',
                           'Flip Along Axis of y and z',
                           'Flip Along Axis of x and z',
                           'Flip Along Axis of x',
                           'Flip Along Axis of y',
                           'Flip Along Axis of x, y and z',
                           'Flip Along Axis of z'])
        comboBox.setCurrentIndex(self.interface.choiceOfFlip.currentIndex)
        layout.addWidget(comboBox)
        layout.setAlignment(0x84)
        layout.setContentsMargins(0, 0, 0, 0)
        widget.setLayout(layout)
        table.setCellWidget(row, 1, widget)

        row = row + 1

  # Function to create a CSV file containing all the SPHARM mesh output files
  # that the user wants to display in ShapePopultaionViewer in order to check the flip
  def creationCSVFileForSPV(self, filepathCVS):
    # Creation a CSV file with a header 'VTK Files'
    file = open(filepathCVS, 'w')
    cw = csv.writer(file, delimiter=',')
    cw.writerow(['VTK Files'])
    # Add the path of the SPHARM mesh output files
    outputDirectory = self.interface.GroupProjectOutputDirectory.directory.encode('utf-8')
    SPHARMMeshOutputDirectory = outputDirectory + "/SPHARMMesh/"
    for basename in self.InputCases:
      filepath = SPHARMMeshOutputDirectory + os.path.splitext(basename)[0] + "SPHARM.vtk"
      if os.path.exists(filepath):
        cw.writerow([filepath])

  def clearFlipOptionsTable(self):
    table = self.interface.tableWidget_ChoiceOfFlip
    table.clear()
    table.setColumnCount(2)
    table.setHorizontalHeaderLabels([' Files ', ' Choice of Flip '])
    table.setColumnWidth(0, 400)
    horizontalHeader = table.horizontalHeader()
    horizontalHeader.setStretchLastSection(False)
    horizontalHeader.setResizeMode(0, qt.QHeaderView.Stretch)
    horizontalHeader.setResizeMode(1, qt.QHeaderView.ResizeToContents)
    table.verticalHeader().setVisible(False)

#
# ShapeAnalysisModuleMRMLUtility
#
'''
This class harbors all the utility functions to load, add, save and remove mrml nodes
'''

class ShapeAnalysisModuleMRMLUtility(object):

  @staticmethod
  def loadMRMLNode(file_path, file_type ):
    properties = {}
    if file_type == 'LabelMapVolumeFile':
      file_type = 'VolumeFile'
      properties['labelmap'] = True
    node = slicer.util.loadNodeFromFile(file_path, file_type, properties, returnNode=True)
    node = node[1]
    return node

  @staticmethod
  def addnewMRMLNode(node_name, node_type):
    node = slicer.mrmlScene.AddNode(node_type)
    node.SetName(node_name)
    return node

  @staticmethod
  def saveMRMLNode(node, filepath):
    slicer.util.saveNode(node, filepath)

  @staticmethod
  def removeMRMLNode(node):
    slicer.mrmlScene.RemoveNode(node)

#
# ShapeAnalysisModuleNode
#
class ShapeAnalysisModuleNode(object):
  nodes = [None]
  save = [False]
  delete = [True]
  filepaths = [" "]

#
# ShapeAnalysisModulePipeline
#
class ShapeAnalysisModulePipeline(VTKObservationMixin):
  def __init__(self, pipelineID, CaseInput, option_pipeline, interface):

    VTKObservationMixin.__init__(self)

    self.interface = interface

    self.pipelineID = pipelineID
    self.strPipelineID = "_" + str(self.pipelineID)
    self.CaseInput = CaseInput
    self.option_pipeline = option_pipeline
    self.skip_segPostProcess = False
    self.skip_genParaMesh = False
    self.skip_paraToSPHARMMesh = False

    # Pipeline computation time
    self.pipelineStartTime = 0
    self.pipelineEndTime = 0

    # Status
    self.StatusModifiedEvent = slicer.vtkMRMLCommandLineModuleNode().StatusModifiedEvent
    self.Node = slicer.vtkMRMLCommandLineModuleNode()
    self.Node.SetStatus(self.Node.Idle)
    self.Node.SetName('Case ' + str(self.pipelineID))
    self.currentCLINode = None
    self.ProgressBar = slicer.qSlicerCLIProgressBar()
    self.ProgressBar.setCommandLineModuleNode(self.Node)
    self.ProgressBar.setNameVisibility(slicer.qSlicerCLIProgressBar.AlwaysVisible)
    self.ErrorMessage = 'Unexpected error'

  # Return pipeline computation time
  def getPipelineComputationTime(self):
    return (self.pipelineEndTime - self.pipelineStartTime)

  def setupSkipCLIs(self):
    outputDirectory = self.interface.GroupProjectOutputDirectory.directory.encode('utf-8')

    # Skip SegPostProcess ?
    if self.option_pipeline == None:
      if not self.interface.OverwriteGenParaMesh.checkState():
        PostProcessDirectory = outputDirectory + "/PostProcess"
        PostProcessOutputFilepath = PostProcessDirectory + "/" + os.path.splitext(self.CaseInput)[0] + "_pp.gipl"
        if os.path.exists(PostProcessOutputFilepath):
          self.skip_segPostProcess = True

      # Skip GenParaMesh ?
      if not self.interface.OverwriteGenParaMesh.checkState():
        GenParaMeshOutputDirectory = outputDirectory + "/MeshParameters"
        ParaOutputFilepath = GenParaMeshOutputDirectory + "/" + os.path.splitext(self.CaseInput)[0] + "_para.vtk"
        SurfOutputFilepath = GenParaMeshOutputDirectory + "/" + os.path.splitext(self.CaseInput)[0] + "_surf.vtk"
        if os.path.exists(ParaOutputFilepath) and os.path.exists(SurfOutputFilepath):
          self.skip_genParaMesh = True

      # Skip ParaToSPHARMMesh ?
      SPHARMMeshOutputDirectory = outputDirectory + "/SPHARMMesh"
      SPHARMMeshFilepath = SPHARMMeshOutputDirectory + "/" + os.path.splitext(self.CaseInput)[0]
      SPHARMMeshDirectory = os.path.dirname(SPHARMMeshFilepath)
      SPHARMMeshBasename = os.path.basename(SPHARMMeshFilepath)
      if not self.interface.OverwriteParaToSPHARMMesh.checkState():
        if os.path.exists(SPHARMMeshDirectory):
          for file in os.listdir(SPHARMMeshDirectory):
            if not file.find(SPHARMMeshBasename) == -1:
                self.skip_paraToSPHARMMesh = True

    elif self.option_pipeline == 'changeFlip':
      self.skip_segPostProcess = True
      self.skip_genParaMesh = True
      self.skip_paraToSPHARMMesh = False

  def setupGlobalVariables(self):
    # Modules
    self.ID = -1
    self.slicerModule = {}
    self.moduleParameters = {}

    # Nodes
    self.nodeDictionary = {}

  def setupModule(self, module, cli_parameters):
    self.slicerModule[self.ID] = module
    self.moduleParameters[self.ID] = cli_parameters

  def setupNode(self, id, cli_nodes, cli_filepaths, cli_saveOutput, cli_deleteOutput):
    self.nodeDictionary[id] = ShapeAnalysisModuleNode()
    self.nodeDictionary[id].nodes = cli_nodes
    self.nodeDictionary[id].filepaths = cli_filepaths
    self.nodeDictionary[id].save = cli_saveOutput
    self.nodeDictionary[id].delete = cli_deleteOutput

  def setup(self):
    # Initialization of global variables
    self.setupGlobalVariables()
    self.setupSkipCLIs()

    inputDirectory = self.interface.GroupProjectInputDirectory.directory.encode('utf-8')
    outputDirectory = self.interface.GroupProjectOutputDirectory.directory.encode('utf-8')

    ## Post Processed Segmentation
    cli_nodes = list() # list of the nodes used in the Post Processed Segmentation step

    cli_filepaths = list() # list of the node filepaths used in the Post Processed Segmentation step
    PostProcessDirectory = outputDirectory + "/PostProcess"
    PostProcessOutputFilepath = PostProcessDirectory + "/" + os.path.splitext(self.CaseInput)[0] + "_pp.gipl"

    if not self.skip_segPostProcess:
      # Setup of the parameters od the CLI
      self.ID += 1

      cli_parameters = {}
      inputFilepath = inputDirectory + '/' + self.CaseInput

      input_node = ShapeAnalysisModuleMRMLUtility.loadMRMLNode(inputFilepath, 'LabelMapVolumeFile')
      cli_parameters["fileName"] = input_node

      pp_output_node = ShapeAnalysisModuleMRMLUtility.addnewMRMLNode("output_PostProcess", slicer.vtkMRMLLabelMapVolumeNode())
      cli_parameters["outfileName"] = pp_output_node.GetID()

      if self.interface.RescaleSegPostProcess.checkState():
        cli_parameters["scaleOn"] = True
        cli_parameters["spacing_vect"] = str(self.interface.sx.value) + "," + str(self.interface.sy.value) + "," + str(self.interface.sz.value)
      cli_parameters["label"] = self.interface.ValueLabelNumber.value
      if self.interface.Debug.checkState():
        cli_parameters["debug"] = True

      #    Advanced parameters
      if self.interface.GaussianFiltering.checkState():
        cli_parameters["gaussianOn"] = True
        cli_parameters["variance_vect"] = str(self.interface.VarianceX.value) + "," + str(self.interface.VarianceY.value) + "," + str(self.interface.VarianceZ.value)

      self.setupModule(slicer.modules.segpostprocessclp, cli_parameters)

      # Setup of the nodes created by the CLI
      #    Creation of a folder in the output folder : PostProcess
      if not os.path.exists(PostProcessDirectory):
        os.makedirs(PostProcessDirectory)

      cli_nodes.append(input_node)
      cli_nodes.append(pp_output_node)
      cli_filepaths.append(inputFilepath)
      cli_filepaths.append(PostProcessOutputFilepath)

      self.setupNode(1, cli_nodes, cli_filepaths, [False,True], [True,True])

    else:
      # Setup of the nodes which will be use by the next CLI
      pp_output_node = ShapeAnalysisModuleMRMLUtility.loadMRMLNode(PostProcessOutputFilepath, 'LabelMapVolumeFile')

      cli_filepaths.append(PostProcessOutputFilepath)
      cli_nodes.append(pp_output_node)

      self.setupNode(1, cli_nodes, cli_filepaths, [False], [True])


    ## Generate Mesh Parameters
    cli_nodes = list() # list of the nodes used in the Generate Mesh Parameters step

    cli_filepaths = list() # list of the node filepaths used in the Generate Mesh Parameters step
    GenParaMeshOutputDirectory = outputDirectory + "/MeshParameters"
    ParaOutputFilepath = GenParaMeshOutputDirectory + "/" + os.path.splitext(self.CaseInput)[0] + "_para.vtk"
    SurfOutputFilepath = GenParaMeshOutputDirectory + "/" + os.path.splitext(self.CaseInput)[0] + "_surf.vtk"

    if not self.skip_genParaMesh:
      # Setup of the parameters od the CLI

      self.ID += 1

      cli_parameters = {}
      cli_parameters["infile"] = pp_output_node

      para_output_model = ShapeAnalysisModuleMRMLUtility.addnewMRMLNode("output_para", slicer.vtkMRMLModelNode())
      cli_parameters["outParaName"] = para_output_model

      surfmesh_output_model = ShapeAnalysisModuleMRMLUtility.addnewMRMLNode("output_surfmesh", slicer.vtkMRMLModelNode())
      cli_parameters["outSurfName"] = surfmesh_output_model

      cli_parameters["numIterations"] = 5 #self.interface.NumberofIterations.value
      if self.interface.Debug.checkState():
        cli_parameters["debug"] = True

      self.setupModule(slicer.modules.genparameshclp, cli_parameters)


      # Setup of the nodes created by the CLI
      #    Creation of a folder in the output folder : GenerateMeshParameters
      if not os.path.exists(GenParaMeshOutputDirectory):
        os.makedirs(GenParaMeshOutputDirectory)

      cli_nodes.append(para_output_model)
      cli_nodes.append(surfmesh_output_model)
      cli_filepaths.append(ParaOutputFilepath)
      cli_filepaths.append(SurfOutputFilepath)

      self.setupNode(2, cli_nodes, cli_filepaths, [True,True], [True,True])


    else:
      # Setup of the nodes which will be use by the next CLI
      para_output_model = ShapeAnalysisModuleMRMLUtility.loadMRMLNode(ParaOutputFilepath, 'ModelFile')
      surfmesh_output_model = ShapeAnalysisModuleMRMLUtility.loadMRMLNode(SurfOutputFilepath, 'ModelFile')

      cli_nodes.append(para_output_model)
      cli_nodes.append(surfmesh_output_model)
      cli_filepaths.append(ParaOutputFilepath)
      cli_filepaths.append(SurfOutputFilepath)

      self.setupNode(2, cli_nodes, cli_filepaths, [False, False], [True, True])


    ##  Parameters to SPHARM Mesh
    SPHARMMeshOutputDirectory = outputDirectory + "/SPHARMMesh"
    SPHARMMeshFilepath = SPHARMMeshOutputDirectory + "/" + os.path.splitext(self.CaseInput)[0]
    if not self.skip_paraToSPHARMMesh:

      # Setup of the parameters od the CLI
      self.ID += 1

      cli_parameters = {}

      cli_parameters["inParaFile"] = para_output_model

      cli_parameters["inSurfFile"] = surfmesh_output_model

      #    Creation of a folder in the output folder : SPHARMMesh
      if not os.path.exists(SPHARMMeshOutputDirectory):
        os.makedirs(SPHARMMeshOutputDirectory)
      cli_parameters["outbase"] = SPHARMMeshFilepath

      cli_parameters["subdivLevel"] = self.interface.SubdivLevelValue.value
      cli_parameters["spharmDegree"] = self.interface.SPHARMDegreeValue.value
      cli_parameters["thetaIteration"] = self.interface.thetaIterationValue.value
      cli_parameters["phiIteration"] = self.interface.phiIterationValue.value
      if self.interface.medialMesh.checkState():
        cli_parameters["medialMesh"] = True
      if self.interface.Debug.checkState():
        cli_parameters["debug"] = True

      #   Advanced parameters
      if self.option_pipeline == 'changeFlip':
        # Recovery of the flip choisen by the user
        row = self.pipelineID
        widget = self.interface.tableWidget_ChoiceOfFlip.cellWidget(row, 1)
        tuple = widget.children()
        comboBox = qt.QComboBox()
        comboBox = tuple[1]
        cli_parameters["finalFlipIndex"] = comboBox.currentIndex
      else:
        cli_parameters["finalFlipIndex"] = self.interface.choiceOfFlip.currentIndex  # 1 = flip along axes of x &amp; y,
        # 2 = flip along y &amp; z,
        # 3 = flip along x &amp; z
        # 4 = flip along x,
        # 5 = flip along y,
        # 6 = flip along x &amp; y &amp; z,
        # 7 = flip along z  where y is the smallest, x is the second smallest and z is the long axis of the ellipsoid

      self.setupModule(slicer.modules.paratospharmmeshclp, cli_parameters)

  # Run the CLI for the current module
  def runCLIModule(self):
    # Create and run CLI
    module = self.slicerModule[self.ID]
    cli_node = self.createCLINode(module)
    self.setCurrentCLINode(cli_node) # ProgressBar
    self.addObserver(cli_node, self.StatusModifiedEvent, self.onCLIModuleModified)
    slicer.cli.run(module, cli_node, self.moduleParameters[self.ID], wait_for_completion = False)

  # Call the next module
  def runNextCLIModule(self):
    self.ID += 1
    self.runCLIModule()

  # Start the pipeline
  def runFirstCLIModule(self):
    if len(self.slicerModule) > 0:
      self.ID = 0
      self.pipelineStartTime = time.time()
      self.runCLIModule()
    else:
      logging.info('Slicer Module queue is empty, nothing to do')
      self.Node.SetStatus(self.Node.Completed)

  def createCLINode(self, module):
    cli_node = slicer.cli.createNode(module)
    cli_node_name = "Case " + str(self.pipelineID) + ": " + os.path.splitext(self.CaseInput)[0] + " - step " + str(self.ID) + ": " + module.title
    cli_node.SetName(cli_node_name)
    return cli_node

  def onCLIModuleModified(self, cli_node, event):
    logging.info('-- %s : %s', cli_node.GetStatusString(), cli_node.GetName())
    if not cli_node.IsBusy():
      if platform.system() != 'Windows':
        self.removeObserver(cli_node, self.StatusModifiedEvent,
                            self.onCLIModuleModified)
        statusForNode =  None
      if cli_node.GetStatusString() == 'Completed':
        if self.ID == len(self.slicerModule) - 1:
          cli_node_name = "Case " + str(self.pipelineID) + ": " + os.path.splitext(self.CaseInput)[0]
          cli_node.SetName(cli_node_name)
          self.setCurrentCLINode(cli_node)
          statusForNode = cli_node.GetStatus()

      elif cli_node.GetStatusString() == 'Cancelled':
        self.ErrorMessage = cli_node.GetName() + " Cancelled (" + self.CaseInput + "), cancelling the pipeline"
        #Interrupt pipeline
        statusForNode = cli_node.GetStatus()

      else:
        # Create Error Message
        if cli_node.GetStatusString() == 'Completed with errors':
          self.ErrorMessage = cli_node.GetName() + " completed with errors"
          # Interrompt Pipeline
          statusForNode = cli_node.GetStatus()

      if statusForNode is None:
        self.runNextCLIModule()
      else:
        self.saveNodes()
        self.deleteNodes()
        self.pipelineEndTime = time.time()
        self.Node.SetStatus(statusForNode)

  def saveNodes(self):
    for id in self.nodeDictionary.keys():
      for i in range(len(self.nodeDictionary[id].save)):
        save = self.nodeDictionary[id].save
        node = self.nodeDictionary[id].nodes
        filepaths = self.nodeDictionary[id].filepaths
        if save[i] == True:
          ShapeAnalysisModuleMRMLUtility.saveMRMLNode( node[i], filepaths[i] )

  def deleteNodes(self):
    for id in self.nodeDictionary.keys():
      for i in range(len(self.nodeDictionary[id].save)):
        delete = self.nodeDictionary[id].delete
        node = self.nodeDictionary[id].nodes
        if delete[i] == True:
          ShapeAnalysisModuleMRMLUtility.removeMRMLNode( node[i] )

  def setCurrentCLINode(self, cli_node):
    self.ProgressBar.setCommandLineModuleNode(cli_node)
    self.currentCLINode = cli_node

  def Cancel(self):
    if self.currentCLINode:
      self.currentCLINode.Cancel()

class ShapeAnalysisModuleTest(ScriptedLoadableModuleTest):
  """
  This is the test case for your scripted module.
  Uses ScriptedLoadableModuleTest base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setUp(self):
    """ Do whatever is needed to reset the state - typically a scene clear will be enough.
    """
    slicer.mrmlScene.Clear(0)

  def runTest(self):
    """Run as few or as many tests as needed here.
    """
    self.setUp()
    self.delayDisplay(' Tests Passed! ')