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
    self.progressbars_layout = None

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
    self.CollapsibleButton_GroupProjectIO = self.getWidget('CollapsibleButton_GroupProjectIO')
    self.GroupProjectInputDirectory = self.getWidget('DirectoryButton_GroupProjectInputDirectory')
    self.GroupProjectOutputDirectory = self.getWidget('DirectoryButton_GroupProjectOutputDirectory')
    self.Debug = self.getWidget('checkBox_Debug')
    #   Post Processed Segmentation
    self.CollapsibleButton_SegPostProcess = self.getWidget('CollapsibleButton_SegPostProcess')
    self.OverwriteSegPostProcess = self.getWidget('checkBox_OverwriteSegPostProcess')
    self.label_RescaleSegPostProcess = self.getWidget('label_RescaleSegPostProcess')
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
    self.CollapsibleButton_GenParaMesh = self.getWidget('CollapsibleButton_GenParaMesh')
    self.OverwriteGenParaMesh = self.getWidget('checkBox_OverwriteGenParaMesh')
    self.NumberofIterations = self.getWidget('SliderWidget_NumberofIterations')
    #   Parameters to SPHARM Mesh
    self.CollapsibleButton_ParaToSPHARMMesh = self.getWidget('CollapsibleButton_ParaToSPHARMMesh')
    self.OverwriteParaToSPHARMMesh = self.getWidget('checkBox_OverwriteParaToSPHARMMesh')
    self.SubdivLevelValue = self.getWidget('SliderWidget_SubdivLevelValue')
    self.SPHARMDegreeValue = self.getWidget('SliderWidget_SPHARMDegreeValue')
    self.thetaIterationValue = self.getWidget('spinBox_thetaIterationValue')
    self.phiIterationValue = self.getWidget('spinBox_phiIterationValue')
    self.medialMesh = self.getWidget('checkBox_medialMesh')
    #   Advanced Post Processed Segmentation
    self.CollapsibleButton_AdvancedPostProcessedSegmentation = self.getWidget('CollapsibleButton_AdvancedPostProcessedSegmentation')
    self.GaussianFiltering = self.getWidget('checkBox_GaussianFiltering')
    self.label_VarianceX = self.getWidget('label_VarianceX')
    self.VarianceX = self.getWidget('SliderWidget_VarianceX')
    self.label_VarianceY = self.getWidget('label_VarianceY')
    self.VarianceY = self.getWidget('SliderWidget_VarianceY')
    self.label_VarianceZ = self.getWidget('label_VarianceZ')
    self.VarianceZ = self.getWidget('SliderWidget_VarianceZ')
    #   Advanced Parameters to SPHARM Mesh
    self.CollapsibleButton_AdvancedParametersToSPHARMMesh = self.getWidget('CollapsibleButton_AdvancedParametersToSPHARMMesh')
    self.useRegTemplate = self.getWidget('checkBox_useRegTemplate')
    self.label_regTemplate = self.getWidget('label_regTemplate')
    self.regTemplate = self.getWidget('PathLineEdit_regTemplate')
    self.useFlipTemplate = self.getWidget('checkBox_useFlipTemplate')
    self.label_flipTemplate = self.getWidget('label_flipTemplate')
    self.flipTemplate = self.getWidget('PathLineEdit_flipTemplate')
    self.choiceOfFlip = self.getWidget('comboBox_choiceOfFlip')
    self.sameFlipForAll = self.getWidget('checkBox_sameFlipForAll')
    self.tableWidget_ChoiceOfFlip = self.getWidget('tableWidget_ChoiceOfFlip')
    #   Visualization
    self.CollapsibleButton_Visualization = self.getWidget('CollapsibleButton_Visualization')
    self.visualizationInSPV = self.getWidget('pushButton_visualizationInSPV')
    self.CheckableComboBox_visualization = self.getWidget('CheckableComboBox_visualization')
    self.tableWidget_visualization = self.getWidget('tableWidget_visualization')
    #   Apply CLIs
    self.ApplyButton = self.getWidget('applyButton')
    self.progress_layout = self.getWidget('progress_layout')

    # Connections
    #   Group Project IO
    self.CollapsibleButton_GroupProjectIO.connect('clicked()',
                                                   lambda: self.onSelectedCollapsibleButtonOpen(
                                                     self.CollapsibleButton_GroupProjectIO))
    self.GroupProjectInputDirectory.connect('directoryChanged(const QString &)', self.onInputDirectoryChanged)


    #   Post Processed Segmentation
    self.CollapsibleButton_SegPostProcess.connect('clicked()',
                                                  lambda: self.onSelectedCollapsibleButtonOpen(
                                                    self.CollapsibleButton_SegPostProcess))
    self.OverwriteSegPostProcess.connect('clicked(bool)', self.onOverwriteFilesSegPostProcess)
    self.RescaleSegPostProcess.connect('stateChanged(int)', self.onSelectSpacing)
    self.LabelState.connect('clicked(bool)', self.onSelectValueLabelNumber)
    #   Generate Mesh Parameters
    self.CollapsibleButton_GenParaMesh.connect('clicked()',
                                                  lambda: self.onSelectedCollapsibleButtonOpen(
                                                    self.CollapsibleButton_GenParaMesh))

    self.OverwriteGenParaMesh.connect('clicked(bool)', self.onOverwriteFilesGenParaMesh)
    #   Parameters to SPHARM Mesh
    self.CollapsibleButton_ParaToSPHARMMesh.connect('clicked()',
                                                  lambda: self.onSelectedCollapsibleButtonOpen(
                                                    self.CollapsibleButton_ParaToSPHARMMesh))
    self.OverwriteParaToSPHARMMesh.connect('clicked(bool)', self.onOverwriteFilesParaToSPHARMMesh)
    #   Advanced Post Processed Segmentation
    self.CollapsibleButton_AdvancedPostProcessedSegmentation.connect('clicked()',
                                                  lambda: self.onSelectedCollapsibleButtonOpen(
                                                    self.CollapsibleButton_AdvancedPostProcessedSegmentation))
    self.GaussianFiltering.connect('clicked(bool)', self.onSelectGaussianVariance)
    #   Advanced Parameters to SPHARM Mesh
    self.CollapsibleButton_AdvancedParametersToSPHARMMesh.connect('clicked()',
                                                  lambda: self.onSelectedCollapsibleButtonOpen(
                                                    self.CollapsibleButton_AdvancedParametersToSPHARMMesh))
    self.useRegTemplate.connect('clicked(bool)', self.onEnableRegTemplate)
    self.useFlipTemplate.connect('clicked(bool)', self.onEnableFlipTemplate)
    self.sameFlipForAll.connect('clicked(bool)', self.onEnableFlipChoices)

    #   Visualization
    self.CollapsibleButton_Visualization.connect('clicked()',
                                                  lambda: self.onSelectedCollapsibleButtonOpen(
                                                    self.CollapsibleButton_Visualization))
    self.CheckableComboBox_visualization.connect('checkedIndexesChanged()', self.onCheckableComboBoxValueChanged)
    self.visualizationInSPV.connect('clicked(bool)', self.onPreviewFlips)

    #   Apply CLIs
    self.ApplyButton.connect('clicked(bool)', self.onApplyButton)

    # Widget Configuration
    #     Table for the Flip Options
    self.tableWidget_ChoiceOfFlip.setColumnCount(2)
    self.tableWidget_ChoiceOfFlip.setHorizontalHeaderLabels([' Input Files ', ' Choice of Flip '])
    self.tableWidget_ChoiceOfFlip.setColumnWidth(0, 400)
    horizontalHeader = self.tableWidget_ChoiceOfFlip.horizontalHeader()
    horizontalHeader.setStretchLastSection(False)
    horizontalHeader.setResizeMode(0, qt.QHeaderView.Stretch)
    horizontalHeader.setResizeMode(1, qt.QHeaderView.ResizeToContents)
    self.tableWidget_ChoiceOfFlip.verticalHeader().setVisible(False)

    #     Progress Bar
    self.progress_layout.addWidget(self.Logic.ProgressBar)

    #     Table for the visualization in SPV
    self.tableWidget_visualization.setColumnCount(2)
    self.tableWidget_visualization.setHorizontalHeaderLabels([' VTK Files ', ' Visualization '])
    self.tableWidget_visualization.setColumnWidth(0, 400)
    horizontalHeader = self.tableWidget_visualization.horizontalHeader()
    horizontalHeader.setStretchLastSection(False)
    horizontalHeader.setResizeMode(0, qt.QHeaderView.Stretch)
    horizontalHeader.setResizeMode(1, qt.QHeaderView.ResizeToContents)
    self.tableWidget_visualization.verticalHeader().setVisible(False)


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

  # Only one tab can be display at the same time:
  #   When one tab is opened all the other tabs are closed
  def onSelectedCollapsibleButtonOpen(self, selectedCollapsibleButton):
    if selectedCollapsibleButton.isChecked():
      collapsibleButtonList = [self.CollapsibleButton_GroupProjectIO,
                               self.CollapsibleButton_SegPostProcess,
                               self.CollapsibleButton_GenParaMesh,
                               self.CollapsibleButton_ParaToSPHARMMesh,
                               self.CollapsibleButton_AdvancedPostProcessedSegmentation,
                               self.CollapsibleButton_AdvancedParametersToSPHARMMesh,
                               self.CollapsibleButton_Visualization]
      for collapsibleButton in collapsibleButtonList:
        collapsibleButton.setChecked(False)
      selectedCollapsibleButton.setChecked(True)

  #
  #   Group Project IO
  #
  def onInputDirectoryChanged(self):
    #  Possible extension
    exts = [".gipl", ".gipl.gz", ".mgh", ".mgh,gz", ".nii", ".nii.gz",".nrrd", "vtk", "vtp" "hdr", "mhd"]

    # Search cases
    self.Logic.InputCases = []
    inputDirectory = self.GroupProjectInputDirectory.directory.encode('utf-8')
    for file in os.listdir(inputDirectory):
      for ext in exts:
        if file.endswith(ext):
          self.Logic.InputCases.append(file)
          if file.endswith(".nii") or file.endswith(".nii.gz"):
            self.RescaleSegPostProcess.setCheckState(qt.Qt.Unchecked)
            self.label_RescaleSegPostProcess.enabled = False
            self.RescaleSegPostProcess.enabled = False

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

  def onOverwriteFilesSegPostProcess(self):
    if self.OverwriteSegPostProcess.checkState():
      slicer.util.errorDisplay("Applying the overwrite option to Post Processed Segmentation step will also apply to the next steps")
      self.OverwriteGenParaMesh.blockSignals(True)
      self.OverwriteGenParaMesh.setCheckState(qt.Qt.Checked)
      self.OverwriteParaToSPHARMMesh.setCheckState(qt.Qt.Checked)
      self.OverwriteGenParaMesh.blockSignals(False)

  #
  #   Generate Mesh Parameters
  #
  def onOverwriteFilesGenParaMesh(self):
    if not self.OverwriteGenParaMesh.checkState():
      if self.OverwriteSegPostProcess.checkState():
        self.OverwriteGenParaMesh.blockSignals(True)
        self.OverwriteGenParaMesh.setCheckState(qt.Qt.Checked)
        slicer.util.errorDisplay("The overwrite option need to be applied to this step as it is set for the previous step")
        self.OverwriteGenParaMesh.blockSignals(False)
    else:
      slicer.util.errorDisplay("Applying the overwrite option to Generate Mesh Parameters step will also apply to the next steps")
      self.OverwriteParaToSPHARMMesh.setCheckState(qt.Qt.Checked)

  #
  #   Parameters to SPHARM Mesh
  #
  def onOverwriteFilesParaToSPHARMMesh(self):
    if not self.OverwriteParaToSPHARMMesh.checkState():
      if self.OverwriteSegPostProcess.checkState() or self.OverwriteGenParaMesh.checkState():
        self.OverwriteSegPostProcess.blockSignals(True)
        self.OverwriteGenParaMesh.blockSignals(True)
        slicer.util.errorDisplay("The overwrite option need to be applied to this step as it is set for the previous step")
        self.OverwriteParaToSPHARMMesh.setCheckState(qt.Qt.Checked)
        self.OverwriteSegPostProcess.blockSignals(False)
        self.OverwriteGenParaMesh.blockSignals(False)

  #
  #   Advanced Post Processed Segmentation
  #
  def onSelectGaussianVariance(self):
    self.label_VarianceX.enabled = self.GaussianFiltering.checkState()
    self.VarianceX.enabled = self.GaussianFiltering.checkState()
    self.label_VarianceY.enabled = self.GaussianFiltering.checkState()
    self.VarianceY.enabled = self.GaussianFiltering.checkState()
    self.label_VarianceZ.enabled = self.GaussianFiltering.checkState()
    self.VarianceZ.enabled = self.GaussianFiltering.checkState()

  #
  #   Advanced Parameters to SPHARM Mesh
  #
  def onEnableFlipChoices(self):
    self.choiceOfFlip.enabled = self.sameFlipForAll.checkState()
    self.tableWidget_ChoiceOfFlip.enabled = not self.sameFlipForAll.checkState()
    if not self.sameFlipForAll.checkState():
      self.Logic.fillTableForFlipOptions()

  def onEnableRegTemplate(self):
    self.label_regTemplate.enabled = self.useRegTemplate.checkState()
    self.regTemplate.enabled = self.useRegTemplate.checkState()

  def onEnableFlipTemplate(self):
    self.label_flipTemplate.enabled = self.useFlipTemplate.checkState()
    self.flipTemplate.enabled = self.useFlipTemplate.checkState()

  #
  #   Apply CLIs
  #
  def onApplyButton(self):
    # Run workflow
    if not self.Logic.Node.IsBusy():

      # Check the registration template file
      if self.useRegTemplate.checkState():
        if not os.path.exists(self.regTemplate.currentPath) or not self.regTemplate.currentPath.endswith(".vtk"):
          slicer.util.errorDisplay("Invalid registration template file in Advanced Parameters to SPHARM Mesh Tab")
          return

      # Check the flip template file
      if self.useFlipTemplate.checkState():
        if not os.path.exists(self.flipTemplate.currentPath) or not self.flipTemplate.currentPath.endswith(".coef"):
          slicer.util.errorDisplay("Invalid flip template file in Advanced Parameters to SPHARM Mesh Tab")
          return

      # Empty the output folders if the options overwrite is checked
      self.Logic.cleanOutputFolders()

      logging.info('Widget: Running ShapeAnalysisModule')
      self.ApplyButton.setText("Cancel")

      self.Logic.addObserver(self.Logic.Node, slicer.vtkMRMLCommandLineModuleNode().StatusModifiedEvent,
                             self.onLogicModified)

      self.Logic.Node.SetStatus(self.Logic.Node.Scheduled)
      self.Logic.allCaseStartTime = time.time()

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
      if status == 'Completed with errors' or status == 'Cancelled':
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

    # if running, create toolbars
    elif status == 'Running':
      self.Logic.ProgressBar.show()
      if self.progressbars_layout:
        self.CLIProgressBars.hide()
      self.CLIProgressBars = ctk.ctkCollapsibleGroupBox()
      self.CLIProgressBars.setTitle('Detail')
      self.progress_layout.addWidget(self.CLIProgressBars)
      self.progressbars_layout = qt.QVBoxLayout(self.CLIProgressBars)
      for i in range(len(self.Logic.pipeline)):
        self.progressbars_layout.addWidget(self.Logic.pipeline[i].ProgressBar)

  def onCheckableComboBoxValueChanged(self):
    currentText = self.CheckableComboBox_visualization.currentText
    currentIndex = self.CheckableComboBox_visualization.currentIndex
    currentItem = self.CheckableComboBox_visualization.model().item(currentIndex, 0)

    # ******* Update the CheckableComboBox ******* #
    self.CheckableComboBox_visualization.blockSignals(True)
    list = self.CheckableComboBox_visualization.model()
    allIndex = self.CheckableComboBox_visualization.findText("All Models")
    allItem = list.item(allIndex, 0)
    allSPHARMIndex = self.CheckableComboBox_visualization.findText("All SPHARM Models")
    allSPHARMItem = list.item(allSPHARMIndex, 0)
    allSPHARMEllalignIndex = self.CheckableComboBox_visualization.findText("All SPHARM Ellipse Aligned Models")
    allSPHARMEllalignItem = list.item(allSPHARMEllalignIndex, 0)
    allSPHARMMesdialMeshesIndex = self.CheckableComboBox_visualization.findText("All SPHARM Medial Meshes")
    if not allSPHARMMesdialMeshesIndex == -1:
      allSPHARMMesdialMeshesItem = list.item(allSPHARMMesdialMeshesIndex, 0)
    allSPHARMProcrustesAffineModelsIndex = self.CheckableComboBox_visualization.findText("All SPHARM Procrustes Affine Models")
    if not allSPHARMProcrustesAffineModelsIndex == -1:
      allSPHARMProcrustesAffineModelsItem = list.item(allSPHARMProcrustesAffineModelsIndex, 0)

    c = True

    if currentText == "All Models":
      for i in range(list.rowCount()):
        item = list.item(i, 0)
        item.setCheckState(currentItem.checkState())

    elif currentText == "All SPHARM Models":
      for i in range(3, list.rowCount()):
        item = list.item(i, 0)
        if not item.text().find("SPHARM Models") == -1:
          item.setCheckState(currentItem.checkState())
      if allSPHARMEllalignItem.checkState():
        if allSPHARMMesdialMeshesIndex == -1 and allSPHARMProcrustesAffineModelsIndex == -1:
            allItem.setCheckState(currentItem.checkState())
        elif not allSPHARMMesdialMeshesIndex == -1 and not allSPHARMProcrustesAffineModelsIndex == -1:
          if allSPHARMMesdialMeshesItem.checkState() and allSPHARMProcrustesAffineModelsItem.checkState():
            allItem.setCheckState(currentItem.checkState())
        elif not allSPHARMMesdialMeshesIndex == -1 and allSPHARMProcrustesAffineModelsIndex == -1:
          if allSPHARMMesdialMeshesItem.checkState():
            allItem.setCheckState(currentItem.checkState())
        elif allSPHARMMesdialMeshesIndex == -1 and not allSPHARMProcrustesAffineModelsIndex == -1:
          if allSPHARMProcrustesAffineModelsItem.checkState():
            allItem.setCheckState(currentItem.checkState())

    elif currentText == "All SPHARM Ellipse Aligned Models":
      for i in range(3,list.rowCount()):
        item = list.item(i, 0)
        if not item.text().find("SPHARM Ellipse Aligned Models") == -1:
          item.setCheckState(currentItem.checkState())
      if allSPHARMItem.checkState():
        if allSPHARMMesdialMeshesIndex == -1 and allSPHARMProcrustesAffineModelsIndex == -1:
            allItem.setCheckState(currentItem.checkState())
        elif not allSPHARMMesdialMeshesIndex == -1 and not allSPHARMProcrustesAffineModelsIndex == -1:
          if allSPHARMMesdialMeshesItem.checkState() and allSPHARMProcrustesAffineModelsItem.checkState():
            allItem.setCheckState(currentItem.checkState())
        elif not allSPHARMMesdialMeshesIndex == -1 and allSPHARMProcrustesAffineModelsIndex == -1:
          if allSPHARMMesdialMeshesItem.checkState():
            allItem.setCheckState(currentItem.checkState())
        elif allSPHARMMesdialMeshesIndex == -1 and not allSPHARMProcrustesAffineModelsIndex == -1:
          if allSPHARMProcrustesAffineModelsItem.checkState():
            allItem.setCheckState(currentItem.checkState())

    elif currentText == "All SPHARM Medial Meshes":
      for i in range(3,list.rowCount()):
        item = list.item(i, 0)
        if not item.text().find("SPHARM Medial Meshes") == -1:
          item.setCheckState(currentItem.checkState())
      if allSPHARMItem.checkState() and allSPHARMEllalignItem.checkState():
        if not allSPHARMProcrustesAffineModelsIndex == -1:
          if allSPHARMProcrustesAffineModelsItem.checkState():
            allItem.setCheckState(currentItem.checkState())

    elif currentText == "All SPHARM Procrustes Affine Models":
      for i in range(3, list.rowCount()):
        item = list.item(i, 0)
        if not item.text().find("SPHARM Procrustes Affine Models") == -1:
          item.setCheckState(currentItem.checkState())
      if allSPHARMItem.checkState() and allSPHARMEllalignItem.checkState():
        if not allSPHARMMesdialMeshesIndex == -1:
          if allSPHARMMesdialMeshesItem.checkState():
            allItem.setCheckState(currentItem.checkState())

    elif not currentText.find("SPHARM Ellipse Aligned Models") == -1:
      if not currentItem.checkState():
        allSPHARMEllalignItem.setCheckState(qt.Qt.Unchecked)
        allItem.setCheckState(qt.Qt.Unchecked)
      else:
        for i in range(3,list.rowCount()):
          item = list.item(i, 0)
          if not item.text().find(" - SPHARM Ellipse Aligned Models") == -1:
            if not item.checkState():
              c = False
        if c:
          allSPHARMEllalignItem.setCheckState(qt.Qt.Checked)
          if allSPHARMItem.checkState():
            if allSPHARMMesdialMeshesIndex == -1 and allSPHARMProcrustesAffineModelsIndex == -1:
              allItem.setCheckState(currentItem.checkState())
            elif not allSPHARMMesdialMeshesIndex == -1 and not allSPHARMProcrustesAffineModelsIndex == -1:
              if allSPHARMMesdialMeshesItem.checkState() and allSPHARMProcrustesAffineModelsItem.checkState():
                allItem.setCheckState(currentItem.checkState())
            elif not allSPHARMMesdialMeshesIndex == -1 and allSPHARMProcrustesAffineModelsIndex == -1:
              if allSPHARMMesdialMeshesItem.checkState():
                allItem.setCheckState(currentItem.checkState())
            elif allSPHARMMesdialMeshesIndex == -1 and not allSPHARMProcrustesAffineModelsIndex == -1:
              if allSPHARMProcrustesAffineModelsItem.checkState():
                allItem.setCheckState(currentItem.checkState())

    elif not currentText.find("SPHARM Models") == -1:
      if not currentItem.checkState():
        allSPHARMItem.setCheckState(qt.Qt.Unchecked)
        allItem.setCheckState(qt.Qt.Unchecked)
      else:
        for i in range(3,list.rowCount()):
          item = list.item(i, 0)
          if not item.text().find(" - SPHARM Models") == -1:
            if not item.checkState():
              c = False
        if c:
          allSPHARMItem.setCheckState(qt.Qt.Checked)
          if allSPHARMEllalignItem.checkState():
            if allSPHARMMesdialMeshesIndex == -1 and allSPHARMProcrustesAffineModelsIndex == -1:
              allItem.setCheckState(currentItem.checkState())
            elif not allSPHARMMesdialMeshesIndex == -1 and not allSPHARMProcrustesAffineModelsIndex == -1:
              if allSPHARMMesdialMeshesItem.checkState() and allSPHARMProcrustesAffineModelsItem.checkState():
                allItem.setCheckState(currentItem.checkState())
            elif not allSPHARMMesdialMeshesIndex == -1 and allSPHARMProcrustesAffineModelsIndex == -1:
              if allSPHARMMesdialMeshesItem.checkState():
                allItem.setCheckState(currentItem.checkState())
            elif allSPHARMMesdialMeshesIndex == -1 and not allSPHARMProcrustesAffineModelsIndex == -1:
              if allSPHARMProcrustesAffineModelsItem.checkState():
                allItem.setCheckState(currentItem.checkState())

    elif not currentText.find("SPHARM Medial Meshes") == -1:
      if not currentItem.checkState():
        allSPHARMMesdialMeshesItem.setCheckState(qt.Qt.Unchecked)
        allItem.setCheckState(qt.Qt.Unchecked)
      else:
        for i in range(4,list.rowCount()):
          item = list.item(i, 0)
          if not item.text().find(" - SPHARM Medial Meshes") == -1:
            if not item.checkState():
              c = False
        if c:
          allSPHARMMesdialMeshesItem.setCheckState(qt.Qt.Checked)
          if allSPHARMEllalignItem.checkState() and allSPHARMItem.checkState():
            if not allSPHARMProcrustesAffineModelsIndex == -1:
              if allSPHARMProcrustesAffineModelsItem.checkState():
                allItem.setCheckState(currentItem.checkState())

    elif not currentText.find("SPHARM Procrustes Affine Models") == -1:
      if not currentItem.checkState():
        allSPHARMProcrustesAffineModelsItem.setCheckState(qt.Qt.Unchecked)
        allItem.setCheckState(qt.Qt.Unchecked)
      else:
        for i in range(4,list.rowCount()):
          item = list.item(i, 0)
          if not item.text().find(" - SPHARM Procrustes Affine Models") == -1:
            if not item.checkState():
              c = False
        if c:
          allSPHARMProcrustesAffineModelsItem.setCheckState(qt.Qt.Checked)
          if allSPHARMEllalignItem.checkState() and allSPHARMItem.checkState():
            if not allSPHARMMesdialMeshesIndex == -1:
              if allSPHARMMesdialMeshesItem.checkState():
                allItem.setCheckState(currentItem.checkState())

    self.CheckableComboBox_visualization.blockSignals(False)

    # ******* Update the checkboxes in the table ******* #
    for row in range(0, self.tableWidget_visualization.rowCount):
      actionOnCheckBox = False
      label = self.tableWidget_visualization.cellWidget(row, 0)
      outputBasename = label.text
      if currentText == "All Models":
        actionOnCheckBox = True

      elif currentText == "All SPHARM Models":
        if not outputBasename.find("SPHARM") == -1 and outputBasename.find("SPHARM_ellalign") == -1 and outputBasename.find("SPHARMMedialMesh") == -1 and outputBasename.find("SPHARM_procalign") == -1:
          actionOnCheckBox = True

      elif currentText == "All SPHARM Ellipse Aligned Models":
        if not outputBasename.find("SPHARM_ellalign") == -1:
          actionOnCheckBox = True

      elif currentText == "All SPHARM Medial Meshes":
        if not outputBasename.find("SPHARMMedialMesh") == -1:
          actionOnCheckBox = True

      elif currentText == "All SPHARM Procrustes Affine Models":
        if not outputBasename.find("SPHARM_procalign") == -1:
          actionOnCheckBox = True

      else:
        for inputFilename in self.Logic.InputCases:
          inputBasename = inputFilename.split('/')[-1].split('.')[0]
          if not currentText.find(inputBasename) == -1:
            if not currentText.find("SPHARM Models") == -1:
              if not outputBasename.find(inputBasename) == -1 and not outputBasename.find("SPHARM") == -1 and outputBasename.find("SPHARM_ellalign") == -1 and outputBasename.find("SPHARMMedialMesh") == -1 and outputBasename.find("SPHARM_procalign") == -1:
                actionOnCheckBox = True
            elif not currentText.find("SPHARM Ellipse Aligned Models") == -1:
              if not outputBasename.find(inputBasename) == -1 and not outputBasename.find("SPHARM_ellalign") == -1:
                actionOnCheckBox = True
            elif not currentText.find("SPHARM Medial Meshes") == -1:
              if not outputBasename.find(inputBasename) == -1 and not outputBasename.find("SPHARMMedialMesh") == -1:
                actionOnCheckBox = True
            elif not currentText.find("SPHARM Procrustes Affine Models") == -1:
              if not outputBasename.find(inputBasename) == -1 and not outputBasename.find("SPHARM_procalign") == -1:
                actionOnCheckBox = True

      # check/uncheck the checkBox in (row,1)
      if actionOnCheckBox:
          widget = self.tableWidget_visualization.cellWidget(row, 1)
          tuple = widget.children()
          checkBox = tuple[1]
          checkBox.blockSignals(True)
          item = self.CheckableComboBox_visualization.model().item(currentIndex, 0)
          if item.checkState():
            checkBox.setChecked(True)
          else:
            checkBox.setChecked(False)
          checkBox.blockSignals(False)

  # Function to manage the checkbox in the table of visualization
  def onCheckBoxTableValueChanged(self):
    self.CheckableComboBox_visualization.blockSignals(True)
    list = self.CheckableComboBox_visualization.model()
    table = self.tableWidget_visualization

    allIndex = self.CheckableComboBox_visualization.findText("All Models")
    allItem = list.item(allIndex, 0)
    allSPHARMIndex = self.CheckableComboBox_visualization.findText("All SPHARM Models")
    allSPHARMItem = list.item(allSPHARMIndex, 0)
    allSPHARMEllalignIndex = self.CheckableComboBox_visualization.findText("All SPHARM Ellipse Aligned Models")
    allSPHARMEllalignItem = list.item(allSPHARMEllalignIndex, 0)
    allSPHARMMesdialMeshesIndex = self.CheckableComboBox_visualization.findText("All SPHARM Medial Meshes")
    if not allSPHARMMesdialMeshesIndex == -1:
      allSPHARMMesdialMeshesItem = list.item(allSPHARMMesdialMeshesIndex, 0)
    allSPHARMProcrustesAffineModelsIndex = self.CheckableComboBox_visualization.findText("All SPHARM Procrustes Affine Models")
    if not allSPHARMProcrustesAffineModelsIndex == -1:
      allSPHARMProcrustesAffineModelsItem = list.item(allSPHARMProcrustesAffineModelsIndex, 0)

    allSPHARMModelsChecked = True
    allSPHARMEllalignModelsChecked = True
    allSPHARMMedialMeshesChecked = True
    allSPHARMProcrustesAffineModelsChecked = True

    for i in range(len(self.Logic.InputCases)):
      allCaseSPHARMModelsChecked = True
      allCaseSPHARMEllalignModelsChecked = True
      allCaseSPHARMMedialMeshesChecked = True
      allCaseSPHARMProcrustesAffineModelsChecked = True

      inputBasename = self.Logic.InputCases[i].split('/')[-1].split('.')[0]
      for row in range(0,table.rowCount):
        label = table.cellWidget(row, 0)
        outputBasename = label.text
        if not outputBasename.find(inputBasename) == -1:
          widget = table.cellWidget(row, 1)
          tuple = widget.children()
          checkBox = tuple[1]
          if not checkBox.checkState():
            if not outputBasename.find("SPHARM") == -1 and outputBasename.find("SPHARM_ellalign") == -1 and outputBasename.find("SPHARMMedialMesh") == -1 and outputBasename.find("SPHARM_procalign") == -1:
              allCaseSPHARMModelsChecked = False
              allSPHARMModelsChecked = False
              text = "Case " + str(i) + ": " + inputBasename + " - SPHARM Models"
              item = list.findItems(text)[0]
              if item.checkState():
                item.setCheckState(qt.Qt.Unchecked)
              if allSPHARMItem.checkState():
                allSPHARMItem.setCheckState(qt.Qt.Unchecked)
                allItem.setCheckState(qt.Qt.Unchecked)

            if not outputBasename.find("SPHARM_ellalign") == -1:
              allCaseSPHARMEllalignModelsChecked = False
              allSPHARMEllalignModelsChecked = False
              text = "Case " + str(i) + ": " + inputBasename + " - SPHARM Ellipse Aligned Models"
              item = list.findItems(text)[0]
              if item.checkState():
                item.setCheckState(qt.Qt.Unchecked)
              if allSPHARMEllalignItem.checkState():
                allSPHARMEllalignItem.setCheckState(qt.Qt.Unchecked)
                allItem.setCheckState(qt.Qt.Unchecked)

            if not allSPHARMMesdialMeshesIndex == -1:
              if not outputBasename.find("SPHARMMedialMesh") == -1:
                allCaseSPHARMMedialMeshesChecked = False
                allSPHARMMedialMeshesChecked = False
                text = "Case " + str(i) + ": " + inputBasename + " - SPHARM Medial Meshes"
                item = list.findItems(text)[0]
                if item.checkState():
                  item.setCheckState(qt.Qt.Unchecked)
                if allSPHARMMesdialMeshesItem.checkState():
                  allSPHARMMesdialMeshesItem.setCheckState(qt.Qt.Unchecked)
                  allItem.setCheckState(qt.Qt.Unchecked)

            if not allSPHARMProcrustesAffineModelsIndex == -1:
              if not outputBasename.find("SPHARM_procalign") == -1:
                allCaseSPHARMProcrustesAffineModelsChecked = False
                allSPHARMProcrustesAffineModelsChecked = False
                text = "Case " + str(i) + ": " + inputBasename + " - SPHARM Procrustes Affine Models"
                item = list.findItems(text)[0]
                if item.checkState():
                  item.setCheckState(qt.Qt.Unchecked)
                if allSPHARMProcrustesAffineModelsItem.checkState():
                  allSPHARMProcrustesAffineModelsItem.setCheckState(qt.Qt.Unchecked)
                  allItem.setCheckState(qt.Qt.Unchecked)

      if allCaseSPHARMModelsChecked:
        text = "Case " + str(i) + ": " + inputBasename + " - SPHARM Models"
        item = list.findItems(text)[0]
        if not item.checkState():
          item.setCheckState(qt.Qt.Checked)
      if allCaseSPHARMEllalignModelsChecked:
        text = "Case " + str(i) + ": " + inputBasename + " - SPHARM Ellipse Aligned Models"
        item = list.findItems(text)[0]
        if not item.checkState():
          item.setCheckState(qt.Qt.Checked)
      if not allSPHARMMesdialMeshesIndex == -1:
        if allCaseSPHARMMedialMeshesChecked:
          text = "Case " + str(i) + ": " + inputBasename + " - SPHARM Medial Meshes"
          item = list.findItems(text)[0]
          if not item.checkState():
            item.setCheckState(qt.Qt.Checked)
      if not allSPHARMProcrustesAffineModelsIndex == -1:
        if allCaseSPHARMProcrustesAffineModelsChecked:
          text = "Case " + str(i) + ": " + inputBasename + " - SPHARM Procrustes Affine Models"
          item = list.findItems(text)[0]
          if not item.checkState():
            item.setCheckState(qt.Qt.Checked)

    if allSPHARMModelsChecked:
      if not allSPHARMItem.checkState():
        allSPHARMItem.setCheckState(qt.Qt.Checked)
    if allSPHARMEllalignModelsChecked:
      if not allSPHARMEllalignItem.checkState():
        allSPHARMEllalignItem.setCheckState(qt.Qt.Checked)
    if not allSPHARMMesdialMeshesIndex == -1:
      if allSPHARMMedialMeshesChecked:
        if not allSPHARMMesdialMeshesItem.checkState():
          allSPHARMMesdialMeshesItem.setCheckState(qt.Qt.Checked)
    if not allSPHARMProcrustesAffineModelsIndex == -1:
      if allSPHARMProcrustesAffineModelsChecked:
        if not allSPHARMProcrustesAffineModelsItem.checkState():
          allSPHARMProcrustesAffineModelsItem.setCheckState(qt.Qt.Checked)

    if allSPHARMItem.checkState() and allSPHARMEllalignItem.checkState():
      if allSPHARMMesdialMeshesIndex == -1 and allSPHARMProcrustesAffineModelsIndex == -1:
        allItem.setCheckState(qt.Qt.Checked)
      elif not allSPHARMMesdialMeshesIndex == -1 and not allSPHARMProcrustesAffineModelsIndex == -1:
        if allSPHARMMesdialMeshesItem.checkState() and allSPHARMProcrustesAffineModelsItem.checkState():
          allItem.setCheckState(qt.Qt.Checked)
      elif not allSPHARMMesdialMeshesIndex == -1 and allSPHARMProcrustesAffineModelsIndex == -1:
        if allSPHARMMesdialMeshesItem.checkState():
          allItem.setCheckState(qt.Qt.Checked)
      elif allSPHARMMesdialMeshesIndex == -1 and not allSPHARMProcrustesAffineModelsIndex == -1:
        if allSPHARMProcrustesAffineModelsItem.checkState():
          allItem.setCheckState(qt.Qt.Checked)

    self.CheckableComboBox_visualization.blockSignals(False)

  def onPreviewFlips(self):
    # Creation of a CSV file to load the vtk files in ShapePopulationViewer
    filePathCSV = slicer.app.temporaryPath + '/' + 'PreviewForVisualizationInSPV.csv'
    self.Logic.creationCSVFileForSPV(filePathCSV)

    # Launch the CLI ShapePopulationViewer
    parameters = {}
    parameters["CSVFile"] = filePathCSV
    launcherSPV = slicer.modules.launcher
    slicer.cli.run(launcherSPV, None, parameters, wait_for_completion=True)

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
    self.ProgressBar.setNameVisibility(slicer.qSlicerCLIProgressBar.AlwaysVisible)
    self.ErrorMessage = 'Unexpected error'

  def ShapeAnalysisCases(self):

    # No cases
    if not len(self.InputCases) > 0:
      inputDirectory = self.interface.GroupProjectInputDirectory.directory.encode('utf-8')
      self.ErrorMessage = "No cases found in " + inputDirectory
      self.Node.SetStatus(self.Node.CompletedWithErrors)
      return -1

    # Create pipelines
    else:
      logging.info('%d case(s) found', len(self.InputCases))
      # Init
      for i in range(len(self.InputCases)):
        self.completed[i] = False
        self.pipeline[i] = ShapeAnalysisModulePipeline(i, self.InputCases[i], self.interface)

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

  def Cancel(self):
    self.Node.SetStatus(self.Node.Cancelling)
    for i in range(len(self.pipeline)):
      self.pipeline[i].Cancel()

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
        self.ErrorMessage = current_pipeline.ErrorMessage
        logging.error(current_pipeline.ErrorMessage)
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
          self.configurationVisualization()


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

  # Empty the output folder if the overwrite option is checked
  def cleanOutputFolders(self):
    outputDirectory = self.interface.GroupProjectOutputDirectory.directory.encode('utf-8')

    if self.interface.OverwriteSegPostProcess.checkState():
      PostProcessDirectory = outputDirectory + "/PostProcess"
      if os.path.exists(PostProcessDirectory):
        for filename in os.listdir(PostProcessDirectory):
          os.remove(os.path.join(PostProcessDirectory, filename))

    if self.interface.OverwriteGenParaMesh.checkState():
      GenParaMeshOutputDirectory = outputDirectory + "/MeshParameters"
      if os.path.exists(GenParaMeshOutputDirectory):
        for filename in os.listdir(GenParaMeshOutputDirectory):
          os.remove(os.path.join(GenParaMeshOutputDirectory, filename))

    if self.interface.OverwriteParaToSPHARMMesh.checkState():
      SPHARMMeshOutputDirectory = outputDirectory + "/SPHARMMesh"
      if os.path.exists(SPHARMMeshOutputDirectory):
        for filename in os.listdir(SPHARMMeshOutputDirectory):
          os.remove(os.path.join(SPHARMMeshOutputDirectory, filename))

  # Function to fill the table of the flip options for all the SPHARM mesh output
  #    - Column 0: filename of the SPHARM mesh output vtk file
  #    - Column 1: combobox with the flip corresponding to the output file
  def fillTableForFlipOptions(self):
    table = self.interface.tableWidget_ChoiceOfFlip
    row = 0

    for basename in self.InputCases:
      table.setRowCount(row + 1)
      # Column 0:
      filename = basename.split('/')[-1].split('.')[0]
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
                         'Flip Along Axis of z',
                         'All'])
      comboBox.setCurrentIndex(self.interface.choiceOfFlip.currentIndex)
      layout.addWidget(comboBox)
      layout.setAlignment(0x84)
      layout.setContentsMargins(0, 0, 0, 0)
      widget.setLayout(layout)
      table.setCellWidget(row, 1, widget)

      row = row + 1

  def configurationVisualization(self):
    # Configuration of the checkable comboBox
    checkableComboBox = self.interface.CheckableComboBox_visualization
    #   clean the checkable comboBox
    list = checkableComboBox.model()
    list.clear()
    checkableComboBox.addItem("All Models")
    checkableComboBox.addItem("All SPHARM Models")
    checkableComboBox.addItem("All SPHARM Ellipse Aligned Models")
    if self.interface.medialMesh.checkState():
      checkableComboBox.addItem("All SPHARM Medial Meshes")
    if self.interface.useRegTemplate.checkState():
      checkableComboBox.addItem("All SPHARM Procrustes Affine Models")
    #   Fill the checkable comboBox
    for i in range(len(self.InputCases)):
      checkableComboBox.addItem("Case " + str(i) + ": " + self.InputCases[i].split('/')[-1].split('.')[0] + " - SPHARM Models")
      checkableComboBox.addItem("Case " + str(i) + ": " + self.InputCases[i].split('/')[-1].split('.')[0] + " - SPHARM Ellipse Aligned Models")
      if self.interface.medialMesh.checkState():
        checkableComboBox.addItem("Case " + str(i) + ": " + self.InputCases[i].split('/')[-1].split('.')[0] + " - SPHARM Medial Meshes")
      if self.interface.useRegTemplate.checkState():
        checkableComboBox.addItem("Case " + str(i) + ": " + self.InputCases[i].split('/')[-1].split('.')[0] + " - SPHARM Procrustes Affine Models")

    # Configuration of the table
    table = self.interface.tableWidget_visualization
    outputDirectory = self.interface.GroupProjectOutputDirectory.directory.encode('utf-8')
    SPHARMMeshOutputDirectory = outputDirectory + "/SPHARMMesh/"
    row = 0

    for filename in os.listdir(SPHARMMeshOutputDirectory):
      if filename.endswith(".vtk") and not filename.endswith("_para.vtk") and not filename.endswith("SPHARMMedialAxis.vtk"):
        table.setRowCount(row + 1)
        # Column 0:
        labelVTKFile = qt.QLabel(os.path.splitext(filename)[0])
        labelVTKFile.setAlignment(0x84)
        table.setCellWidget(row, 0, labelVTKFile)

        # Column 1:
        widget = qt.QWidget()
        layout = qt.QHBoxLayout(widget)
        checkBox = qt.QCheckBox()
        layout.addWidget(checkBox)
        layout.setAlignment(0x84)
        layout.setContentsMargins(0, 0, 0, 0)
        widget.setLayout(layout)
        table.setCellWidget(row, 1, widget)
        checkBox.connect('stateChanged(int)', self.interface.onCheckBoxTableValueChanged)

        row = row + 1

  # Function to create a CSV file containing all the SPHARM mesh output files
  # that the user wants to display in ShapePopultaionViewer in order to check the flip
  def creationCSVFileForSPV(self, filepathCSV):
    table = self.interface.tableWidget_visualization

    # Creation a CSV file with a header 'VTK Files'
    file = open(filepathCSV, 'w')
    cw = csv.writer(file, delimiter=',')
    cw.writerow(['VTK Files'])

    # Add the filepath of the vtk file checked in the table
    outputDirectory = self.interface.GroupProjectOutputDirectory.directory.encode('utf-8')
    SPHARMMeshOutputDirectory = outputDirectory + "/SPHARMMesh/"

    # Add the path of the vtk files if the users selected it
    for row in range(0, table.rowCount):
      # check the checkBox
      widget = table.cellWidget(row, 1)
      tuple = widget.children()
      checkBox = tuple[1]
      if checkBox.isChecked():
        # Recovery of the vtk filename
        qlabel = table.cellWidget(row, 0)
        vtkBasename = qlabel.text
        VTKfilepath = SPHARMMeshOutputDirectory + vtkBasename + ".vtk"
        if os.path.exists(VTKfilepath):
          cw.writerow([VTKfilepath])
    file.close()

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
  def __init__(self, pipelineID, CaseInput, interface):

    VTKObservationMixin.__init__(self)

    self.interface = interface

    self.pipelineID = pipelineID
    self.strPipelineID = "_" + str(self.pipelineID)
    self.CaseInput = CaseInput

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

  def setupGlobalVariables(self):
    # Modules
    self.ID = -1
    self.slicerModule = {}
    self.moduleParameters = {}

    # Nodes
    self.nodeDictionary = {}

    # Input filename and extension
    filepathSplit = self.CaseInput.split('/')[-1].split('.')
    self.inputFilename = filepathSplit[0]
    self.inputExtension = filepathSplit[1]
    if len(filepathSplit) == 3:
      self.inputExtension = self.inputExtension + "." + filepathSplit[2]

  def setupSkipCLIs(self):

    self.skip_meshToLabelMap = False
    self.skip_segPostProcess = False
    self.skip_genParaMesh = False
    self.skip_paraToSPHARMMesh = False
    outputDirectory = self.interface.GroupProjectOutputDirectory.directory.encode('utf-8')

    # Skip MeshToLabelMap?
    if not self.inputExtension == "vtk" and not self.inputExtension == "vtp":
      self.skip_meshToLabelMap = True
    else:
      LabelMapDirectory = outputDirectory + "/LabelMap"
      LabelMapOutputFilepath = LabelMapDirectory + "/" + self.inputFilename + ".nrrd"
      if os.path.exists(LabelMapOutputFilepath):
        self.inputExtension = "nrrd"
        self.skip_meshToLabelMap = True

    # If MeshToLabelMap is not skip, do not skip the next CLIs: SegPostProcess, GenParaMesh and ParaToSPHARMMesh
    if self.skip_meshToLabelMap == False:
      return

    # Skip SegPostProcess ?
    if not self.interface.OverwriteSegPostProcess.checkState():
      PostProcessDirectory = outputDirectory + "/PostProcess"
      PostProcessOutputFilepath = PostProcessDirectory + "/" + self.inputFilename + "_pp." + self.inputExtension
      if os.path.exists(PostProcessOutputFilepath):
        self.skip_segPostProcess = True

    # If SegPostProcess is not skip, do not skip the next CLIs: GenParaMesh and ParaToSPHARMMesh
    if self.skip_segPostProcess == False:
      return

    # Skip GenParaMesh ?
    if not self.interface.OverwriteGenParaMesh.checkState():
      GenParaMeshOutputDirectory = outputDirectory + "/MeshParameters"
      ParaOutputFilepath = GenParaMeshOutputDirectory + "/" + self.inputFilename + "_para.vtk"
      SurfOutputFilepath = GenParaMeshOutputDirectory + "/" + self.inputFilename + "_surf.vtk"
      if os.path.exists(ParaOutputFilepath) and os.path.exists(SurfOutputFilepath):
        self.skip_genParaMesh = True

    # If GenParaMesh is not skip, do not skip the next CLI: ParaToSPHARMMesh
    if self.skip_genParaMesh == False:
      return

    # Skip ParaToSPHARMMesh ?
    if not self.interface.OverwriteParaToSPHARMMesh.checkState():
      SPHARMMeshOutputDirectory = outputDirectory + "/SPHARMMesh"
      SPHARMMeshFilepath = SPHARMMeshOutputDirectory + "/" + self.inputFilename
      SPHARMMeshDirectory = os.path.dirname(SPHARMMeshFilepath)
      SPHARMMeshBasename = os.path.basename(SPHARMMeshFilepath)
      if os.path.exists(SPHARMMeshDirectory):
        for file in os.listdir(SPHARMMeshDirectory):
          if not file.find(SPHARMMeshBasename) == -1:
              self.skip_paraToSPHARMMesh = True

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

    ## Mesh To Label Map: Transform input model in label map
    cli_nodes = list() # list of the nodes used in the Mesh to Label Map step
    cli_filepaths = list() # list of the node filepaths used in the Mesh to Label Map step
    LabelMapDirectory = outputDirectory + "/LabelMap"
    LabelMapOutputFilepath = LabelMapDirectory + "/" + self.inputFilename + ".nrrd"
    if not self.skip_meshToLabelMap:
      # Setup of the parameters od the CLI
      self.ID += 1

      cli_parameters = {}

      inputFilepath = inputDirectory + '/' + self.CaseInput
      model_input_node = ShapeAnalysisModuleMRMLUtility.loadMRMLNode(inputFilepath, 'ModelFile')

      cli_parameters["mesh"] = model_input_node

      meshtolabelmap_output_node = ShapeAnalysisModuleMRMLUtility.addnewMRMLNode("output_MeshToLabelMap", slicer.vtkMRMLLabelMapVolumeNode())
      cli_parameters["labelMap"] = meshtolabelmap_output_node

      cli_parameters["spacingVec"] = "0.1,0.1,0.1"

      self.inputExtension = "nrrd"

      self.setupModule(slicer.modules.meshtolabelmap, cli_parameters)

      # Setup of the nodes created by the CLI
      #    Creation of a folder in the output folder : LabelMap
      if not os.path.exists(LabelMapDirectory):
        os.makedirs(LabelMapDirectory)

      cli_nodes.append(model_input_node)
      cli_nodes.append(meshtolabelmap_output_node)
      cli_filepaths.append(inputFilepath)
      cli_filepaths.append(LabelMapOutputFilepath)

      self.setupNode(0, cli_nodes, cli_filepaths, [False, True], [True, True])
    else:
      if os.path.exists(LabelMapOutputFilepath):
        # Setup of the nodes which will be use by the next CLI
        meshtolabelmap_output_node = ShapeAnalysisModuleMRMLUtility.loadMRMLNode(LabelMapOutputFilepath, 'LabelMapVolumeFile')

        cli_nodes.append(meshtolabelmap_output_node)
        cli_filepaths.append(LabelMapOutputFilepath)

        self.setupNode(0, cli_nodes, cli_filepaths, [False], [True])

    ## Post Processed Segmentation
    cli_nodes = list() # list of the nodes used in the Post Processed Segmentation step
    cli_filepaths = list() # list of the node filepaths used in the Post Processed Segmentation step
    PostProcessDirectory = outputDirectory + "/PostProcess"
    PostProcessOutputFilepath = PostProcessDirectory + "/" + self.inputFilename + "_pp." + self.inputExtension

    if not self.skip_segPostProcess:
      # Setup of the parameters of the CLI
      self.ID += 1

      cli_parameters = {}

      if self.skip_meshToLabelMap: # IF the input given was already a label map
        inputFilepath = inputDirectory + '/' + self.CaseInput
        labelmap_input_node = ShapeAnalysisModuleMRMLUtility.loadMRMLNode(inputFilepath, 'LabelMapVolumeFile')
      else:                        # ELSE the input given was a model which has been transformed by MeshToLabelMap and store in the folder LabelMap
        labelmap_input_node = meshtolabelmap_output_node
        inputFilepath = LabelMapOutputFilepath

      cli_parameters["fileName"] = labelmap_input_node

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

      cli_nodes.append(labelmap_input_node)
      cli_nodes.append(pp_output_node)
      cli_filepaths.append(inputFilepath)
      cli_filepaths.append(PostProcessOutputFilepath)

      self.setupNode(1, cli_nodes, cli_filepaths, [False,True], [True,True])

    else:
      # Setup of the nodes which will be use by the next CLI
      pp_output_node = ShapeAnalysisModuleMRMLUtility.loadMRMLNode(PostProcessOutputFilepath, 'LabelMapVolumeFile')

      cli_nodes.append(pp_output_node)
      cli_filepaths.append(PostProcessOutputFilepath)

      self.setupNode(1, cli_nodes, cli_filepaths, [False], [True])


    ## Generate Mesh Parameters
    cli_nodes = list() # list of the nodes used in the Generate Mesh Parameters step
    cli_filepaths = list() # list of the node filepaths used in the Generate Mesh Parameters step
    GenParaMeshOutputDirectory = outputDirectory + "/MeshParameters"
    ParaOutputFilepath = GenParaMeshOutputDirectory + "/" + self.inputFilename + "_para.vtk"
    SurfOutputFilepath = GenParaMeshOutputDirectory + "/" + self.inputFilename + "_surf.vtk"

    if not self.skip_genParaMesh:
      # Setup of the parameters of the CLI
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
    cli_nodes = list()  # list of the nodes used in the Parameters To SPHARM Mesh step
    cli_filepaths = list()  # list of the node filepaths used in the Parameters To SPHARM Mesh step
    SPHARMMeshOutputDirectory = outputDirectory + "/SPHARMMesh"
    if not self.skip_paraToSPHARMMesh:

      # Search of the flip to apply:
      #       1 = flip along axes of x &amp; y,
      #       2 = flip along y &amp; z,
      #       3 = flip along x &amp; z
      #       4 = flip along x,
      #       5 = flip along y,
      #       6 = flip along x &amp; y &amp; z,
      #       7 = flip along z  where y is the smallest, x is the second smallest and z is the long axis of the ellipsoid
      #       8 = All the flips
      if not self.interface.sameFlipForAll.checkState():
        # Recovery of the flip chosen by the user
        row = self.pipelineID
        widget = self.interface.tableWidget_ChoiceOfFlip.cellWidget(row, 1)
        tuple = widget.children()
        comboBox = qt.QComboBox()
        comboBox = tuple[1]
        flipIndexToApply = comboBox.currentIndex
      else:
        flipIndexToApply = self.interface.choiceOfFlip.currentIndex

      # Only one flip to apply
      if flipIndexToApply < 8:
        L = [1]
      # All the flip to apply
      else:
        L = range(1,8)

      for i in L:

        # Setup of the parameters of the CLI
        self.ID += 1

        cli_parameters = {}

        cli_parameters["inParaFile"] = para_output_model

        cli_parameters["inSurfFile"] = surfmesh_output_model

        #    Creation of a folder in the output folder : SPHARMMesh
        if not os.path.exists(SPHARMMeshOutputDirectory):
          os.makedirs(SPHARMMeshOutputDirectory)
        if flipIndexToApply < 8:
          SPHARMMeshFilepath = SPHARMMeshOutputDirectory + "/" + self.inputFilename
          cli_parameters["outbase"] = SPHARMMeshFilepath
        #   For each flip creation of a output filename
        else:
          flipName = ['AlongXY', 'AlongYZ', 'AlongXZ', 'AlongX', 'AlongY', 'AlongXYZ', 'AlongZ']
          SPHARMMeshFilepath = SPHARMMeshOutputDirectory + "/" + self.inputFilename + "_flip" + flipName[i - 1]
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
        if self.interface.useRegTemplate.checkState():
          cli_parameters["regTemplateFileOn"] = True
          regtemplate_filepath = self.interface.regTemplate.currentPath
          regtemplate_model = ShapeAnalysisModuleMRMLUtility.loadMRMLNode(regtemplate_filepath, 'ModelFile')
          cli_parameters["regTemplateFile"] = regtemplate_model

          cli_nodes.append(regtemplate_model)
          cli_filepaths.append(regtemplate_filepath)

          self.setupNode(i + 2, cli_nodes, cli_filepaths, [False], [True])
        if self.interface.useFlipTemplate.checkState():
          cli_parameters["flipTemplateFileOn"] = True
          cli_parameters["flipTemplateFile"] = self.interface.flipTemplate.currentPath

        if flipIndexToApply < 8 :
            cli_parameters["finalFlipIndex"] = flipIndexToApply
        else:
          cli_parameters["finalFlipIndex"] = i

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
      self.deleteNodes()
      logging.info('Slicer Module queue is empty, nothing to do')
      self.Node.SetStatus(self.Node.Completed)

  def onCLIModuleModified(self, cli_node, event):
    logging.info('-- %s : %s', cli_node.GetStatusString(), cli_node.GetName())
    if not cli_node.IsBusy():
      if platform.system() != 'Windows':
        self.removeObserver(cli_node, self.StatusModifiedEvent, self.onCLIModuleModified)
        statusForNode =  None

      if cli_node.GetStatusString() == 'Completed':
        if self.ID == len(self.slicerModule) - 1:
          cli_node_name = "Case " + str(self.pipelineID) + ": " + self.inputFilename
          cli_node.SetName(cli_node_name)
          self.setCurrentCLINode(cli_node)
          statusForNode = cli_node.GetStatus()

      elif cli_node.GetStatusString() == 'Cancelled':
        self.ErrorMessage = cli_node.GetName() + " Cancelled: cancelling the pipeline"
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

  def createCLINode(self, module):
    cli_node = slicer.cli.createNode(module)
    cli_node_name = "Case " + str(self.pipelineID) + ": " + self.inputFilename + " - step " + str(self.ID) + ": " + module.title
    cli_node.SetName(cli_node_name)
    return cli_node

  def saveNodes(self):
    for id in self.nodeDictionary.keys():
      save = self.nodeDictionary[id].save
      node = self.nodeDictionary[id].nodes
      for i in range(len(self.nodeDictionary[id].save)):
        filepaths = self.nodeDictionary[id].filepaths
        if save[i] == True:
          ShapeAnalysisModuleMRMLUtility.saveMRMLNode( node[i], filepaths[i] )

  def deleteNodes(self):
    for id in self.nodeDictionary.keys():
      delete = self.nodeDictionary[id].delete
      node = self.nodeDictionary[id].nodes
      for i in range(len(self.nodeDictionary[id].delete)):
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