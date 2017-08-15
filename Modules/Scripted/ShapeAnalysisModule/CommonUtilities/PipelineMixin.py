import vtk, ctk, slicer
import logging
import time
import gc
import platform
import os

from slicer.util import VTKObservationMixin
from CommonUtilities import MRMLUtility


#
# ShapeAnalysisModuleNode
#
class ShapeAnalysisModuleNode(object):
  nodes = [None]
  save = [False]
  delete = [True]
  filepaths = [" "]

#
# PipelineMixin
#
'''
This class is a base class for pipelines.
A pipeline class deriving from the Mixin should impletent a setup  function.
'''
class PipelineMixin(VTKObservationMixin):
  def __init__(self, pipelineID, CaseInput, interface):

    VTKObservationMixin.__init__(self)

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
    self.inputRootname = filepathSplit[0]
    self.inputExtension = filepathSplit[1]
    if len(filepathSplit) == 3:
      self.inputExtension = self.inputExtension + "." + filepathSplit[2]

  def setupModule(self, module, cli_parameters):
    self.slicerModule[self.ID] = module
    self.moduleParameters[self.ID] = cli_parameters

  def setupNode(self, id, cli_nodes, cli_filepaths, cli_saveOutput, cli_deleteOutput):
    self.nodeDictionary[id] = ShapeAnalysisModuleNode()
    self.nodeDictionary[id].nodes = cli_nodes
    self.nodeDictionary[id].filepaths = cli_filepaths
    self.nodeDictionary[id].save = cli_saveOutput
    self.nodeDictionary[id].delete = cli_deleteOutput

  def onCLIModuleModified(self, cli_node, event):
    logging.info('-- %s : %s', cli_node.GetStatusString(), cli_node.GetName())
    statusForNode = None
    if not cli_node.IsBusy():
      if platform.system() != 'Windows':
        self.removeObserver(cli_node, self.StatusModifiedEvent, self.onCLIModuleModified)
        statusForNode =  None

      if cli_node.GetStatusString() == 'Completed':
        if self.ID == len(self.slicerModule) - 1:
          cli_node_name = "Case " + str(self.pipelineID) + ": " + self.inputRootname
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

  # Run the CLI for the current module
  def runCLIModule(self):
    # Create and run CLI
    module = self.slicerModule[self.ID]
    cli_node = self.createCLINode(module)
    self.setCurrentCLINode(cli_node) # ProgressBar
    self.addObserver(cli_node, self.StatusModifiedEvent, self.onCLIModuleModified)
    slicer.cli.run(module, cli_node, self.moduleParameters[self.ID], wait_for_completion = self.interface.waitForCompletion)

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

  def saveNodes(self):
    for id in self.nodeDictionary.keys():
      save = self.nodeDictionary[id].save
      node = self.nodeDictionary[id].nodes
      for i in range(len(self.nodeDictionary[id].save)):
        filepaths = self.nodeDictionary[id].filepaths
        if save[i] == True:
          MRMLUtility.saveMRMLNode( node[i], filepaths[i] )

  def deleteNodes(self):
    for id in self.nodeDictionary.keys():
      delete = self.nodeDictionary[id].delete
      node = self.nodeDictionary[id].nodes
      for i in range(len(self.nodeDictionary[id].delete)):
        if delete[i] == True:
          MRMLUtility.removeMRMLNode( node[i] )

  def createCLINode(self, module):
    cli_node = slicer.cli.createNode(module)
    cli_node_name = "Case " + str(self.pipelineID) + ": " + self.inputRootname + " - step " + str(self.ID) + ": " + module.title
    cli_node.SetName(cli_node_name)
    return cli_node

  def setCurrentCLINode(self, cli_node):
    self.ProgressBar.setCommandLineModuleNode(cli_node)
    self.currentCLINode = cli_node

  def Cancel(self):
    if self.currentCLINode:
      self.currentCLINode.Cancel()