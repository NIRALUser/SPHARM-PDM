import vtk, ctk, slicer
import logging
import time
import gc

from slicer.util import VTKObservationMixin
from slicer.ScriptedLoadableModule import *


class LogicMixin(ScriptedLoadableModuleLogic, VTKObservationMixin):
  def __init__(self, ModuleNodeName):
    VTKObservationMixin.__init__(self)

    self.InputCases = list()
    self.allCaseStartTime = 0

    # Dictionaries
    self.pipeline = {}
    self.completed = {}

    # Status
    self.Node = slicer.vtkMRMLCommandLineModuleNode()
    self.Node.SetStatus(self.Node.Idle)
    self.Node.SetName(ModuleNodeName)
    self.ProgressBar = slicer.qSlicerCLIProgressBar()
    self.ProgressBar.setCommandLineModuleNode(self.Node)
    self.ProgressBar.setNameVisibility(slicer.qSlicerCLIProgressBar.AlwaysVisible)
    self.ErrorMessage = 'Unexpected error'

  def startPipeline(self, id):
    self.pipeline[id].setup()
    self.pipeline[id].Node.SetStatus(self.Node.Scheduled)
    self.pipeline[id].runFirstCLIModule()

  def onPipelineModified(self, pipeline_node, event):

    pipeline_id = None
    current_pipeline = None
    for key, pipeline in self.pipeline.items():
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