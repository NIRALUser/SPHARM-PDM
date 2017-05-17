import os, sys
import vtk, ctk, slicer

#
# MRMLUtility
#

'''
This class harbors all the utility functions to load, add, save and remove mrml nodes
'''

class MRMLUtility(object):

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