import sys, os
from ShapeAnalysisModule import ShapeAnalysisModuleWrapper
from ConfigParser import SafeConfigParser

if len(sys.argv) != 2:
  print("\n\nUsage: ./SlicerSALT --no-main-window --python-script " + sys.argv[0] + " path/to/SPHARM-PDM/parameters/file\n\n")
  sys.exit(1)

if not os.path.exists(sys.argv[1]):
    print "\n\nThe SPHARM-PDM parameters file " + sys.argv[1] + " doesn't exit!\n\n"
    sys.exit(1)

parser = SafeConfigParser()
parser.read(sys.argv[1])

inputDirectoryPath = parser.get('DirectoryPath', 'inputDirectoryPath')
if not os.path.exists(inputDirectoryPath):
    print "/!\ The input directory: " + inputDirectoryPath + " doesn't exit! /!\ "
    sys.exit(1)

outputDirectoryPath = parser.get('DirectoryPath', 'outputDirectoryPath')
if not os.path.exists(outputDirectoryPath):
    print "/!\ The output directory path: " + outputDirectoryPath + " doesn't exit! /!\ "
    sys.exit(1)

if parser.get('SegPostProcess', 'rescale') == 'True' or parser.get('SegPostProcess', 'rescale') == 'true':
    RescaleSegPostProcess = True
elif parser.get('SegPostProcess', 'rescale') == 'False' or parser.get('SegPostProcess', 'rescale') == 'false':
    RescaleSegPostProcess = False
else:
    print "/!\ The paremeter 'rescale' can only take as value True or False /!\ "
    sys.exit(1)

sx = float(parser.get('SegPostProcess', 'space').split(',')[0])
if sx < 0  or sx > 1:
    print "/!\ The paremeter 'space' for the x direction can only be a float between 0 and 1 /!\ "
    sys.exit(1)
sy = float(parser.get('SegPostProcess', 'space').split(',')[1])
if sy < 0  or sy > 1:
    print "/!\ The paremeter 'space' for the y direction can only be a float between 0 and 1 /!\ "
    sys.exit(1)
sz = float(parser.get('SegPostProcess', 'space').split(',')[2])
if sz < 0  or sz > 1:
    print "/!\ The paremeter 'space' for the z direction can only be a float between 0 and 1 /!\ "
    sys.exit(1)

if int(parser.get('SegPostProcess', 'label')) < 0 or int(parser.get('SegPostProcess', 'label')) > 100:
    print "/!\ The paremeter 'label' can only be a integer between 0 and 100 /!\ "
    sys.exit(1)
labelNumber = int(parser.get('SegPostProcess', 'label'))

if parser.get('SegPostProcess', 'gauss') == 'True' or parser.get('SegPostProcess', 'gauss') == 'true':
    GaussianFiltering = True
elif parser.get('SegPostProcess', 'gauss') == 'False' or parser.get('SegPostProcess', 'gauss') == 'false':
    GaussianFiltering = False
else:
    print "/!\ The paremeter 'gauss' can only take as value True or False /!\ "
    sys.exit(1)

VarianceX = int(parser.get('SegPostProcess', 'var').split(',')[0])
if VarianceX < 0  or VarianceX > 100:
    print "/!\ The paremeter 'var' for the x direction can only be a integer between 0 and 100 /!\ "
    sys.exit(1)
VarianceY = int(parser.get('SegPostProcess', 'var').split(',')[1])
if VarianceY < 0  or VarianceY > 100:
    print "/!\ The paremeter 'var' for the y direction can only be a integer between 0 and 100 /!\ "
    sys.exit(1)
VarianceZ = int(parser.get('SegPostProcess', 'var').split(',')[2])
if VarianceZ < 0  or VarianceZ > 100:
    print "/!\ The paremeter 'var' for the z direction can only be a integer between 0 and 100 /!\ "
    sys.exit(1)

if int(parser.get('GenParaMesh', 'iter')) < 0 or int(parser.get('GenParaMesh', 'iter')) > 10000:
    print "/!\ The paremeter 'iter' can only be a integer between 0 and 10 000 /!\ "
    sys.exit(1)
numberofIterations = int(parser.get('GenParaMesh', 'iter'))

if int(parser.get('ParaToSPHARMMesh', 'subdivLevel')) < 0 or int(parser.get('ParaToSPHARMMesh', 'subdivLevel')) > 100:
    print "/!\ The paremeter 'subdivLevel' can only be a integer between 0 and 100 /!\ "
    sys.exit(1)
SubdivLevelValue = int(parser.get('ParaToSPHARMMesh', 'subdivLevel'))

if int(parser.get('ParaToSPHARMMesh', 'spharmDegree')) < 0 or int(parser.get('ParaToSPHARMMesh', 'spharmDegree')) > 100:
    print "/!\ The paremeter 'spharmDegree' can only be a integer between 0 and 100 /!\ "
    sys.exit(1)
SPHARMDegreeValue = int(parser.get('ParaToSPHARMMesh', 'spharmDegree'))

if parser.get('ParaToSPHARMMesh', 'medialMesh') == 'True' or parser.get('ParaToSPHARMMesh', 'medialMesh') == 'true':
    medialMesh = True
elif parser.get('ParaToSPHARMMesh', 'medialMesh') == 'False' or parser.get('ParaToSPHARMMesh', 'medialMesh') == 'false':
    medialMesh = False
else:
    print "/!\ The paremeter 'gauss' can only take as value True or False /!\ "
    sys.exit(1)

if int(parser.get('ParaToSPHARMMesh', 'thetaIteration')) < 0 or int(parser.get('ParaToSPHARMMesh', 'thetaIteration')) > 21474836:
    print "/!\ The paremeter 'spharmDegree' can only be a integer between 0 and 21474836 /!\ "
    sys.exit(1)
thetaIterationValue = parser.get('ParaToSPHARMMesh', 'thetaIteration')

if int(parser.get('ParaToSPHARMMesh', 'phiIteration')) < 0 or int(parser.get('ParaToSPHARMMesh', 'phiIteration')) > 21474836:
    print "/!\ The paremeter 'spharmDegree' can only be a integer between 0 and 21474836 /!\ "
    sys.exit(1)
phiIterationValue = int(parser.get('ParaToSPHARMMesh', 'phiIteration'))

if parser.get('ParaToSPHARMMesh', 'regParaTemplateFileOn') == 'True' or parser.get('ParaToSPHARMMesh', 'regParaTemplateFileOn') == 'true':
    useRegTemplate = True
    regTemplate = parser.get('ParaToSPHARMMesh', 'regParaTemplate')
    if not os.path.exists(regTemplate):
        print "/!\ The registration template: " + regTemplate + " doesn't exit! /!\ "
        sys.exit(1)
    if not regTemplate.endswith('.vtk') and not regTemplate.endswith('.vtp'):
        print "/!\ The registration template: " + regTemplate + " is not a VTK file or VTP file! /!\ "
        sys.exit(1)
elif parser.get('ParaToSPHARMMesh', 'regParaTemplateFileOn') == 'False' or parser.get('ParaToSPHARMMesh', 'regParaTemplateFileOn') == 'false':
    useRegTemplate = False
    regTemplate = " "
else:
    print "/!\ The paremeter 'regParaTemplateFileOn' can only take as value True or False /!\ "
    sys.exit(1)

if parser.get('ParaToSPHARMMesh', 'flipTemplateOn') == 'True' or parser.get('ParaToSPHARMMesh', 'flipTemplateOn') == 'true':
    useFlipTemplate = True
    flipTemplate = parser.get('ParaToSPHARMMesh', 'flipTemplate')
    if not os.path.exists(flipTemplate):
        print "/!\ The flip template: " + flipTemplate + " doesn't exit! /!\ "
        sys.exit(1)
    if not flipTemplate.endswith('.coef'):
        print "/!\ The flip template: " + regTemplate + " is not a COEF file! /!\ "
        sys.exit(1)
    flipTemplate = parser.get('ParaToSPHARMMesh', 'flipTemplate')
elif parser.get('ParaToSPHARMMesh', 'flipTemplateOn') == 'False' or parser.get('ParaToSPHARMMesh', 'flipTemplateOn') == 'false':
    useFlipTemplate = False
    flipTemplate = " "
else:
    print "/!\ The paremeter 'flipTemplateOn' can only take as value True or False /!\ "
    sys.exit(1)

if int(parser.get('ParaToSPHARMMesh', 'flip')) < 0 or int(parser.get('ParaToSPHARMMesh', 'flip')) > 8:
    print "/!\ The paremeter 'flip' can only be a integer between 0 and 8 /!\ "
    sys.exit(1)
choiceOfFlip = int(parser.get('ParaToSPHARMMesh', 'flip'))

ShapeAnalysisModuleInstance = ShapeAnalysisModuleWrapper(inputDirectoryPath, outputDirectoryPath,
                                                         RescaleSegPostProcess, sx, sy, sz, labelNumber,
                                                         GaussianFiltering, VarianceX, VarianceY, VarianceZ,
                                                         numberofIterations,
                                                         SubdivLevelValue, SPHARMDegreeValue,
                                                         medialMesh, thetaIterationValue, phiIterationValue,
                                                         useRegTemplate, regTemplate,
                                                         useFlipTemplate, flipTemplate, choiceOfFlip)
ShapeAnalysisModuleInstance.startProcessing()
