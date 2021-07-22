import os
import gemmi
from privateer import privateer_core as pvt

# Method to add LINK record into input files: 
# https://project-gemmi.github.io/python-api/gemmi.Structure.html 
# https://project-gemmi.github.io/python-api/gemmi.Structure.html#connections
# https://project-gemmi.github.io/python-api/gemmi.Connection.html
# Combine this method with privateer.privateer_core.CarbohydrateStructure.get_sugar_linkage_info()

inputPath = '/home/harold/Dev/privateer_python/project_alliance/glycampdbfiles/confConvertedPDB/'
outputPath = '/home/harold/Dev/privateer_python/project_alliance/glycampdbfiles/confConvertedmmCIF/'

singular_mmCIF_output = gemmi.cif.Document()
for root, dirs, files in os.walk(inputPath, topdown=False):
  for name in files:
    head, tail = os.path.split(root)
    outputroot = os.path.join(outputPath, tail)
    
    if not os.path.exists(outputroot):
      os.makedirs(outputroot)
      
    inputFilePath = os.path.join(root, name)
    print(inputFilePath)
    
    clusterName = name.replace(".pdb", "")
    name_of_file = name.replace(".pdb", ".mmCIF")
    outputFilePath = os.path.join(outputroot, name_of_file)
    
    per_input_mmCIF = gemmi.cif.Document()
    glycosylation = pvt.GlycosylationComposition(inputFilePath)
    inputGlycan = glycosylation.get_glycan(0)
    glycanWURCS = inputGlycan.get_wurcs_notation()
    gemmiStructure = gemmi.read_structure(inputFilePath)
    gemmiDocument = gemmiStructure.make_mmcif_document()
    gemmiBlock = gemmiDocument.sole_block()
    gemmiBlock.set_pair("_wurcs", glycanWURCS)
        
    # outputTail = tail + clusterName
    # singular_mmCIF_output.add_new_block(outputTail)
    # singular_output_block = singular_mmCIF_output.find_block(outputTail)
    
    # for item in gemmiBlock:
    #   singular_output_block.add_item(item)
    
    with open(outputFilePath, mode="w") as newfile:
      gemmiDocument.write_file(outputFilePath)

# with open("/home/harold/Dev/privateer_python/project_alliance/single.mmCIF", mode="w") as newfile:
#   singular_mmCIF_output.write_file("/home/harold/Dev/privateer_python/project_alliance/single.mmCIF")

# mmcif = structure.make_mmcif_document()
# mmcif.write_file("Cluster1.mmCIF")
