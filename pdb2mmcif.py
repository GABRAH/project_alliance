import os
import shutil
import gemmi
from privateer import (
    privateer_core as pvt,
)  # https://github.com/glycojones/privateer/tree/privateerpython

# Utility functions
def search_list_of_dicts(key, value, list):
    return [element for element in list if element[key] == value][0]


def CreateFolder(path):
    if not os.path.exists(path):
        os.makedirs(path)
    else:
        shutil.rmtree(path)  # Removes all the subdirectories!
        os.makedirs(path)


def getResidueSummaryInChain(chain):
    summary = []

    chainLength = len(chain)

    for i in range(chainLength):
        residue = chain[i]
        currentSeqID = residue.seqid.num
        residueName = residue.name

        gemmiResidueInfo = {
            "index": i,
            "seqID": currentSeqID,
            "residueName": residueName,
        }

        summary.append(gemmiResidueInfo)

    return summary


# Privateer/gemmi functions
def getMetadataFromPrivateer(inputFilePath, privateerJSON):
    glycosylation = pvt.GlycosylationComposition(inputFilePath)
    inputGlycan = glycosylation.get_glycan(0)
    # glycanWURCS = inputGlycan.get_wurcs_notation()
    glycomics = inputGlycan.query_offline_database(privateerJSON, False, False)
    glycanWURCS = glycomics["wurcs"]
    glytoucanID = glycomics["glytoucan_id"]

    totalNumberOfSugars = inputGlycan.get_total_number_of_sugars()

    linkages = []
    for i in range(totalNumberOfSugars):
        currentSugar = inputGlycan.get_monosaccharide(i)
        # currentSugarPDBID = currentSugar.get_sugar_pdb_id()
        currentSugarIndexInVector = i
        currentSugarPDBID = int(currentSugar.get_sugar_pdb_id())
        currentSugarLinkageInfo = currentSugar.get_sugar_linkage_info()
        outputLinkageInfo = {
            "index": currentSugarIndexInVector,
            "pdb_id": currentSugarPDBID,
            "linkage_info": currentSugarLinkageInfo,
        }

        linkages.append(outputLinkageInfo)

    output = {
        "glycanWURCS": glycanWURCS,
        "glytoucanID": glytoucanID,
        "sugar_connections": linkages,
    }

    return output


def addGemmiConnectionsBetweenSugars(gemmiStructure, privateerMetaData):
    outputGemmiStructure = gemmiStructure
    sugarConnectionMetaData = privateerMetaData["sugar_connections"]

    model = outputGemmiStructure[0]
    chain = model[0]
    gemmiResidueInfo = getResidueSummaryInChain(chain)
    chainLength = len(chain)
    for i in range(chainLength):
        residue = chain[i]
        currentPDBID = residue.seqid.num
        currentResidueName = residue.name

        privateerConnectionsSummary = search_list_of_dicts(
            "pdb_id", currentPDBID, sugarConnectionMetaData
        )

        currentResidueLinkageInfo = privateerConnectionsSummary["linkage_info"]
        for connection in currentResidueLinkageInfo:
            connectedToSugarIndex = connection["connectedToSugarID"]
            connectedToSugarIndexPrivateerList = search_list_of_dicts(
                "index", connectedToSugarIndex, sugarConnectionMetaData
            )
            connectedToPDBID = connectedToSugarIndexPrivateerList["pdb_id"]
            connectedToInfoSummary = search_list_of_dicts(
                "seqID", connectedToPDBID, gemmiResidueInfo
            )
            connectedToResidueName = connectedToInfoSummary["residueName"]

            currentResidueAtomName = "C" + connection["hostLinkagePosition"]
            connectedToResidueAtomName = "C" + connection["linkagePositionForeign"]

            new_connection = gemmi.Connection()
            new_connection.name = (
                str(currentPDBID)
                + currentResidueName
                + "_"
                + currentResidueAtomName
                + "-"
                + connectedToResidueAtomName
                + "_"
                + connectedToResidueName
                + str(connectedToPDBID)
            )
            new_connection.type = gemmi.ConnectionType.Covale
            new_connection.asu = gemmi.Asu.Same
            residue_Alpha = chain[str(currentPDBID)][currentResidueName]
            residue_Bravo = chain[str(connectedToPDBID)][connectedToResidueName]

            new_connection.partner1 = gemmi.make_address(
                chain, residue_Alpha, residue_Alpha.sole_atom(currentResidueAtomName)
            )
            new_connection.partner2 = gemmi.make_address(
                chain,
                residue_Bravo,
                residue_Bravo.sole_atom(connectedToResidueAtomName),
            )
            outputGemmiStructure.connections.append(new_connection)

    return outputGemmiStructure


def convertSinglePDBtoCIF(inputFilePath, outputFilePath, privateerJSON):
    privateerMetaData = getMetadataFromPrivateer(inputFilePath, privateerJSON)
    gemmiStructure = gemmi.read_structure(inputFilePath)
    gemmiStructureWithSugarLinks = addGemmiConnectionsBetweenSugars(
        gemmiStructure, privateerMetaData
    )

    gemmiDocument = gemmiStructureWithSugarLinks.make_mmcif_document()
    gemmiBlock = gemmiDocument.sole_block()
    gemmiBlock.set_pair("_WURCS", gemmi.cif.quote(privateerMetaData["glycanWURCS"]))
    gemmiBlock.set_pair("_GlyTouCan", gemmi.cif.quote(privateerMetaData["glytoucanID"]))
    gemmiDocument.write_file(outputFilePath)


inputPath = "/home/harold/Dev/privateer_python/project_alliance/glycampdbfiles/VolumeConvertedPDB/"
outputPath = "/home/harold/Dev/privateer_python/project_alliance/glycampdbfiles/VolumeConvertedmmCIF/"
privateerJSON = pvt.OfflineDatabase()

singular_mmCIF_output = gemmi.cif.Document()
CreateFolder(outputPath)
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

        convertSinglePDBtoCIF(inputFilePath, outputFilePath, privateerJSON)

        # outputTail = tail + clusterName
        # singular_mmCIF_output.add_new_block(outputTail)
        # singular_output_block = singular_mmCIF_output.find_block(outputTail)

        # for item in gemmiBlock:
        #   singular_output_block.add_item(item)


# with open("/home/harold/Dev/privateer_python/project_alliance/single.mmCIF", mode="w") as newfile:
#   singular_mmCIF_output.write_file("/home/harold/Dev/privateer_python/project_alliance/single.mmCIF")

# mmcif = structure.make_mmcif_document()
# mmcif.write_file("Cluster1.mmCIF")
