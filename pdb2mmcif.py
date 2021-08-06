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


# Privateer/gemmi functions
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

            currentResidueAtomName = "O" + connection["donorPosition"]
            connectedToResidueAtomName = "C" + connection["acceptorPosition"]

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


def convertSinglePDBtoSingleCIF(inputFilePath, outputFilePath, privateerJSON):
    root, fileName = os.path.split(inputFilePath)
    trash, glycanName = os.path.split(root)

    outputName = ""
    if fileName == "Cluster1.pdb":
        outputName = "Phi=66_Psi=-179_Omega=-177"
    elif fileName == "Cluster2.pdb":
        outputName = "Phi=66_Psi=-179_Omega=55.1"

    privateerMetaData = getMetadataFromPrivateer(inputFilePath, privateerJSON)
    gemmiStructure = gemmi.read_structure(inputFilePath)
    gemmiStructureWithSugarLinks = addGemmiConnectionsBetweenSugars(
        gemmiStructure, privateerMetaData
    )

    gemmiDocument = gemmiStructureWithSugarLinks.make_mmcif_document()
    gemmiBlock = gemmiDocument.sole_block()
    gemmiBlock.set_pair("_WURCS", gemmi.cif.quote(privateerMetaData["glycanWURCS"]))
    gemmiBlock.set_pair("_GlyTouCan", gemmi.cif.quote(privateerMetaData["glytoucanID"]))
    gemmiBlock.name = glycanName + "/" + outputName
    gemmiBlock.set_pair("_entry.id", glycanName)
    gemmiBlock.set_pair("_cell.entry_id", glycanName)
    gemmiBlock.set_pair("_symmetry.entry_id", glycanName)
    gemmiDocument.write_file(outputFilePath)


def convertAllPDBtoSingleCIF(
    inputFilePath, singular_mmCIF_output, privateerJSON, blockNameList
):
    root, fileName = os.path.split(inputFilePath)
    trash, glycanName = os.path.split(root)

    outputName = ""
    if fileName == "Cluster1.pdb":
        outputName = "Phi=66_Psi=-179_Omega=-177"
    elif fileName == "Cluster2.pdb":
        outputName = "Phi=66_Psi=-179_Omega=55.1"

    privateerMetaData = getMetadataFromPrivateer(inputFilePath, privateerJSON)
    gemmiStructure = gemmi.read_structure(inputFilePath)
    gemmiStructureWithSugarLinks = addGemmiConnectionsBetweenSugars(
        gemmiStructure, privateerMetaData
    )

    gemmiDocument = gemmiStructureWithSugarLinks.make_mmcif_document()
    gemmiBlock = gemmiDocument.sole_block()
    gemmiBlock.set_pair("_WURCS", gemmi.cif.quote(privateerMetaData["glycanWURCS"]))
    gemmiBlock.set_pair("_GlyTouCan", gemmi.cif.quote(privateerMetaData["glytoucanID"]))

    gemmiBlock.name = glycanName + "/" + outputName
    blockNameList.append(gemmiBlock.name)

    gemmiBlock.set_pair("_entry.id", glycanName)
    gemmiBlock.set_pair("_cell.entry_id", glycanName)
    gemmiBlock.set_pair("_symmetry.entry_id", glycanName)
    singular_mmCIF_output.add_copied_block(gemmiBlock, pos=-1)


inputPath = "/home/harold/Dev/privateer_python/project_alliance/glycampdbfiles/VolumeConvertedPDB/"
# inputPath = "/Users/haroldas/Dev/privateer_python/project_alliance/glycampdbfiles/VolumeConvertedPDB/"
outputPath = "/home/harold/Dev/privateer_python/project_alliance/glycampdbfiles/VolumeConvertedmmCIF/"
# outputPath = "/Users/haroldas/Dev/privateer_python/project_alliance/glycampdbfiles/VolumeConvertedmmCIF/"
privateerJSON = pvt.OfflineDatabase()

single_mmCIF_output_path = os.path.join(outputPath, "compilation.mmCIF")
singular_mmCIF_output = gemmi.cif.Document()
blockNameList = []
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

        convertSinglePDBtoSingleCIF(inputFilePath, outputFilePath, privateerJSON)
        convertAllPDBtoSingleCIF(
            inputFilePath, singular_mmCIF_output, privateerJSON, blockNameList
        )

indexBlock = singular_mmCIF_output.add_new_block(name="index", pos=0)
indexLoop = indexBlock.init_loop("_index.", ["id", "block_name"])
for count, item in enumerate(blockNameList):
    indexLoop.add_row([str(count), item])

singular_mmCIF_output.write_file(single_mmCIF_output_path)
