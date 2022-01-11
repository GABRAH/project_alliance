import os
import re
import sys
import argparse
import requests
from privateer import privateer_core as pvtcore
from privateer import privateer_modelling as pvtmodelling


defaultDonorLocation = "input/glycanblocks/man5/cluster1.pdb"
defaultInputModelLocation = "input/receiving_model"
defaultOutputModelLocation = "output"

# The plan for tomorrow in terms of improvements:
# 1. User input as a JSON file for specific grafting instructions
# 1.1. for every defined grafting location, define donor
# 2. Create functions that would allow user to input receiver or donor models to extract grafting location indices in some sort of nice output.
# 3. Check how to make it obvious at local run that users for -output argument need to put output file name.
# 4. If enough time, maybe take a look at TRP mannosylation cuz that is special.


def get_working_directory_path(scriptfilepath):
    return os.path.dirname(scriptfilepath)


def import_list_of_uniprotIDs_to_glycosylate(inputFilePath):
    output = []
    with open(inputFilePath) as file:
        lines = file.readlines()
    for line in lines:
        lineSplitCommaList = line.split(",")
        if lineSplitCommaList:
            for split_line in lineSplitCommaList:
                cleanLine = re.sub("\W+", "", split_line)
                output.append(cleanLine)
        else:
            cleanLine = re.sub("\W+", "", line)
            output.append(cleanLine)

    return output


def download_and_prepare_alphafoldDB_model(uniprotID, downloadLocation):
    outputFileName = uniprotID + ".pdb"
    outputFilePath = os.path.join(downloadLocation, outputFileName)
    requestURL = f"https://alphafold.ebi.ac.uk/files/AF-{uniprotID}-F1-model_v1.pdb"
    query = requests.get(requestURL, allow_redirects=True)

    outputLines = []
    downloadedLines = query.iter_lines()
    for line in downloadedLines:
        decodedLine = line.decode("utf-8")
        if decodedLine[:5] != "MODEL":
            outputLines.append(decodedLine)

    with open(outputFilePath, "w") as file:
        file.writelines("%s\n" % l for l in outputLines)

    print(
        f"Successfully downloaded model from AlphaFoldDB with UniProt ID: '{uniprotID}' to {outputFilePath}"
    )
    return outputFilePath


def query_uniprot_for_glycosylation_locations(uniprotID):
    uniprotRequestURL = f"https://www.ebi.ac.uk/proteins/api/proteins/{uniprotID}"
    uniprotResponse = requests.get(
        uniprotRequestURL, headers={"Accept": "application/json"}
    )
    if not uniprotResponse.ok:
        uniprotResponse.raise_for_status()
        sys.exit()
    uniprotResponseJSON = uniprotResponse.json()
    uniprotSequence = uniprotResponseJSON["sequence"]
    uniprotFeatures = uniprotResponseJSON["features"]
    uniprotGlycosylations = []
    for item in uniprotFeatures:
        if item["type"] == "CARBOHYD":
            uniprotGlycosylations.append(item)

    outputSequence = uniprotSequence["sequence"]
    outputSequenceLength = uniprotSequence["length"]
    output = {
        "sequenceLength": outputSequenceLength,
        "sequence": outputSequence,
        "glycosylations": uniprotGlycosylations,
    }

    return output


def get_sequences_in_receiving_model(receiverpath):
    builder_sequence_only = pvtmodelling.Builder(receiverpath, True)
    receiver_sequence = builder_sequence_only.get_receiving_model_sequence_info()

    return receiver_sequence


def get_NGlycosylation_targets_via_consensus_seq(sequences):
    NGlycosylationConsensus = "[N][^P][ST]|[N][A-Z][C]"
    output = []

    for item in sequences:
        currentChainIndex = item["index"]
        currentChainID = item["ChainID"]
        currentSequence = item["Sequence"]

        glycosylationTargets = []

        for match in re.finditer(NGlycosylationConsensus, currentSequence):
            if (
                currentSequence[match.start()]
                == item["Residues"][match.start()]["residueCode"]
            ):
                glycosylationTargets.append(
                    {"start": match.start(), "end": match.end(), "match": match.group()}
                )

        output.append(
            {
                "Sequence": currentSequence,
                "chainIndex": currentChainIndex,
                "currentChainID": currentChainID,
                "glycosylationTargets": glycosylationTargets,
            }
        )

    return output


def glycosylate_receiving_model_using_consensus_seq(
    receiverpath,
    donorpath,
    outputpath,
    glycosylationTargets,
    enableUserMessages,
    trimGlycanIfClashesDetected,
):
    builder = pvtmodelling.Builder(
        receiverpath,
        donorpath,
        -1,
        trimGlycanIfClashesDetected,
        True,
        enableUserMessages,
        False,
    )
    for item in glycosylationTargets:
        chainIndex = item["chainIndex"]
        targets = item["glycosylationTargets"]
        for target in targets:
            currentTargetIndex = target["start"]
            builder.graft_glycan_to_receiver(0, chainIndex, currentTargetIndex)

    graftedGlycanSummary = builder.get_summary_of_grafted_glycans()
    builder.export_grafted_model(outputpath)

    return graftedGlycanSummary


def glycosylate_receiving_model_using_uniprot_info(
    receiverpath,
    donorpath,
    outputpath,
    targets,
    enableUserMessages,
    trimGlycanIfClashesDetected,
):
    builder = pvtmodelling.Builder(
        receiverpath,
        donorpath,
        -1,
        trimGlycanIfClashesDetected,
        True,
        enableUserMessages,
        False,
    )
    for currentTarget in targets:
        chainIndex = 0
        builder.graft_glycan_to_receiver(0, chainIndex, currentTarget)

    graftedGlycanSummary = builder.get_summary_of_grafted_glycans()
    builder.export_grafted_model(outputpath)

    return graftedGlycanSummary


def print_grafted_glycans_summary(graftedGlycans):
    for idx, graft in enumerate(graftedGlycans):
        proteinChainID = graft["receiving_protein_residue_chain_PDBID"]
        proteinPDBID = graft["receiving_protein_residue_monomer_PDBID"]
        proteinResidueType = graft["receiving_protein_residue_monomer_type"]
        graftedGlycanChainID = graft["glycan_grafted_as_chainID"]

        if len(graft["ClashingResidues"]):
            averageTotalAtomicDistance = graft["AvgTotalAtomicDistance"]
            numberOfClashingResidues = len(graft["ClashingResidues"])
            print(
                f"{idx+1}/{len(graftedGlycans)}: Grafted donor glycan as chain {graftedGlycanChainID} to {proteinChainID}/{proteinResidueType}-{proteinPDBID}. The graft has resulted in {numberOfClashingResidues} clashes with an average atomic distance of from detected clashing residues: {averageTotalAtomicDistance}."
            )
        else:
            print(
                f"{idx+1}/{len(graftedGlycans)}: Grafted donor glycan as chain {graftedGlycanChainID} to {proteinChainID}/{proteinResidueType}-{proteinPDBID}. The graft did not produce any clashes."
            )


def local_input_model_pipeline(receiverpath, donorpath, outputpath, uniprotID):
    sequences = get_sequences_in_receiving_model(receiverpath)
    if uniprotID is not None:
        uniprotQuery = query_uniprot_for_glycosylation_locations(uniprotID)
        uniprotSequence = uniprotQuery["sequence"]
        receiverModelSequence = sequences[0]["Sequence"]
        if receiverModelSequence != uniprotSequence:
            raise ValueError(
                "Receiving model sequence does not match the sequence retrieved from UniProt. Please graft glycans using consensus sequence glycan grafting method"
            )
        else:
            uniprotGlycosylations = uniprotQuery["glycosylations"]
            targets = []
            for item in uniprotGlycosylations:
                if (
                    item["description"][0]
                    == "N"
                    # or item["description"][0] == "O"
                    # or item["description"][0] == "S"
                    # or item["description"][0] == "C"
                    # or item["description"][0] == "P"
                ):
                    targets.append(int(item["begin"]) - 1)
            graftedGlycans = glycosylate_receiving_model_using_uniprot_info(
                receiverpath, donorpath, outputpath, targets, True, False
            )
            print_grafted_glycans_summary(graftedGlycans)

    else:
        targets = get_NGlycosylation_targets_via_consensus_seq(sequences)
        graftedGlycans = glycosylate_receiving_model_using_consensus_seq(
            receiverpath, donorpath, outputpath, targets, True, False
        )
        print_grafted_glycans_summary(graftedGlycans)


# P27918 is a good test for TRP, as TRP is a bit more unique and requires different approach. C-mannosylation currently no worky.
def online_input_model_pipeline(
    uniprotID, donorpath, defaultInputModelPath, outputLocation
):
    outputFileName = uniprotID + ".pdb"
    outputpath = os.path.join(outputLocation, outputFileName)
    receiverpath = download_and_prepare_alphafoldDB_model(
        uniprotID, defaultInputModelPath
    )
    uniprotGlycosylationQuery = query_uniprot_for_glycosylation_locations(uniprotID)
    uniprotGlycosylations = uniprotGlycosylationQuery["glycosylations"]
    targets = []
    for item in uniprotGlycosylations:
        if (
            item["description"][0]
            == "N"
            # or item["description"][0] == "O"
            # or item["description"][0] == "S"
            # or item["description"][0] == "C"
            # or item["description"][0] == "P"
        ):
            targets.append(int(item["begin"]) - 1)
    graftedGlycans = glycosylate_receiving_model_using_uniprot_info(
        receiverpath, donorpath, outputpath, targets, True, False
    )
    print_grafted_glycans_summary(graftedGlycans)


scriptFilePath = os.path.abspath(__file__)
workingDirectoryPath = get_working_directory_path(scriptFilePath)
defaultDonorPath = os.path.join(workingDirectoryPath, defaultDonorLocation)
defaultInputModelDirectory = os.path.join(
    workingDirectoryPath, defaultInputModelLocation
)
defaultOutputModelPath = os.path.join(workingDirectoryPath, defaultOutputModelLocation)
defaultuniprotIDsListPath = os.path.join(scriptFilePath, "uniprotIDinputs.txt")

defaultUniprotID = "P29016"


parser = argparse.ArgumentParser(
    prog="grafter.py",
    usage="%(prog)s [options]. Most convenient usage: python grafter.py -import_uniprotIDs_from_file uniprotIDinputs.txt",
    description=f"Graft Glycans to AlphaFoldDB models using Privateer Modelling module.",
    epilog=f"If -local_receiver_path or -uniprotID are not provided, the script will default to using UniProtID: {defaultUniprotID} as default input. Will download the PDB from AlphaFoldDB and N-glycosylate according to UniProt data.",
)
parser.add_argument(
    "-uniprotID",
    action="store",
    default=None,
    dest="user_uniprotID",
    help="If used with -local_receiver_path, N-glycosylate according to UniProt targets. If used without -local_receiver_path, this variable is used in the download of AlphaFoldDB .pdb file and N-Glycosylation according to UniProt targets.",
)
parser.add_argument(
    "-local_receiver_path",
    action="store",
    default=None,
    dest="user_localReceiverPath",
    help=f"Path to locally saved AlpfaFoldDB model on the computer. If -uniprotID is not provided, will carry out N-glycosylation according to regex consensus sequence of '[N][^P][ST]|[N][A-Z][C]'. The argument overrides default behaviour of downloading AlpfaFoldDB model from the server. WARNING: Ensure that \"MODEL 0\" line is deleted in the local file, as otherwise Privateer's MMBD dependency will not be able to import the model!",
)
parser.add_argument(
    "-donor_path",
    action="store",
    default=None,
    dest="user_donorPath",
    help=f"Path to the glycan that is to be grafted throughout AlphaFoldDB model. If not specified, the script will default to using glycan located in '{defaultDonorPath}'",
)
parser.add_argument(
    "-download_path",
    action="store",
    default=None,
    dest="user_inputModelDirectory",
    help=f"Specify download directory where original AlpfaFoldDB models downloaded from the server should be saved. If unspecified, the script will default to '{defaultInputModelDirectory}'",
)
parser.add_argument(
    "-output_path",
    action="store",
    default=None,
    dest="user_outputPath",
    help=f"Specify output directory where AlpfaFoldDB models with grafted glycans should be saved. If unspecified, the script will default to '{defaultOutputModelPath}'",
)
parser.add_argument(
    "-import_uniprotIDs_from_file",
    action="store",
    default=None,
    dest="user_uniprotIDsList",
    help=f"Glycosylate multiple AlphaFoldDB models from a list of UniProtIDs. Example file is located in '{defaultuniprotIDsListPath}' By default will download files from the server and save them localy in specified or default directory locations.",
)

args = parser.parse_args()


if args.user_uniprotID is not None:
    uniprotID = args.user_uniprotID
else:
    uniprotID = defaultUniprotID
if args.user_donorPath is not None:
    donorPath = args.user_donorPath
else:
    donorPath = defaultDonorPath
if args.user_outputPath is not None:
    outputPath = args.user_outputPath
else:
    outputPath = defaultOutputModelPath
if args.user_inputModelDirectory is not None:
    inputModelDirectory = args.user_inputModelDirectory
else:
    inputModelDirectory = defaultInputModelDirectory

if args.user_uniprotIDsList is not None:
    uniprotIDListPath = args.user_uniprotIDsList


if args.user_localReceiverPath is not None and args.user_uniprotID is None:
    uniprotID = None
    local_input_model_pipeline(
        args.user_localReceiverPath, donorPath, outputPath, uniprotID
    )
elif args.user_localReceiverPath is not None and args.user_uniprotID is not None:
    local_input_model_pipeline(
        args.user_localReceiverPath, donorPath, outputPath, uniprotID
    )
elif args.user_uniprotIDsList is not None:
    uniprotIDList = import_list_of_uniprotIDs_to_glycosylate(uniprotIDListPath)
    for idx, uniprotID in enumerate(uniprotIDList):
        online_input_model_pipeline(
            uniprotID, donorPath, inputModelDirectory, outputPath
        )
        print(
            f"\n{idx+1}/{len(uniprotIDList)}: Successfully finished processing AlphaFoldDB model with UniProt ID of {uniprotID}.\n"
        )
else:
    online_input_model_pipeline(uniprotID, donorPath, inputModelDirectory, outputPath)
