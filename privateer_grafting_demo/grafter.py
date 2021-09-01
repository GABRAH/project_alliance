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


def get_working_directory_path(scriptfilepath):
    return os.path.dirname(scriptfilepath)


def import_list_of_uniprotIDs_to_glycosylate(inputFilePath):
    pass


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
                if item["description"][0] == "N":
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


def online_input_model_pipeline(uniprotID, donorpath, outputLocation):
    outputFileName = uniprotID + ".pdb"
    outputpath = os.path.join(outputLocation, outputFileName)
    receiverpath = download_and_prepare_alphafoldDB_model(
        uniprotID, defaultInputModelPath
    )
    uniprotGlycosylationQuery = query_uniprot_for_glycosylation_locations(uniprotID)
    uniprotGlycosylations = uniprotGlycosylationQuery["glycosylations"]
    targets = []
    for item in uniprotGlycosylations:
        if item["description"][0] == "N":
            targets.append(int(item["begin"]) - 1)
    graftedGlycans = glycosylate_receiving_model_using_uniprot_info(
        receiverpath, donorpath, outputpath, targets, True, False
    )
    print_grafted_glycans_summary(graftedGlycans)


scriptFilePath = os.path.abspath(__file__)
workingDirectoryPath = get_working_directory_path(scriptFilePath)
defaultDonorPath = os.path.join(workingDirectoryPath, defaultDonorLocation)
defaultInputModelPath = os.path.join(workingDirectoryPath, defaultInputModelLocation)
defaultOutputModelPath = os.path.join(workingDirectoryPath, defaultOutputModelLocation)

localReceiverPath = f"/home/harold/Dev/privateer_python/project_alliance/privateer_grafting_demo/input/receiving_model/O15552.pdb"

uniprotID = "O15552"


# local_input_model_pipeline(localReceiverPath, defaultDonorPath, outputpath, uniprotID)
online_input_model_pipeline(uniprotID, defaultDonorPath, defaultOutputModelPath)

# Implement argparse here to complete demo for Zenodo release.
