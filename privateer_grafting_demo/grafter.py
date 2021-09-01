import os
import re
import sys
import json
import argparse
import requests
from privateer import privateer_core as pvtcore
from privateer import privateer_modelling as pvtmodelling


# regex for N glycosylation - [N][^P][ST]|[N][A-Z][C]
donorPDBID = "MAN6"
receiverPDBID = "FOLD"
outputPDBID = "t35t"

donorpath = "/home/harold/Dev/privateer_python/project_alliance/privateer_grafting_demo/input/glycanblocks/man9/cluster1.pdb"
# donorpath = (
#     f"/home/harold/Dev/privateer_python/project_alliance/builder_test/{donorPDBID}.pdb"
# )

receiverpath = f"/home/harold/Dev/privateer_python/project_alliance/privateer_grafting_demo/input/receiving_model/P29016.pdb"
# receiverpath = f"/home/harold/Dev/privateer_python/project_alliance/builder_test/{receiverPDBID}.pdb"

outputpath = "/home/harold/Dev/privateer_python/project_alliance/privateer_grafting_demo/output/P29016.pdb"
# outputpath = (
#     f"/home/harold/Dev/privateer_python/project_alliance/builder_test/{outputPDBID}.pdb"
# )

uniprotID = "P29016"


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


def local_input_model_pipeline(receiverpath, donorpath, outputpath, uniprotID):
    sequences = get_sequences_in_receiving_model(receiverpath)
    if uniprotID is not None:
        uniprotQuery = query_uniprot_for_glycosylation_locations(uniprotID)
        # function to align uniprot sequence with input model
        # then do the grafting
    else:
        targets = get_NGlycosylation_targets_via_consensus_seq(sequences)
        graftedGlycans = glycosylate_receiving_model_using_consensus_seq(
            receiverpath, donorpath, outputpath, targets, True, False
        )

        for idx, graft in enumerate(graftedGlycans):
            index = graft["index"]
            proteinChainID = graft["receiving_protein_residue_chain_PDBID"]
            proteinPDBID = graft["receiving_protein_residue_monomer_PDBID"]
            proteinResidueType = graft["receiving_protein_residue_monomer_type"]
            graftedGlycanChainID = graft["glycan_grafted_as_chainID"]

            if len(graft["ClashingResidues"]):
                averageTotalAtomicDistance = graft["AvgTotalAtomicDistance"]
                numberOfClashingResidues = len(graft["ClashingResidues"])
                print(
                    f"{idx+1}/{len(graftedGlycans)}: Grafted donor glycan as chain {graftedGlycanChainID} to {proteinChainID}/{proteinResidueType}-{proteinPDBID}. The graft has resulted in {numberOfClashingResidues} clashes with an average atomic distance of: {averageTotalAtomicDistance}."
                )
            else:
                print(
                    f"{idx+1}/{len(graftedGlycans)}: Grafted donor glycan as chain {graftedGlycanChainID} to {proteinChainID}/{proteinResidueType}-{proteinPDBID}. The graft did not produce any clashes."
                )


local_input_model_pipeline(receiverpath, donorpath, outputpath, uniprotID)
