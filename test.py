import os
import requests
import json
from privateer import privateer_core as pvt
from privateer import privateer_modelling as pvtmodelling

receiverPDBID = "lmao"
donorPDBID = "5fjj"
donorPDBID = "Cluster2"
outputPDBID = "t35t"
# 5tv4, 3hzk, 3v0w
# 3v0w - very good test for linkage detection.
# 5cpw - 2-8 linkages
# 3u75 - 1-2 linkages
# 5oyj - 1-1 linkages

# donorpath = (
#     f"/home/harold/Dev/privateer_python/project_sigma/select_cases/{donorPDBID}.pdb"
# )

donorpath = f"/home/harold/Dev/privateer_python/project_alliance/glycampdbfiles/VolumeConvertedPDB/man5/{donorPDBID}.pdb"
receiverpath = (
    f"/home/harold/Dev/privateer_python/project_sigma/select_cases/{receiverPDBID}.pdb"
)

outputpath = (
    f"/home/harold/Dev/privateer_python/project_alliance/builder_test/{outputPDBID}.pdb"
)

# path = f"/Users/haroldas/Dev/privateer_python/project_sigma/select_cases/{currentPDBID}.pdb"
print(donorpath)


privateerJSON = pvt.OfflineDatabase()

glycosylation = pvt.GlycosylationComposition_memsafe(donorpath)
numGlycans = glycosylation.get_number_of_glycan_chains_detected()
for i in range(numGlycans):
    glycan = glycosylation.get_glycan(i)
    rootInfo = glycan.get_root_info()
    queryResult = glycan.query_offline_database(privateerJSON, False, False)
    glycanWURCS = queryResult["wurcs"]
    glytoucanID = queryResult["glytoucan_id"]
    rootSugarChainID = rootInfo["RootSugarChainID"]
    rootSugarType = glycan.get_monosaccharide(0).get_name_short()
    rootSugarID = glycan.get_monosaccharide(0).get_sugar_pdb_id().replace(" ", "")
    proteinBackboneChainID = rootInfo["ProteinChainID"]
    proteinResidueType = rootInfo["ProteinResidueType"]
    proteinResidueID = rootInfo["ProteinResidueID"]
    print(
        f"{i} - Detected {glycanWURCS} with GlyTouCanID of {glytoucanID}. In protein backbone {proteinBackboneChainID}/{proteinResidueType}-{proteinResidueID} with root sugar {rootSugarChainID}/{rootSugarType}-{rootSugarID}."
    )

builder = pvtmodelling.Builder(receiverpath, donorpath)
sequence = builder.get_receiving_model_sequence_info()
builder.graft_glycan_to_receiver(0, 0, 37)
builder.export_grafted_model(outputpath)
