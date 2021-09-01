import os
import requests
import json
from privateer import privateer_core as pvt
from privateer import privateer_modelling as pvtmodelling

donorPDBID = "MAN6"
receiverPDBID = "FOLD"
outputPDBID = "t35t"

donorpath = (
    f"/home/harold/Dev/privateer_python/project_alliance/builder_test/{donorPDBID}.pdb"
)
receiverpath = f"/home/harold/Dev/privateer_python/project_alliance/builder_test/{receiverPDBID}.pdb"

outputpath = (
    f"/home/harold/Dev/privateer_python/project_alliance/builder_test/{outputPDBID}.pdb"
)


# path = f"/Users/haroldas/Dev/privateer_python/project_sigma/select_cases/{currentPDBID}.pdb"
print(donorpath)


builder = pvtmodelling.Builder(receiverpath, donorpath, False, True, True, False)
sequence = builder.get_receiving_model_sequence_info()
print(sequence)
builder.graft_glycan_to_receiver(0, 0, 37)
builder.graft_glycan_to_receiver(0, 0, 74)
builder.export_grafted_model(outputpath)
grafted_glycan_summary = builder.get_summary_of_grafted_glycans()

print(grafted_glycan_summary)

# with 5fjj as donor to FOLD model.
# donorPDBID = "5fjj"
# receiverPDBID = "FOLD"
# outputPDBID = "t35t"

# donorpath = (
#     f"/home/harold/Dev/privateer_python/project_alliance/builder_test/{donorPDBID}.pdb"
# )
# receiverpath = f"/home/harold/Dev/privateer_python/project_alliance/builder_test/{receiverPDBID}.pdb"

# outputpath = (
#     f"/home/harold/Dev/privateer_python/project_alliance/builder_test/{outputPDBID}.pdb"
# )

# # path = f"/Users/haroldas/Dev/privateer_python/project_sigma/select_cases/{currentPDBID}.pdb"
# print(donorpath)

# donor = pvt.GlycosylationComposition(donorpath)
# # glycanSummary = donor.get_summary_of_detected_glycans()
# # print(glycanSummary)
# builder = pvtmodelling.Builder(receiverpath, donorpath, True, True, True)
# sequence = builder.get_receiving_model_sequence_info()
# # print(sequence)
# builder.graft_glycan_to_receiver(0, 0, 37)
# builder.export_grafted_model(outputpath)
