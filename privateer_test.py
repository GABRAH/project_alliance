from privateer import privateer_core as pvt

# glycosylation = pvt.GlycosylationComposition("/home/harold/Dev/privateer_python/tests/test_data/2h6o_carbremediation.pdb")
glycosylation = pvt.GlycosylationComposition("/home/harold/Dev/privateer_python/project_alliance/glycampdbfiles/VolumeConverted/man8.2/Cluster2.pdb")


listOfDetectedGlycans = glycosylation.get_summary_of_detected_glycans()

print("\n")
for entry in listOfDetectedGlycans:
    for key, value in entry.items():
        print('{}: {}'.format(key, value))
    print("_______________________")

glycan = glycosylation.get_glycan(0)

print("Summary of the glycan: " + str(glycan.get_glycan_summary()))

sugar = glycan.get_monosaccharide(6)
sugarSummary = sugar.get_sugar_summary()

print(sugarSummary)
print("\n")

print("is_sane: " + str(sugar.is_sane()))
print("If one of the following is False, then is_sane will also return False")
print("\tsugar_diag_puckering: " + str(sugar.ok_with_puckering()))
print("\tsugar_diag_anomer: " + str(sugar.ok_with_anomer()))
print("\tsugar_diag_chirality: " + str(sugar.ok_with_chirality()))
print("\tsugar_diag_ring: " + str(sugar.ok_with_ring()))
