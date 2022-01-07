import os
import argparse
from privateer import privateer_core as pvtcore
from prettytable import PrettyTable


def privateerValidation(path):
    outputList = []
    glycosylation = pvtcore.GlycosylationComposition_memsafe(path)
    numGlycans = glycosylation.get_number_of_glycan_chains_detected()
    for glycanIndex in range(numGlycans):
        glycanList = []
        glycan = glycosylation.get_glycan(glycanIndex)
        glycanWURCS = glycan.get_wurcs_notation()
        numSugars = glycan.get_total_number_of_sugars()
        for sugarIndex in range(numSugars):
            sugar = glycan.get_monosaccharide(sugarIndex)
            sugarPDBCode = sugar.get_name_short()
            sugarPDBChain = sugar.get_sugar_chain_id()
            sugarPDBID = int(sugar.get_sugar_pdb_id())
            sugarCremerPople = sugar.get_cremer_pople_params()
            puckeringAmplitude = sugarCremerPople[0]
            Phi = sugarCremerPople[1]
            Theta = sugarCremerPople[2]
            sugarType = sugar.get_denomination()
            sugarConformation = sugar.get_conformation_name()
            sugarBFac = sugar.get_bfactor()
            sugarCtx = sugar.get_glycosylation_type()
            sugarPrivateerDiagnostic = sugar.get_privateer_diagnostic()
            if sugarConformation != "4c1" or sugarPrivateerDiagnostic != "yes":
                glycanList.append(
                    {
                        "Code": sugarPDBCode,
                        "Chain": sugarPDBChain,
                        "ID": sugarPDBID,
                        "Q": puckeringAmplitude,
                        "Phi": Phi,
                        "Theta": Theta,
                        "detectedType": sugarType,
                        "cnf": sugarConformation,
                        "bfac": sugarBFac,
                        "ctx": sugarCtx,
                        "diagnostic": sugarPrivateerDiagnostic,
                    }
                )
        if len(glycanList):
            sortedList = sorted(glycanList, key=lambda k: k["ID"])
            outputDict = {
                "WURCS": glycanWURCS,
                "totalSugars": numSugars,
                "problematicSugars": sortedList,
            }
            outputList.append(outputDict)

    if len(outputList):
        return {"path": path, "Glycans": outputList}
    else:
        return None


def print_privateer_validation_results(fileResults):
    totalSugars = 0
    totalProblematicSugars = 0
    table = PrettyTable()
    table.field_names = [
        "Sugar",
        "Q",
        "Phi",
        "Theta",
        "Detected type",
        "Cnf",
        "<Bfac>",
        "Ctx",
        "Ok?",
    ]
    path = fileResults["path"]
    glycans = fileResults["Glycans"]
    print(f"File path: {path}")
    for idx, glycan in enumerate(glycans):
        print(f'Glycans - {idx+1}/{len(glycans)}: Conversion WURCS: {glycan["WURCS"]}')
        problematicSugars = glycan["problematicSugars"]
        totalSugars += glycan["totalSugars"]
        totalProblematicSugars += len(problematicSugars)
        for sugar in problematicSugars:
            sugarString = f'{sugar["Code"]}-{sugar["Chain"]}-{sugar["ID"]}'
            table.add_row(
                [
                    sugarString,
                    sugar["Q"],
                    sugar["Phi"],
                    sugar["Theta"],
                    sugar["detectedType"],
                    sugar["cnf"],
                    sugar["bfac"],
                    sugar["ctx"],
                    sugar["diagnostic"],
                ]
            )
    print(table)
    print(
        f"{totalProblematicSugars}/{totalSugars} sugars in the input file have been detected as having problems.\n\n"
    )


parser = argparse.ArgumentParser(
    prog="privateer_quick_validate.py",
    usage="%(prog)s PATH.",
    description=f"Quickly validate converted files through Privateer.",
)
required = parser.add_argument_group("required arguments")
required.add_argument(
    "-input",
    action="store",
    dest="user_inputPath",
    help="Input path either directly to a single file or a root directory that contains multiple PDB files(no trailing slash should be left)",
    required=True,
)

args = parser.parse_args()
inputpath = args.user_inputPath

currentDirectory = os.getcwd()
completeInputPath = os.path.join(currentDirectory, inputpath)

if completeInputPath[-1] == "/":
    basePath = os.path.dirname(os.path.dirname(completeInputPath))
else:
    basePath = os.path.dirname(completeInputPath)

if os.path.isdir(completeInputPath):
    inputDirectory = os.path.basename(os.path.normpath(completeInputPath))
    for root, dirs, files in os.walk(completeInputPath, topdown=False):
        for name in files:
            head, tail = os.path.split(root)
            fileResults = privateerValidation(os.path.join(root, name))
            if fileResults is not None:
                print_privateer_validation_results(fileResults)
            else:
                print(f"Privateer detected no issues in: {os.path.join(root, name)}")
elif os.path.isdir(completeInputPath) is False:
    fileResults = privateerValidation(completeInputPath)
    if fileResults is not None:
        print_privateer_validation_results(fileResults)
    else:
        print(f"Privateer detected no issues in: {completeInputPath}")
else:
    raise ValueError(
        f"Unable to determine whether {completeInputPath} is a file or a directory!"
    )
