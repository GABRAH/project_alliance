import os
import re
import shutil

# import gemmi
# from privateer import libprivateer as pvt

# gemmi will be used to produce reworked PDB files into mmCIF format. Most likely it will also be used to add LINK information into the final outputs.
# Privateer will be used to validate linkage positions(atomic distances) and generate WURCS strings and cross-validate against GlyConnect entries.

"""
  "aB": { "PDB": "ARB",
          "supported": True,
          "full-name": "beta-L-arabinopyranose", 
          "abbreviation": "Ara" },

aB: a = L-arabinose; B = beta anomer
PDB: ARB = PDB three letter code(contains both residue name and anomer info in a single code)
supported: Whether the specific sugar is present in internal Privateer database(clipper-glyco_data.cpp) used for geometric validations.
full-name: full chemical name used in Google and https://www.ebi.ac.uk/pdbe-srv/pdbechem/chemicalCompound/show/LFR search queries for PDB codes.
abbreviation: Glycam/common abbreviation of a sugar unit in mass spec terminology(cuz they don't care about anomeric info lol)

"""
# Input glycam-PDB structure will conver second and third characters into equivalent PDB codes.
glycamOneLetterToPDBThreeLetterCodeConversion = {
    "AA": {
        "PDB": "64K",
        "supported": False,
        "full-name": "alpha-D-arabinopyranose",
        "abbreviation": "Ara",
    },
    "AB": {
        "PDB": "SEJ",
        "supported": False,
        "full-name": "beta-D-arabinopyranose",
        "abbreviation": "Ara",
    },
    "aA": {
        "PDB": "ARA",
        "supported": True,
        "full-name": "alpha-L-arabinopyranose",
        "abbreviation": "Ara",
    },
    "aB": {
        "PDB": "ARB",
        "supported": True,
        "full-name": "beta-L-arabinopyranose",
        "abbreviation": "Ara",
    },
    "AD": {
        "PDB": "BXY",
        "supported": True,
        "full-name": "alpha-D-arabinofuranose",
        "abbreviation": "Ara",
    },
    "AU": {
        "PDB": "BXX",
        "supported": True,
        "full-name": "beta-D-arabinofuranose",
        "abbreviation": "Ara",
    },
    "aD": {
        "PDB": "AHR",
        "supported": True,
        "full-name": "alpha-L-arabinofuranose",
        "abbreviation": "Ara",
    },
    "aU": {
        "PDB": "FUB",
        "supported": False,
        "full-name": "beta-L-arabinofuranose",
        "abbreviation": "Ara",
    },
    "DA": {
        "PDB": "LDY",
        "supported": True,
        "full-name": "alpha-D-lyxopyranose",
        "abbreviation": "Lyx",
    },
    "DB": {
        "PDB": "Z4W",
        "supported": False,
        "full-name": "beta-D-lyxopyranose",
        "abbreviation": "Lyx",
    },
    "dA": {
        "PDB": "UNK",
        "supported": False,
        "full-name": "alpha-L-lyxopyranose",
        "abbreviation": "Lyx",
    },
    "dB": {
        "PDB": "UNK",
        "supported": False,
        "full-name": "beta-L-lyxopyranose",
        "abbreviation": "Lyx",
    },
    "RA": {
        "PDB": "YYM",
        "supported": False,
        "full-name": "alpha-D-ribopyranose",
        "abbreviation": "Rib",
    },
    "RB": {
        "PDB": "RIP",
        "supported": True,
        "full-name": "beta-D-ribopyranose",
        "abbreviation": "Rib",
    },
    "rA": {
        "PDB": "UNK",
        "supported": False,
        "full-name": "alpha-L-ribopyranose",
        "abbreviation": "Rib",
    },
    "rB": {
        "PDB": "0MK",
        "supported": False,
        "full-name": "beta-L-ribopyranose",
        "abbreviation": "Rib",
    },
    "XA": {
        "PDB": "XYS",
        "supported": True,
        "full-name": "alpha-D-xylopyranose",
        "abbreviation": "Xyl",
    },
    "XB": {
        "PDB": "XYP",
        "supported": True,
        "full-name": "beta-D-xylopyranose",
        "abbreviation": "Xyl",
    },
    "xA": {
        "PDB": "HSY",
        "supported": False,
        "full-name": "alpha-L-xylopyranose",
        "abbreviation": "Xyl",
    },
    "xB": {
        "PDB": "LXC",
        "supported": True,
        "full-name": "beta-L-xylopyranose",
        "abbreviation": "Xyl",
    },
    "NA": {
        "PDB": "AFD",
        "supported": False,
        "full-name": "alpha-D-allopyranose",
        "abbreviation": "All",
    },
    "NB": {
        "PDB": "ALL",
        "supported": False,
        "full-name": "beta-D-allopyranose",
        "abbreviation": "All",
    },
    "nA": {
        "PDB": "UNK",
        "supported": False,
        "full-name": "alpha-L-allopyranose",
        "abbreviation": "All",
    },
    "nB": {
        "PDB": "WOO",
        "supported": False,
        "full-name": "beta-L-allopyranose",
        "abbreviation": "All",
    },
    "EA": {
        "PDB": "SHD",
        "supported": False,
        "full-name": "alpha-D-altropyranose",
        "abbreviation": "Alt",
    },
    "EB": {
        "PDB": "UNK",
        "supported": False,
        "full-name": "beta-D-altropyranose",
        "abbreviation": "Alt",
    },
    "eA": {
        "PDB": "UNK",
        "supported": False,
        "full-name": "alpha-L-altropyranose",
        "abbreviation": "Alt",
    },
    "eB": {
        "PDB": "UNK",
        "supported": False,
        "full-name": "beta-L-altropyranose",
        "abbreviation": "Alt",
    },
    "LA": {
        "PDB": "GLA",
        "supported": True,
        "full-name": "alpha-D-galactopyranose",
        "abbreviation": "Gal",
    },
    "LB": {
        "PDB": "GAL",
        "supported": True,
        "full-name": "beta-D-galactopyranose",
        "abbreviation": "Gal",
    },
    "lA": {
        "PDB": "GXL",
        "supported": True,
        "full-name": "alpha-L-galactopyranose",
        "abbreviation": "Gal",
    },
    "lB": {
        "PDB": "GIV",
        "supported": False,
        "full-name": "beta-L-galactopyranose",
        "abbreviation": "Gal",
    },
    "GA": {
        "PDB": "GLC",
        "supported": True,
        "full-name": "alpha-D-glucopyranose",
        "abbreviation": "Glc",
    },
    "GB": {
        "PDB": "BGC",
        "supported": True,
        "full-name": "beta-D-glucopyranose",
        "abbreviation": "Glc",
    },
    "gA": {
        "PDB": "UNK",
        "supported": False,
        "full-name": "alpha-L-glucopyranose",
        "abbreviation": "Glc",
    },
    "gB": {
        "PDB": "Z8T",
        "supported": False,
        "full-name": "beta-L-glucopyranose",
        "abbreviation": "Glc",
    },
    "KA": {
        "PDB": "4GL",
        "supported": False,
        "full-name": "alpha-D-gulopyranose",
        "abbreviation": "Gul",
    },
    "KB": {
        "PDB": "GL0",
        "supported": True,
        "full-name": "beta-D-gulopyranose",
        "abbreviation": "Gul",
    },
    "kA": {
        "PDB": "GUP",
        "supported": True,
        "full-name": "alpha-L-gulopyranose",
        "abbreviation": "Gul",
    },
    "kB": {
        "PDB": "Z8T",
        "supported": False,
        "full-name": "beta-L-gulopyranose",
        "abbreviation": "Gul",
    },
    "IA": {
        "PDB": "ZCD",
        "supported": False,
        "full-name": "alpha-D-idopyranose",
        "abbreviation": "Ido",
    },
    "IB": {
        "PDB": "UNK",
        "supported": False,
        "full-name": "beta-D-idopyranose",
        "abbreviation": "Ido",
    },
    "iA": {
        "PDB": "UNK",
        "supported": False,
        "full-name": "alpha-L-idopyranose",
        "abbreviation": "Ido",
    },
    "iB": {
        "PDB": "4N2",
        "supported": False,
        "full-name": "beta-L-idopyranose",
        "abbreviation": "Ido",
    },
    "MA": {
        "PDB": "MAN",
        "supported": True,
        "full-name": "alpha-D-mannopyranose",
        "abbreviation": "Man",
    },
    "MB": {
        "PDB": "BMA",
        "supported": True,
        "full-name": "beta-D-mannopyranose",
        "abbreviation": "Man",
    },
    "mA": {
        "PDB": "UNK",
        "supported": False,
        "full-name": "alpha-L-mannopyranose",
        "abbreviation": "Man",
    },
    "mB": {
        "PDB": "UNK",
        "supported": False,
        "full-name": "beta-L-mannopyranose",
        "abbreviation": "Man",
    },
    "TA": {
        "PDB": "UNK",
        "supported": False,
        "full-name": "alpha-D-talopyranose",
        "abbreviation": "Tal",
    },
    "TB": {
        "PDB": "SDY",
        "supported": False,
        "full-name": "beta-D-talopyranose",
        "abbreviation": "Tal",
    },
    "tA": {
        "PDB": "UNK",
        "supported": False,
        "full-name": "alpha-L-talopyranose",
        "abbreviation": "Tal",
    },
    "tB": {
        "PDB": "ZEE",
        "supported": False,
        "full-name": "beta-L-talopyranose",
        "abbreviation": "Tal",
    },
    "CA": {
        "PDB": "UNK",
        "supported": False,
        "full-name": "alpha-D-fructopyranose",
        "abbreviation": "Fru",
    },
    "CB": {
        "PDB": "BDF",
        "supported": True,
        "full-name": "beta-D-fructopyranose",
        "abbreviation": "Fru",
    },
    "cA": {
        "PDB": "UNK",
        "supported": False,
        "full-name": "alpha-L-fructopyranose",
        "abbreviation": "Fru",
    },
    "cB": {
        "PDB": "UNK",
        "supported": False,
        "full-name": "beta-L-fructopyranose",
        "abbreviation": "Fru",
    },
    "CD": {
        "PDB": "Z9N",
        "supported": False,
        "full-name": "alpha-D-fructofuranose",
        "abbreviation": "Fru",
    },
    "CU": {
        "PDB": "FRU",
        "supported": True,
        "full-name": "beta-D-fructofuranose",
        "abbreviation": "Fru",
    },
    "cD": {
        "PDB": "UNK",
        "supported": False,
        "full-name": "alpha-L-fructofuranose",
        "abbreviation": "Fru",
    },
    "cU": {
        "PDB": "LFR",
        "supported": True,
        "full-name": "beta-L-fructofuranose",
        "abbreviation": "Fru",
    },
    "PA": {
        "PDB": "UNK",
        "supported": False,
        "full-name": "alpha-D-psicopyranose",
        "abbreviation": "Psi",
    },
    "PB": {
        "PDB": "UNK",
        "supported": False,
        "full-name": "beta-D-psicopyranose",
        "abbreviation": "Psi",
    },
    "pA": {
        "PDB": "UNK",
        "supported": False,
        "full-name": "alpha-L-psicopyranose",
        "abbreviation": "Psi",
    },
    "pB": {
        "PDB": "UNK",
        "supported": False,
        "full-name": "beta-L-psicopyranose",
        "abbreviation": "Psi",
    },
    "PD": {
        "PDB": "PSV",
        "supported": True,
        "full-name": "alpha-D-psicofuranose",
        "abbreviation": "Psi",
    },
    "PU": {
        "PDB": "TTV",
        "supported": False,
        "full-name": "beta-D-psicofuranose",
        "abbreviation": "Psi",
    },
    "pD": {
        "PDB": "SF6",
        "supported": False,
        "full-name": "alpha-L-psicofuranose",
        "abbreviation": "Psi",
    },
    "pU": {
        "PDB": "SF9",
        "supported": False,
        "full-name": "beta-L-psicofuranose",
        "abbreviation": "Psi",
    },
    "BA": {
        "PDB": "UNK",
        "supported": False,
        "full-name": "alpha-D-sorbopyranose",
        "abbreviation": "Sor",
    },
    "BB": {
        "PDB": "UNK",
        "supported": False,
        "full-name": "beta-D-sorbopyranose",
        "abbreviation": "Sor",
    },
    "bA": {
        "PDB": "SOE",
        "supported": True,
        "full-name": "alpha-L-sorbopyranose",
        "abbreviation": "Sor",
    },
    "bB": {
        "PDB": "UNK",
        "supported": False,
        "full-name": "beta-L-sorbopyranose",
        "abbreviation": "Sor",
    },
    "JA": {
        "PDB": "T6T",
        "supported": False,
        "full-name": "alpha-D-tagatopyranose",
        "abbreviation": "Tag",
    },
    "JB": {
        "PDB": "UNK",
        "supported": False,
        "full-name": "beta-D-tagatopyranose",
        "abbreviation": "Tag",
    },
    "jA": {
        "PDB": "UNK",
        "supported": False,
        "full-name": "alpha-L-tagatopyranose",
        "abbreviation": "Tag",
    },
    "jB": {
        "PDB": "UNK",
        "supported": False,
        "full-name": "beta-L-tagatopyranose",
        "abbreviation": "Tag",
    },
    "FA": {
        "PDB": "FCA",
        "supported": True,
        "full-name": "alpha-D-fucopyranose",
        "abbreviation": "Fuc",
    },
    "FB": {
        "PDB": "FCB",
        "supported": True,
        "full-name": "beta-D-fucopyranose",
        "abbreviation": "Fuc",
    },
    "fA": {
        "PDB": "FUC",
        "supported": True,
        "full-name": "alpha-L-fucopyranose",
        "abbreviation": "Fuc",
    },
    "fB": {
        "PDB": "FUL",
        "supported": True,
        "full-name": "beta-L-fucopyranose",
        "abbreviation": "Fuc",
    },
    "QA": {
        "PDB": "G6D",
        "supported": False,
        "full-name": "alpha-D-quinovopyranose",
        "abbreviation": "Qui",
    },
    "QB": {
        "PDB": "YYK",
        "supported": False,
        "full-name": "beta-D-quinovopyranose",
        "abbreviation": "Qui",
    },
    "qA": {
        "PDB": "UNK",
        "supported": False,
        "full-name": "alpha-L-quinovopyranose",
        "abbreviation": "Qui",
    },
    "qB": {
        "PDB": "UNK",
        "supported": False,
        "full-name": "beta-L-quinovopyranose",
        "abbreviation": "Qui",
    },
    "HA": {
        "PDB": "XXR",
        "supported": False,
        "full-name": "alpha-D-rhamnopyranose",
        "abbreviation": "Rha",
    },
    "HB": {
        "PDB": "UNK",
        "supported": False,
        "full-name": "beta-D-rhamnopyranose",
        "abbreviation": "Rha",
    },
    "hA": {
        "PDB": "RAM",
        "supported": True,
        "full-name": "alpha-L-rhamnopyranose",
        "abbreviation": "Rha",
    },
    "hB": {
        "PDB": "RM4",
        "supported": True,
        "full-name": "beta-L-rhamnopyranose",
        "abbreviation": "Rha",
    },
    "OA": {
        "PDB": "ADA",
        "supported": True,
        "full-name": "alpha-D-galactopyranuronic acid",
        "abbreviation": "GalA",
    },
    "OB": {
        "PDB": "GTR",
        "supported": True,
        "full-name": "beta-D-galactopyranuronic acid",
        "abbreviation": "GalA",
    },
    "oA": {
        "PDB": "UNK",
        "supported": False,
        "full-name": "alpha-L-galactopyranuronic acid",
        "abbreviation": "GalA",
    },
    "oB": {
        "PDB": "UNK",
        "supported": False,
        "full-name": "beta-L-galactopyranuronic acid",
        "abbreviation": "GalA",
    },
    "ZA": {
        "PDB": "GCU",
        "supported": True,
        "full-name": "alpha-D-glucopyranuronic acid",
        "abbreviation": "GlcA",
    },
    "ZB": {
        "PDB": "BDP",
        "supported": True,
        "full-name": "beta-D-glucopyranuronic acid",
        "abbreviation": "GlcA",
    },
    "zA": {
        "PDB": "UNK",
        "supported": False,
        "full-name": "alpha-L-glucopyranuronic acid",
        "abbreviation": "GlcA",
    },
    "zB": {
        "PDB": "UNK",
        "supported": False,
        "full-name": "beta-L-glucopyranuronic acid",
        "abbreviation": "GlcA",
    },
    "UA": {
        "PDB": "UNK",
        "supported": False,
        "full-name": "alpha-D-idopyranuronic acid",
        "abbreviation": "IdoA",
    },
    "UB": {
        "PDB": "UNK",
        "supported": False,
        "full-name": "beta-D-idopyranuronic acid",
        "abbreviation": "IdoA",
    },
    "uA": {
        "PDB": "IDR",
        "supported": True,
        "full-name": "alpha-L-idopyranuronic acid",
        "abbreviation": "IdoA",
    },
    "uB": {
        "PDB": "UNK",
        "supported": False,
        "full-name": "beta-L-idopyranuronic acid",
        "abbreviation": "IdoA",
    },
    "VA": {
        "PDB": "A2G",
        "supported": True,
        "full-name": "2-acetamido-2-deoxy-alpha-D-galactopyranose",
        "abbreviation": "GalNAc",
    },
    "VB": {
        "PDB": "NGA",
        "supported": True,
        "full-name": "2-acetamido-2-deoxy-beta-D-galactopyranose",
        "abbreviation": "GalNAc",
    },
    "vA": {
        "PDB": "YYQ",
        "supported": False,
        "full-name": "2-acetamido-2-deoxy-alpha-L-galactopyranose",
        "abbreviation": "GalNAc",
    },
    "vB": {
        "PDB": "UNK",
        "supported": False,
        "full-name": "2-acetamido-2-deoxy-beta-L-galactopyranose",
        "abbreviation": "GalNAc",
    },
    "YA": {
        "PDB": "NDG",
        "supported": True,
        "full-name": "2-acetamido-2-deoxy-alpha-D-glucoopyranose",
        "abbreviation": "GlcNAc",
    },
    "YB": {
        "PDB": "NAG",
        "supported": True,
        "full-name": "2-acetamido-2-deoxy-beta-D-glucopyranose",
        "abbreviation": "GlcNAc",
    },
    "yA": {
        "PDB": "NGZ",
        "supported": True,
        "full-name": "2-acetamido-2-deoxy-alpha-L-glucoopyranose",
        "abbreviation": "GlcNAc",
    },
    "yB": {
        "PDB": "UNK",
        "supported": False,
        "full-name": "2-acetamido-2-deoxy-beta-L-glucoopyranose",
        "abbreviation": "GlcNAc",
    },
    "WA": {
        "PDB": "BM3",
        "supported": True,
        "full-name": "2-acetamido-2-deoxy-alpha-D-mannopyranose",
        "abbreviation": "ManNAc",
    },
    "WB": {
        "PDB": "BM7",
        "supported": True,
        "full-name": "2-acetamido-2-deoxy-beta-D-mannopyranose",
        "abbreviation": "ManNAc",
    },
    "wA": {
        "PDB": "UNK",
        "supported": False,
        "full-name": "2-acetamido-2-deoxy-alpha-L-mannoopyranose",
        "abbreviation": "ManNAc",
    },
    "wB": {
        "PDB": "UNK",
        "supported": False,
        "full-name": "2-acetamido-2-deoxy-beta-L-mannoopyranose",
        "abbreviation": "ManNAc",
    },
    "SA": {
        "PDB": "SIA",
        "supported": True,
        "full-name": "5-N-ACETYL-ALPHA-D-NEURAMINIC ACID",
        "abbreviation": "NeuNAc" or "Neu5Ac",
    },
    "SB": {
        "PDB": "SLB",
        "supported": True,
        "full-name": "5-N-ACETYL-BETA-D-NEURAMINIC ACID",
        "abbreviation": "NeuNAc" or "Neu5Ac",
    },
    "KNA": {
        "PDB": "KDM",
        "supported": True,
        "full-name": "3-deoxy-D-glycero-alpha-D-galacto-non-2-ulopyranosonic acid",
        "abbreviation": "KDN",
    },
    "KNB": {
        "PDB": "KDN",
        "supported": True,
        "full-name": "3-deoxy-D-glycero-beta-D-galacto-non-2-ulopyranosonic acid",
        "abbreviation": "KDN",
    },
    "KOA"
    or "KO"
    or "KOB": {
        "PDB": "KDO",
        "supported": True,
        "full-name": "3-deoxy-alpha-D-manno-oct-2-ulopyranosonic acid",
        "abbreviation": "KDO",
    },
    "SGA": {
        "PDB": "NGC",
        "supported": False,
        "full-name": "N-glycolyl-alpha-neuraminic acid",
        "abbreviation": "NeuNGc" or "Neu5Gc",
    },
    "SGB": {
        "PDB": "NGE",
        "supported": False,
        "full-name": "N-glycolyl-beta-neuraminic acid",
        "abbreviation": "NeuNGc" or "Neu5Gc",
    },
}


def import_pdb(path):
    outputLines = []
    file = open(path, "r")
    Lines = file.readlines()
    file.close()

    for line in Lines:
        if type(line) is not str:
            decodedLine = line.decode("utf-8")
            if decodedLine[:5] != "MODEL":
                outputLines.append(decodedLine)
        else:
            if line[:5] != "MODEL":
                outputLines.append(line)
    return outputLines


def export_pdb(path, output):
    file = open(path, "a")
    file.writelines(output)
    file.close()


def CreateFolder(path):
    if not os.path.exists(path):
        os.makedirs(path)
    else:
        shutil.rmtree(path)  # Removes all the subdirectories!
        os.makedirs(path)


def replace_char_at_index(org_str, index, replacement):
    """Replace character at index in string org_str with the
    given replacement character."""
    new_str = org_str
    if index < len(org_str):
        new_str = org_str[0:index] + replacement + org_str[index + 1 :]
    return new_str


def find_location_of_glycam_residue_code(line):
    regex = "|".join(
        "([\w]({}))".format(k) for k in glycamOneLetterToPDBThreeLetterCodeConversion
    )
    match = re.search(regex, line)
    if match is not None:
        return {"start": match.start(), "end": match.end(), "match": match.group()}
    else:
        return None


def generateROHReplacementInstructions(Lines):
    for idx, line in enumerate(Lines):
        ROHregex = "(ROH)"
        ROHmatch = re.search(ROHregex, line)
        if ROHmatch is not None:
            while idx < len(Lines):
                glycamCodeMatch = find_location_of_glycam_residue_code(Lines[idx])
                if glycamCodeMatch is not None:
                    splitLine = Lines[idx].split()
                    residueID = splitLine[4]
                    return {
                        "replacement_code": glycamCodeMatch["match"],
                        "replacement_residue_ID": residueID,
                    }
                idx += 1
    return None


def replaceROH(instructionsForROHReplacement, Lines):
    output = []
    if instructionsForROHReplacement == None:
        return Lines
    else:
        replacementResidue = instructionsForROHReplacement["replacement_code"]
        replacement_residue_ID = instructionsForROHReplacement["replacement_residue_ID"]
        for line in Lines:
            if "ROH" in line:
                ROHreplaced = line.replace("ROH", replacementResidue)
                splitLine = re.split(r"(\s+)", ROHreplaced)
                splitLine[8] = re.sub("\d", replacement_residue_ID, splitLine[8])
                outputLine = "".join(splitLine)
                output.append(outputLine)
            else:
                output.append(line)
    return output


def convertGlycamToPDB(glycamOneLetterToPDBThreeLetterCodeConversion, inputPDB, path):
    output = []
    unsupportedByPrivaterCodes = []

    regex = "|".join(
        "([\w]({}))".format(k) for k in glycamOneLetterToPDBThreeLetterCodeConversion
    )

    for line in inputPDB:
        splitLine = re.split(r"(\s+)", line)
        if len(splitLine) > 3:
            match = re.search(regex, splitLine[6])
            if match is not None:
                codeToReplace = splitLine[6][match.start() + 1 : match.end()]
                replacementDictionary = glycamOneLetterToPDBThreeLetterCodeConversion[
                    codeToReplace
                ]
                replacementPDBCode = replacementDictionary["PDB"]
                supportedByPrivateer = replacementDictionary["supported"]
                splitLine[6] = re.sub(match.group(), replacementPDBCode, splitLine[6])
                outputLine = "".join(splitLine)
                output.append(outputLine)

                if supportedByPrivateer == False:
                    unsupportedByPrivaterCodes.append(
                        {"GlycamCode": match.group(), "PDBCode": replacementPDBCode}
                    )
        else:
            output.append(line)

    if len(unsupportedByPrivaterCodes):
        print(
            f"Detected {len(unsupportedByPrivaterCodes)} sugars that are not supported by Privateer in the following structure {path}!"
        )
        for idx, unsupportedCode in enumerate(unsupportedByPrivaterCodes):
            print(
                f'{idx}/{len(unsupportedByPrivaterCodes)}: Glycam ID - {unsupportedCode["GlycamCode"]}\t\tPDB ID: {unsupportedCode["PDBCode"]}'
            )

    return output


def conversionPipeline(path):
    glycamPDB = import_pdb(path)
    instructionsForROHReplacement = generateROHReplacementInstructions(glycamPDB)
    ROH_removed = replaceROH(instructionsForROHReplacement, glycamPDB)
    convertedPDB = convertGlycamToPDB(
        glycamOneLetterToPDBThreeLetterCodeConversion, ROH_removed, path
    )
    return convertedPDB


# inputpath = (
#     "/home/harold/Dev/privateer_python/project_alliance/glycampdbfiles/omannose/"
# )
# outputpath = "/home/harold/Dev/privateer_python/project_alliance/glycampdbfiles/omannoseConvertedPDB/"

inputpath = (
    "/home/harold/Dev/privateer_python/project_alliance/glycampdbfiles/omannose/"
)
outputpath = "/home/harold/Dev/privateer_python/project_alliance/glycampdbfiles/omannoseConvertedPDB/"
CreateFolder(outputpath)

for root, dirs, files in os.walk(inputpath, topdown=False):
    for name in files:
        head, tail = os.path.split(root)
        outputroot = os.path.join(outputpath, tail)
        if not os.path.exists(outputroot):
            os.makedirs(outputroot)
        with open(os.path.join(outputroot, name), mode="w") as newfile:
            newfile.writelines(conversionPipeline(os.path.join(root, name)))