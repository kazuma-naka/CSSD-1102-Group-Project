import doctest

import requests
from rdkit import Chem
from rdkit.Chem import Draw
import json
from collections import Counter
import csv
from io import StringIO

"""
############################################################

Complete functions are provided to you. Please do not make
modifications.

You can find this code: https://github.com/kazuma-naka/CSSD-1102-Group-Project

############################################################
"""


def make_query(url: str) -> str:
    """
    Given a URL, make a requst to the ChemPub online database
    The interpreter complains if the escape character is not used (this is add only for the purpose of accommodating doctest)

    >>> make_query("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/water/cids/txt")
    '962\\n'
    """
    result = requests.get(url)

    # You can find what the code represents here: https://pubchem.ncbi.nlm.nih.gov//docs/pug-rest#section=Status-Codes
    assert result.status_code == 200, f"Error with url {url}\nstatus_code: {result.status_code}\nTo interpret the status code, please visit https://pubchem.ncbi.nlm.nih.gov//docs/pug-rest#section=Status-Codes\ncontent: {result.content}"

    return result.text


def visualize_multiple_smiles(smiles_lst: list[str], cids: list[int], mol_name: str) -> None:
    """
    Given a list of Canonical smiles string and a list of corresponding cids, visualize the structures

    It will create a png file in the same directory of where the script is stored.
    """
    mols = []
    for smiles in smiles_lst:
        mol = Chem.MolFromSmiles(smiles)
        # Identify potential stereo bonds!
        Chem.FindPotentialStereoBonds(mol)
        mols.append(mol)

    mylegends = [f"CID {cid}" for cid in cids]

    img = Draw.MolsToGridImage(
        mols, molsPerRow=2, subImgSize=(400, 400), legends=mylegends)
    img.save(f'{mol_name}.png')


"""
############################################################

You are welcome to add more help functions as needed.
Please add those under this section

############################################################
"""

DRUG_NAMES = [
    "Nicotine",
    "Enzalutamide",
    "Alpelisib",
    "Crizotinib",
    "Axitinib",
    "Vemurafenib",
    "Vandetanib",
    "Dacomitinib",
    "Avapritinib",
    "Selumetinib",
    "Pexidartinib",
]


def get_non_CH_atoms_from_smiles(smiles: str) -> list[str]:
    """
    Extract unique element symbols except C and H from a SMILES string.

    Returns a sorted list of element symbols (e.g., ["N", "O", "Cl"]).
    """
    atoms = set()
    i = 0

    while i < len(smiles):
        ch = smiles[i]

        if ch.isalpha():
            if i + 1 < len(smiles) and smiles[i + 1].islower():
                symbol = smiles[i:i+2]
                i += 2
            else:
                symbol = ch
                i += 1

            if symbol not in ("C", "H"):
                atoms.add(symbol)
        else:
            i += 1

    return sorted(atoms)


def analyze_drug_for_clo5(mol_name: str, top_k: int = 10) -> None:
    print("=" * 80)
    print(f"Analyzing drug: {mol_name}")

    mol_cid = name_to_cid(mol_name)
    print(f"CID of {mol_name}: {mol_cid}")

    similar_cids = find_similar_structures(mol_cid)
    print(f"Number of similar structures: {len(similar_cids)}")

    similar_cids = remove_query_cid(similar_cids, mol_cid)
    print(
        f"After removing same connectivity, {len(similar_cids)} compounds left")

    similar_cids = sort_by_frequency(similar_cids)
    print("Top 10 CIDs after frequency sorting:")
    for cid in similar_cids[:10]:
        print(f"\t{cid}")

    drug_candidate = remove_non_drug(
        similar_cids, save_result=True, mol_name=mol_name
    )
    print(
        f"After removing non-drug-like structures, {len(drug_candidate)} candidates left."
    )

    if not drug_candidate:
        print("No drug-like candidates found. Skipping property analysis.")
        return

    properties = ["MolecularWeight", "XLogP"]

    rows = sort_and_retrieve_property(
        cids=drug_candidate,
        compound_property="MolecularWeight",
        properties_to_include=properties,
        descending=True,
    )

    if not rows:
        print("No property rows retrieved. Skipping.")
        return

    top_rows = rows[:top_k]

    for row in top_rows:
        cid_str = row.get("CID", "").strip()
        if cid_str:
            try:
                smiles = cid_to_smiles(int(cid_str))
            except Exception:
                smiles = ""
        else:
            smiles = ""
        row["CanonicalSMILES"] = smiles

    print("\nTop 10 SMILES (sorted by MolecularWeight):")
    for row in top_rows:
        smiles = row.get("CanonicalSMILES", "").strip()
        xlogp = row.get("XLogP", "").strip()
        mw = row.get("MolecularWeight", "").strip()
        print(f"{smiles}: XlogP - {xlogp}  (MW={mw})")

    min_xlogp = None
    most_polar_smiles: list[str] = []

    for row in top_rows:
        smiles = row.get("CanonicalSMILES", "").strip()
        xlogp_str = row.get("XLogP", "").strip()
        if not smiles or not xlogp_str:
            continue
        try:
            x = float(xlogp_str)
        except ValueError:
            continue

        if min_xlogp is None or x < min_xlogp:
            min_xlogp = x
            most_polar_smiles = [smiles]
        elif x == min_xlogp:
            most_polar_smiles.append(smiles)

    print("\nMost polar by XLogP (SMILES):")
    for s in most_polar_smiles:
        print(s)

    if most_polar_smiles:
        atoms = get_non_CH_atoms_from_smiles(most_polar_smiles[0])
        print("\nNon-C and Non-H atoms:")
        print(", ".join(atoms))

    cids_to_draw = [int(row["CID"]) for row in top_rows]
    smiles_to_draw = [row["CanonicalSMILES"].strip() for row in top_rows]
    visualize_multiple_smiles(smiles_to_draw, cids_to_draw, mol_name)
    print(f"Generated {mol_name}.png for the top {len(top_rows)} candidates.")


"""
############################################################

Please implement functions in this section. Please follow
the design recipe specified for each function

############################################################
"""


def generate_url(input: str, operation: str, output: str, options: str = "") -> str:
    """
    Generates a URL for the PubChem PUG REST API based on the provided input, operation, output format, and options.

    >>> generate_url("compound/name/butadiene", "property/HeavyAtomCount", "txt")
    'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/butadiene/property/HeavyAtomCount/txt'

    >>> generate_url("compound/cid/129825914,129742624,129783988", "sids", "json", "list_return=grouped")
    'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/129825914,129742624,129783988/sids/json?list_return=grouped'
    """
    base = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/"
    url = f"{base}{input}/{operation}/{output}"
    if options:
        url = f"{url}?{options}"
    return url


def name_to_cid(name: str) -> int:
    """
    Given the name of a compound, return the cid

    >>> name_to_cid("water")
    962
    """
    url = generate_url(f"compound/name/{name}", "cids", "txt")
    text = make_query(url)
    cid_str = text.strip()
    return int(cid_str)


def find_similar_structures(cid: int) -> list[int]:
    """
    Given a valid cid, return a list of similar structures based on fastsimilarity_2d

    >>> find_similar_structures(962)
    [962, 961, 24602, 104752, 22247451, 10129877, 157350, 6914120, 123332, 105142, 137349153, 139859, 5460554, 10197601, 16048639, 6914119, 17985087, 21871715, 22146520, 44144404, 57470146, 129631210, 171528, 10214376, 20266850, 57422081, 72499223, 5460604, 12676391, 21338204, 23560202, 57418915, 57634601, 57937911, 58422532, 58609138, 58750500, 58828612, 59045931, 59049839, 59102850, 59367517, 59984341, 78060557, 85595863, 85602483, 139617572, 157621312, 158956564, 160528969, 161115002, 163575618, 139505, 12161503, 13644079, 15593902, 17769776, 19818716, 20037327, 71327550, 71333220, 72700706, 76071083]
    """
    url = generate_url(
        f"compound/fastsimilarity_2d/cid/{cid}",
        "cids",
        "txt",
    )
    text = make_query(url)
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    return [int(line) for line in lines]


def remove_query_cid(cids: list[int], query_cid: int) -> list[int]:
    """
    Given a list of cids with similar structures to the query_cid based on fastidentity, remove the compound with the same connectivity

    >>> remove_query_cid([962, 961, 24602, 104752, 22247451, 157350, 10129877, 6914120, 123332, 105142, 137349153, 139859, 5460554, 10197601, 16048639, 6914119, 17985087, 21871715, 22146520, 44144404, 57470146, 129631210, 171528, 10214376, 20266850, 57422081, 72499223, 5460604, 12676391, 21338204, 23560202, 57418915, 57634601, 57937911, 58422532, 58609138, 58750500, 58828612, 59045931, 59049839, 59102850, 59367517, 59984341, 78060557, 85595863, 85602483, 139617572, 157621312, 158956564, 160528969, 161115002, 163575618, 139505, 12161503, 13644079, 15593902, 17769776, 19818716, 20037327, 71327550, 71333220, 72700706, 76071083], 962)
    [961, 22247451, 157350, 123332, 137349153, 5460554, 6914119, 17985087, 21871715, 22146520, 44144404, 57470146, 129631210, 171528, 20266850, 57422081, 72499223, 21338204, 23560202, 57418915, 57634601, 57937911, 58422532, 58609138, 58750500, 58828612, 59045931, 59049839, 59102850, 59367517, 59984341, 78060557, 85595863, 85602483, 139617572, 157621312, 158956564, 161115002, 163575618, 13644079, 15593902, 17769776, 19818716, 20037327, 71327550, 71333220]
    """
    url = generate_url(
        f"compound/fastidentity/cid/{query_cid}",
        "cids",
        "json",
        "identity_type=same_connectivity",
    )
    text = make_query(url)
    data = json.loads(text)

    same_connectivity_cids = set(
        data.get("IdentifierList", {}).get("CID", [])
    )
    filtered = [cid for cid in cids if cid not in same_connectivity_cids]

    return filtered


def sort_by_frequency(cids: list[int]) -> list[int]:
    """
    Sort the cids by frequency by descending order and remove duplicate items.

    NOTE: Your implementation may result in a list with different orders as in the example below as you are free to handle ties in any way you want
    In the autograder on PrairieLearn, we will NOT include examples with ties.

    >>> sort_by_frequency([961, 22247451, 157350, 123332, 137349153, 5460554, 6914119, 17985087, 21871715, 22146520, 44144404, 57470146, 129631210, 171528, 20266850, 57422081, 72499223, 21338204, 23560202, 57418915, 57634601, 57937911, 58422532, 58609138, 58750500, 58828612, 59045931, 59049839, 59102850, 59367517, 59984341, 78060557, 85595863, 85602483, 139617572, 157621312, 158956564, 161115002, 163575618, 13644079, 15593902, 17769776, 19818716, 20037327, 71327550, 71333220])
    [57422081, 59102850, 58422532, 171528, 23560202, 78060557, 44144404, 158956564, 72499223, 22247451, 137349153, 57418915, 58750500, 139617572, 157350, 57634601, 59045931, 15593902, 13644079, 17769776, 85602483, 71327550, 17985087, 157621312, 961, 57470146, 163575618, 123332, 58828612, 6914119, 5460554, 20037327, 59984341, 85595863, 22146520, 21338204, 59367517, 19818716, 20266850, 21871715, 71333220, 129631210, 59049839, 58609138, 57937911, 161115002]
    """
    if not cids:
        return []

    counts = Counter(cids)

    first_index: dict[int, int] = {}
    for i, cid in enumerate(cids):
        if cid not in first_index:
            first_index[cid] = i

    unique_cids = list(counts.keys())

    unique_cids.sort(key=lambda cid: (-counts[cid], first_index[cid]))

    return unique_cids


def get_chunks(cids: list[int], chunk_size=100) -> list[list[int]]:
    """
    Divide the list into multiple chunks of certain chunk_size.

    >>> get_chunks([57422081, 59102850, 58422532, 171528, 23560202, 78060557, 44144404, 158956564, 72499223, 22247451, 137349153, 57418915, 58750500, 139617572, 157350, 57634601, 59045931, 15593902, 13644079, 17769776, 85602483, 71327550, 17985087, 157621312, 961, 57470146, 163575618, 123332, 58828612, 6914119, 5460554, 20037327, 59984341, 85595863, 22146520, 21338204, 59367517, 19818716, 20266850, 21871715, 71333220, 129631210, 59049839, 58609138, 57937911, 161115002])
    [[57422081, 59102850, 58422532, 171528, 23560202, 78060557, 44144404, 158956564, 72499223, 22247451, 137349153, 57418915, 58750500, 139617572, 157350, 57634601, 59045931, 15593902, 13644079, 17769776, 85602483, 71327550, 17985087, 157621312, 961, 57470146, 163575618, 123332, 58828612, 6914119, 5460554, 20037327, 59984341, 85595863, 22146520, 21338204, 59367517, 19818716, 20266850, 21871715, 71333220, 129631210, 59049839, 58609138, 57937911, 161115002]]
    """
    if chunk_size <= 0:
        raise ValueError("chunk_size must be a positive integer")

    chunks: list[list[int]] = []

    for i in range(0, len(cids), chunk_size):
        chunks.append(cids[i:i + chunk_size])

    return chunks


def remove_non_drug(sorted_cids: list[int], save_result=False, mol_name="") -> list[int]:
    """
    Given a list of sorted cids, remove the non-drug structures based on Lipinski's rule of five
    The queries involved in this function will be time-consuming, you will need to create mulitple chunks.
    Ignore the entries that don't have all five properties available.

    Args:
         sorted_cids - a list of sorted cids
         save_result - whether to save the result to a csv file. The default is False. Please make sure this is set to the False when send your code to the autograder
         mol_name    - the name of the compound. This will be used as the name of the csv file generated. You must provide it when save_result is True

    >>> remove_non_drug([57422081, 59102850, 58422532, 171528, 23560202, 78060557, 44144404, 158956564, 72499223, 22247451, 137349153, 57418915, 58750500, 139617572, 157350, 57634601, 59045931, 15593902, 13644079, 17769776, 85602483, 71327550, 17985087, 157621312, 961, 57470146, 163575618, 123332, 58828612, 6914119, 5460554, 20037327, 59984341, 85595863, 22146520, 21338204, 59367517, 19818716, 20266850, 21871715, 71333220, 129631210, 59049839, 58609138, 57937911, 161115002])
    [58422532, 171528, 137349153, 58750500, 157350, 57634601, 961, 123332, 6914119, 5460554, 59049839, 58609138, 57937911]
    """
    if not sorted_cids:
        return []

    # Lipinski-like conditions
    MAX_MW = 500.0
    MAX_XLOGP = 5.0
    MAX_HBD = 5
    MAX_HBA = 10

    # This dict will store properties for each CID that passes Lipinski rules
    # so that we can optionally write them to a CSV file.
    druglike_props: dict[int, tuple[float, float, int, int]] = {}

    # Split the CID list into chunks to avoid overly large single requests.
    chunks = get_chunks(sorted_cids, chunk_size=100)

    for chunk in chunks:
        if not chunk:
            continue

        # Build a comma-separated list of CIDs for this chunk
        cid_str = ",".join(str(cid) for cid in chunk)

        # Prepare the property list to query from PubChem
        properties = "MolecularWeight,XLogP,HBondDonorCount,HBondAcceptorCount"

        # Example URL:
        # https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/962,910/property/MolecularWeight,XLogP,HBondDonorCount,HBondAcceptorCount/csv
        url = generate_url(
            f"compound/cid/{cid_str}",
            f"property/{properties}",
            "csv",
        )

        # Query PubChem and get the CSV text
        csv_text = make_query(url)

        # Parse the CSV text using csv.DictReader for convenience
        reader = csv.DictReader(StringIO(csv_text))

        for row in reader:
            # Skip rows that are missing any of the required columns
            try:
                cid_value = row.get("CID", "").strip()
                mw_value = row.get("MolecularWeight", "").strip()
                xlogp_value = row.get("XLogP", "").strip()
                hbd_value = row.get("HBondDonorCount", "").strip()
                hba_value = row.get("HBondAcceptorCount", "").strip()

                # If any field is empty or missing, ignore this entry
                if not (cid_value and mw_value and xlogp_value and hbd_value and hba_value):
                    continue

                cid = int(cid_value)
                mw = float(mw_value)
                xlogp = float(xlogp_value)
                # Sometimes these may be represented as floats in the CSV, so cast via float -> int
                hbd = int(float(hbd_value))
                hba = int(float(hba_value))
            except (ValueError, KeyError):
                # If parsing fails for any reason, skip this row
                continue

            # Apply Lipinski's rule of five
            if (
                mw <= MAX_MW
                and xlogp <= MAX_XLOGP
                and hbd <= MAX_HBD
                and hba <= MAX_HBA
            ):
                druglike_props[cid] = (mw, xlogp, hbd, hba)

    # Preserve the original order from sorted_cids,
    # but only keep CIDs that passed the Lipinski filter.
    druglike_cids = [cid for cid in sorted_cids if cid in druglike_props]

    # Optionally write results to a CSV file
    if save_result:
        if not mol_name:
            # A name is required when saving to file, to determine the file name
            raise ValueError(
                "mol_name must be provided when save_result is True")

        filename = f"{mol_name}.csv"
        with open(filename, "w", newline="") as f:
            writer = csv.writer(f)
            # Header row
            writer.writerow(
                ["CID", "MolecularWeight", "XLogP",
                    "HBondDonorCount", "HBondAcceptorCount"]
            )
            # Write properties in the same order as druglike_cids
            for cid in druglike_cids:
                mw, xlogp, hbd, hba = druglike_props[cid]
                writer.writerow([cid, mw, xlogp, hbd, hba])

    return druglike_cids


def cid_to_smiles(cid: int) -> str:
    """
    Given the cid, return the canonical smiles string

    >>> cid_to_smiles(962)
    'O'
    """

    url = generate_url(

        f"compound/cid/{cid}",
        "property/CanonicalSMILES",
        "txt",
    )
    smiles_text = make_query(url)
    smiles = smiles_text.strip()

    return smiles


def sort_and_retrieve_property(cids: list[int], compound_property: str, properties_to_include: list[str], descending=True):
    """
    Given a list of cids, retrieve a few properties (e.g., properties_to_include) and sort them based on a particular compound_property.
    If decending is Ture, then sort them in the decedning order. Otherwise, sort them in the ascending order.
    The queries involved in this function will be time-consuming, you will need to create mulitple chunks.

    This function is used to answer questions. We will not test the correctness of this function in the autograder.

    You are free to return the results in any format. Or you may choose not to return those, but to print it or save the results to file.
    """

    if not cids or not properties_to_include:
        return []

    properties_str = ",".join(properties_to_include)

    all_rows: list[dict[str, str]] = []

    chunks = get_chunks(cids, chunk_size=100)

    for chunk in chunks:
        if not chunk:
            continue

        cid_str = ",".join(str(cid) for cid in chunk)

        url = generate_url(
            f"compound/cid/{cid_str}",
            f"property/{properties_str}",
            "csv",
        )

        csv_text = make_query(url)

        reader = csv.DictReader(StringIO(csv_text))

        for row in reader:
            cid_value = row.get("CID", "").strip()
            sort_value = row.get(compound_property, "").strip()

            if not cid_value or not sort_value:
                continue

            all_rows.append(row)

    def _sort_key(row: dict[str, str]):
        value_str = row.get(compound_property, "").strip()
        try:
            return float(value_str)
        except ValueError:
            return value_str

    all_rows.sort(key=_sort_key, reverse=descending)

    return all_rows


if __name__ == "__main__":
    for mol_name in DRUG_NAMES:
        analyze_drug_for_clo5(mol_name)

    # run doctest
    # doctest.testmod()
