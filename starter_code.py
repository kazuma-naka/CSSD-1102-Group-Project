import doctest

import requests
from rdkit import Chem
from rdkit.Chem import Draw

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
    pass


def name_to_cid(name: str) -> int:
    """
    Given the name of a compound, return the cid

    >>> name_to_cid("water")
    962
    """
    pass


def find_similar_structures(cid: int) -> list[int]:
    """
    Given a valid cid, return a list of similar structures based on fastsimilarity_2d

    >>> find_similar_structures(962)
    [962, 961, 24602, 104752, 22247451, 10129877, 157350, 6914120, 123332, 105142, 137349153, 139859, 5460554, 10197601, 16048639, 6914119, 17985087, 21871715, 22146520, 44144404, 57470146, 129631210, 171528, 10214376, 20266850, 57422081, 72499223, 5460604, 12676391, 21338204, 23560202, 57418915, 57634601, 57937911, 58422532, 58609138, 58750500, 58828612, 59045931, 59049839, 59102850, 59367517, 59984341, 78060557, 85595863, 85602483, 139617572, 157621312, 158956564, 160528969, 161115002, 163575618, 139505, 12161503, 13644079, 15593902, 17769776, 19818716, 20037327, 71327550, 71333220, 72700706, 76071083]
    """
    pass


def remove_query_cid(cids: list[int], query_cid: int) -> list[int]:
    """
    Given a list of cids with similar structures to the query_cid based on fastidentity, remove the compound with the same connectivity

    >>> remove_query_cid([962, 961, 24602, 104752, 22247451, 157350, 10129877, 6914120, 123332, 105142, 137349153, 139859, 5460554, 10197601, 16048639, 6914119, 17985087, 21871715, 22146520, 44144404, 57470146, 129631210, 171528, 10214376, 20266850, 57422081, 72499223, 5460604, 12676391, 21338204, 23560202, 57418915, 57634601, 57937911, 58422532, 58609138, 58750500, 58828612, 59045931, 59049839, 59102850, 59367517, 59984341, 78060557, 85595863, 85602483, 139617572, 157621312, 158956564, 160528969, 161115002, 163575618, 139505, 12161503, 13644079, 15593902, 17769776, 19818716, 20037327, 71327550, 71333220, 72700706, 76071083], 962)
    [961, 22247451, 157350, 123332, 137349153, 5460554, 6914119, 17985087, 21871715, 22146520, 44144404, 57470146, 129631210, 171528, 20266850, 57422081, 72499223, 21338204, 23560202, 57418915, 57634601, 57937911, 58422532, 58609138, 58750500, 58828612, 59045931, 59049839, 59102850, 59367517, 59984341, 78060557, 85595863, 85602483, 139617572, 157621312, 158956564, 161115002, 163575618, 13644079, 15593902, 17769776, 19818716, 20037327, 71327550, 71333220]
    """
    pass


def sort_by_frequency(cids: list[int]) -> list[int]:
    """
    Sort the cids by frequency by descending order and remove duplicate items.

    NOTE: Your implementation may result in a list with different orders as in the example below as you are free to handle ties in any way you want
    In the autograder on PrairieLearn, we will NOT include examples with ties.

    >>> sort_by_frequency([961, 22247451, 157350, 123332, 137349153, 5460554, 6914119, 17985087, 21871715, 22146520, 44144404, 57470146, 129631210, 171528, 20266850, 57422081, 72499223, 21338204, 23560202, 57418915, 57634601, 57937911, 58422532, 58609138, 58750500, 58828612, 59045931, 59049839, 59102850, 59367517, 59984341, 78060557, 85595863, 85602483, 139617572, 157621312, 158956564, 161115002, 163575618, 13644079, 15593902, 17769776, 19818716, 20037327, 71327550, 71333220])
    [57422081, 59102850, 58422532, 171528, 23560202, 78060557, 44144404, 158956564, 72499223, 22247451, 137349153, 57418915, 58750500, 139617572, 157350, 57634601, 59045931, 15593902, 13644079, 17769776, 85602483, 71327550, 17985087, 157621312, 961, 57470146, 163575618, 123332, 58828612, 6914119, 5460554, 20037327, 59984341, 85595863, 22146520, 21338204, 59367517, 19818716, 20266850, 21871715, 71333220, 129631210, 59049839, 58609138, 57937911, 161115002]
    """
    pass


def get_chunks(cids: list[int], chunk_size=100) -> list[list[int]]:
    """
    Divide the list into multiple chunks of certain chunk_size.

    >>> get_chunks([57422081, 59102850, 58422532, 171528, 23560202, 78060557, 44144404, 158956564, 72499223, 22247451, 137349153, 57418915, 58750500, 139617572, 157350, 57634601, 59045931, 15593902, 13644079, 17769776, 85602483, 71327550, 17985087, 157621312, 961, 57470146, 163575618, 123332, 58828612, 6914119, 5460554, 20037327, 59984341, 85595863, 22146520, 21338204, 59367517, 19818716, 20266850, 21871715, 71333220, 129631210, 59049839, 58609138, 57937911, 161115002])
    [[57422081, 59102850, 58422532, 171528, 23560202, 78060557, 44144404, 158956564, 72499223, 22247451, 137349153, 57418915, 58750500, 139617572, 157350, 57634601, 59045931, 15593902, 13644079, 17769776, 85602483, 71327550, 17985087, 157621312, 961, 57470146, 163575618, 123332, 58828612, 6914119, 5460554, 20037327, 59984341, 85595863, 22146520, 21338204, 59367517, 19818716, 20266850, 21871715, 71333220, 129631210, 59049839, 58609138, 57937911, 161115002]]
    """
    pass


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
    pass


def cid_to_smiles(cid: int) -> str:
    """
    Given the cid, return the canonical smiles string

    >>> cid_to_smiles(962)
    'O'
    """
    pass


def sort_and_retrieve_property(cids: list[int], compound_property: str, properties_to_include: list[str], descending=True):
    """
    Given a list of cids, retrieve a few properties (e.g., properties_to_include) and sort them based on a particular compound_property.
    If decending is Ture, then sort them in the decedning order. Otherwise, sort them in the ascending order.
    The queries involved in this function will be time-consuming, you will need to create mulitple chunks.

    This function is used to answer questions. We will not test the correctness of this function in the autograder.

    You are free to return the results in any format. Or you may choose not to return those, but to print it or save the results to file.
    """
    pass


if __name__ == "__main__":
    mol_name = "water"

    # Retrieve CID
    mol_cid = name_to_cid(mol_name)
    print(f"The CID of the compound {mol_name} is {mol_cid}")

    # Identify compounds with similar chemical structures.
    similar_cids = find_similar_structures(mol_cid)
    print(f"There are {len(similar_cids)} similar structures")

    # Filter out compounds with identical connectivity to remove duplicates.
    similar_cids = remove_query_cid(similar_cids, mol_cid)
    print(
        f"After removing compounds with identical connectivity, {len(similar_cids)} compounds left")

    # Sort the remaining compounds by frequency of occurrences.
    similar_cids = sort_by_frequency(similar_cids)
    compounds = '\n'.join([f'\t{cid}' for cid in similar_cids[:10]])
    print(f"After sorting by frequency, the top 10 compounds are\n{compounds}")

    # Exclude non-drug-like structures based on Lipinski’s rule of five and generate reports
    drug_candidate = remove_non_drug(
        similar_cids, save_result=True, mol_name=mol_name)
    print(
        f"After removing non-drug-like structures, {len(drug_candidate)} compounds left.")
    compounds = '\n'.join([f'\t{cid}' for cid in drug_candidate[:10]])
    print(f"The top 10 compounds are\n{compounds}")

    # Visualize the top ten candidates
    cids_to_draw = drug_candidate[:10]
    smiles_to_draw = [cid_to_smiles(cid) for cid in cids_to_draw]
    visualize_multiple_smiles(smiles_to_draw, cids_to_draw, mol_name)

    # Run these functions to answer questions (you need to uncomment them to run them):
    # sort_and_retrieve_property(drug_candidate, "MolecularWeight", ["MolecularWeight", "XLogP"]) # Q1

    # run doctest
    # doctest.testmod()
